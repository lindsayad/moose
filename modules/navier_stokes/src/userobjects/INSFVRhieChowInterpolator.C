//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVRhieChowInterpolator.h"
#include "INSFVAttributes.h"
#include "GatherRCDataElementThread.h"
#include "GatherRCDataFaceThread.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "SystemBase.h"

#include "libmesh/mesh_base.h"
#include "libmesh/elem_range.h"
#include "libmesh/parallel_algebra.h"
#include "metaphysicl/dualsemidynamicsparsenumberarray.h"
#include "metaphysicl/parallel_dualnumber.h"
#include "metaphysicl/parallel_dynamic_std_array_wrapper.h"
#include "metaphysicl/parallel_semidynamicsparsenumberarray.h"
#include "timpi/parallel_sync.h"

using namespace libMesh;

registerMooseObject("NavierStokesApp", INSFVRhieChowInterpolator);

InputParameters
INSFVRhieChowInterpolator::validParams()
{
  auto params = GeneralUserObject::validParams();
  params.set<ExecFlagEnum>("execute_on") = {EXEC_NONLINEAR, EXEC_LINEAR};
  params.addRequiredParam<VariableName>("u", "The x-component of velocity");
  params.addParam<VariableName>("v", "The y-component of velocity");
  params.addParam<VariableName>("w", "The z-component of velocity");
  return params;
}

INSFVRhieChowInterpolator::INSFVRhieChowInterpolator(const InputParameters & params)
  : GeneralUserObject(params),
    _moose_mesh(_subproblem.mesh()),
    _mesh(_moose_mesh.getMesh()),
    _u(_subproblem.getVariable(0, getParam<VariableName>("u"))),
    _v(isParamValid("v") ? &_subproblem.getVariable(0, getParam<VariableName>("v")) : nullptr),
    _w(isParamValid("w") ? &_subproblem.getVariable(0, getParam<VariableName>("w")) : nullptr),
    _zero(0)
{
  _var_numbers.insert(_u.number());
  if (_v)
    _var_numbers.insert(_v->number());
  if (_w)
    _var_numbers.insert(_w->number());
}

void
INSFVRhieChowInterpolator::addToA(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  if (elem->processor_id() != this->processor_id())
    _elements_to_push_pull.insert(elem);

  _a[elem->id()](component) += value;
}

void
INSFVRhieChowInterpolator::addToB(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  mooseAssert(elem->processor_id() == this->processor_id(), "Sources should be local");

  // We have our users write their RC data imagining that they've moved all terms to the LHS, but
  // the balance in Moukalled assumes that the body forces are on the RHS with positive sign, e.g.
  // 0 = -\nabla p + \mathbf{B}, so we must apply a minus sign here
  _b[elem->id()](component) -= value;
}

const VectorValue<ADReal> &
INSFVRhieChowInterpolator::rcCoeff(const libMesh::Elem * const elem) const
{
  auto it = _a.find(elem->id());
  mooseAssert(it != _a.end(), "Could not find the requested element with id " << elem->id());
  return it->second;
}

void
INSFVRhieChowInterpolator::initialize()
{
  _elements_to_push_pull.clear();

  _a.clear();
  _b.clear();
  _b1.clear();
  _b2.clear();
  _b3.clear();
}

void
INSFVRhieChowInterpolator::execute()
{
  PARALLEL_TRY
  {
    GatherRCDataElementThread et(_fe_problem, _var_numbers);
    Threads::parallel_reduce(*_moose_mesh.getActiveLocalElementRange(), et);
  }
  PARALLEL_CATCH;

  PARALLEL_TRY
  {
    using FVRange = StoredRange<std::vector<const FaceInfo *>::const_iterator, const FaceInfo *>;
    GatherRCDataFaceThread<FVRange> fvr(_fe_problem, _var_numbers);
    FVRange faces(_fe_problem.mesh().faceInfo().begin(), _fe_problem.mesh().faceInfo().end());
    Threads::parallel_reduce(faces, fvr);
  }
  PARALLEL_CATCH;
}

void
INSFVRhieChowInterpolator::finalizeAData()
{
  using Datum = std::pair<dof_id_type, VectorValue<ADReal>>;
  std::unordered_map<processor_id_type, std::vector<Datum>> push_data;
  std::unordered_map<processor_id_type, std::vector<dof_id_type>> pull_requests;
  static const VectorValue<ADReal> example;

  for (auto * const elem : _elements_to_push_pull)
  {
    const auto id = elem->id();
    const auto pid = elem->processor_id();
    auto it = _a.find(id);
    mooseAssert(it != _a.end(), "We definitely should have found something");
    push_data[pid].push_back(std::make_pair(id, it->second));
    pull_requests[pid].push_back(id);
  }

  // First push
  {
    auto action_functor = [this](const processor_id_type libmesh_dbg_var(pid),
                                 const std::vector<Datum> & sent_data) {
      mooseAssert(pid != this->processor_id(), "We do not send messages to ourself here");
      for (const auto & pr : sent_data)
        _a[pr.first] += pr.second;
    };
    TIMPI::push_parallel_vector_data(_communicator, push_data, action_functor);
  }

  // Then pull
  {
    auto gather_functor = [this](const processor_id_type libmesh_dbg_var(pid),
                                 const std::vector<dof_id_type> & elem_ids,
                                 std::vector<VectorValue<ADReal>> & data_to_fill) {
      mooseAssert(pid != this->processor_id(), "We shouldn't be gathering from ourselves.");
      data_to_fill.resize(elem_ids.size());
      for (const auto i : index_range(elem_ids))
      {
        const auto id = elem_ids[i];
        auto it = _a.find(id);
        mooseAssert(it != _a.end(), "We should hold the value for this locally");
        data_to_fill[i] = it->second;
      }
    };

    auto action_functor = [this](const processor_id_type libmesh_dbg_var(pid),
                                 const std::vector<dof_id_type> & elem_ids,
                                 const std::vector<VectorValue<ADReal>> & filled_data) {
      mooseAssert(pid != this->processor_id(), "The requst filler shouldn't have been ourselves");
      mooseAssert(elem_ids.size() == filled_data.size(), "I think these should be the same size");
      for (const auto i : index_range(elem_ids))
      {
        const auto id = elem_ids[i];
        auto it = _a.find(id);
        mooseAssert(it != _a.end(), "We requested this so we must have it in the map");
        it->second = filled_data[i];
      }
    };
    TIMPI::pull_parallel_vector_data(
        _communicator, pull_requests, gather_functor, action_functor, &example);
  }
}

void
INSFVRhieChowInterpolator::computeFirstAndSecondOverBars()
{
  const auto & all_fi = _moose_mesh.allFaceInfo();
  // reserve to avoid re-allocating all the time
  _b2.reserve(_fe_problem.getEvaluableElementRange().size());

  for (const auto & fi : all_fi)
  {
    if (fi.neighborPtr() == remote_elem)
      // Let's assume that this face information object is at our ghosted elements boundary and in
      // that case we don't care about its face computation
      continue;

    const auto elem_id = fi.elem().id();
    // if there are no forces then keys will not already exist so we take advantage of operator[]
    // zero instertion here
    const auto & b_elem = _b[elem_id];
    const auto & b_neighbor = fi.neighborPtr() ? _b[fi.neighborPtr()->id()] : b_elem;
    const auto it =
        _b1.emplace(
               std::make_pair(&fi, Moose::FV::linearInterpolation(b_elem, b_neighbor, fi, true)))
            .first;

    Real coord;
    coordTransformFactor(_subproblem, fi.elem().subdomain_id(), fi.faceCentroid(), coord);
    const Point surface_vector = fi.normal() * fi.faceArea() * coord;

    // Begin of equation 15.211 in Moukalled. I honestly don't know what to do when we are in an RZ
    // coordinate system since this is supposed to be mimicking the gradient of pressure
    auto product = (it->second * fi.dCF()) * surface_vector;
    coordTransformFactor(_subproblem, fi.elem().subdomain_id(), fi.elemCentroid(), coord);
    _b2[elem_id] += product * fi.gC() / (coord * fi.elemVolume());

    if (fi.neighborPtr())
    {
      coordTransformFactor(
          _subproblem, fi.neighborPtr()->subdomain_id(), fi.neighborCentroid(), coord);
      _b2[fi.neighborPtr()->id()] +=
          std::move(product) * (1. - fi.gC()) / (coord * fi.neighborVolume());
    }
  }

  // We now no longer need to store _b so we can drop its memory
  _b.clear();
}

void
INSFVRhieChowInterpolator::computeThirdOverBar()
{
  const auto & local_fi = _moose_mesh.faceInfo();

  for (auto * const fi : local_fi)
  {
    const auto elem_id = fi->elem().id();
    const auto & b_elem = libmesh_map_find(_b2, elem_id);
    const auto & b_neighbor =
        fi->neighborPtr() ? libmesh_map_find(_b2, fi->neighborPtr()->id()) : b_elem;

    _b3.emplace(std::make_pair(fi, Moose::FV::linearInterpolation(b_elem, b_neighbor, *fi, true)));
  }
}

void
INSFVRhieChowInterpolator::finalizeBData()
{
  // We do not have to push _b data because all that data should initially be
  // local, e.g. we only loop over active local elements for FVElementalKernels

  std::unordered_map<processor_id_type, std::vector<dof_id_type>> pull_requests;

  for (auto * const elem : _fe_problem.getEvaluableElementRange())
  {
    const auto pid = elem->processor_id();
    if (pid != this->processor_id())
      pull_requests[pid].push_back(elem->id());
  }

  auto gather_functor = [this](const processor_id_type libmesh_dbg_var(pid),
                               const std::vector<dof_id_type> & elem_ids,
                               std::vector<VectorValue<ADReal>> & data_to_fill) {
    mooseAssert(pid != this->processor_id(), "We shouldn't be gathering from ourselves.");
    data_to_fill.resize(elem_ids.size());
    for (const auto i : index_range(elem_ids))
    {
      const auto id = elem_ids[i];
      // It's possible that there are no sources in which case we actually want "accidental"
      // insertion of a 0 vector
      data_to_fill[i] = _b[id];
    }
  };

  auto action_functor = [this](const processor_id_type libmesh_dbg_var(pid),
                               const std::vector<dof_id_type> & elem_ids,
                               const std::vector<VectorValue<ADReal>> & filled_data) {
    mooseAssert(pid != this->processor_id(), "The requst filler shouldn't have been ourselves");
    mooseAssert(elem_ids.size() == filled_data.size(), "I think these should be the same size");
    for (const auto i : index_range(elem_ids))
    {
      const auto id = elem_ids[i];
      _b[id] = filled_data[i];
    }
  };
  TIMPI::pull_parallel_vector_data(
      _communicator, pull_requests, gather_functor, action_functor, &_zero);

  // We can proceed to all the overbar operations for _b
  computeFirstAndSecondOverBars();
  computeThirdOverBar();
}

void
INSFVRhieChowInterpolator::finalize()
{
  finalizeAData();
  finalizeBData();
}
