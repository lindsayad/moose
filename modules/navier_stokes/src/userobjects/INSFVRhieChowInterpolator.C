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
#include "NS.h"

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
  params += TaggingInterface::validParams();
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_PRE_KERNELS);
  exec_enum = {EXEC_PRE_KERNELS};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  params.addRequiredParam<VariableName>("u", "The x-component of velocity");
  params.addParam<VariableName>("v", "The y-component of velocity");
  params.addParam<VariableName>("w", "The z-component of velocity");
  params.addParam<bool>("standard_body_forces", false, "Whether to just apply normal body forces");
  return params;
}

INSFVRhieChowInterpolator::INSFVRhieChowInterpolator(const InputParameters & params)
  : GeneralUserObject(params),
    TaggingInterface(this),
    _moose_mesh(UserObject::_subproblem.mesh()),
    _mesh(_moose_mesh.getMesh()),
    _sys(*getCheckedPointerParam<SystemBase *>("_sys")),
    _u(UserObject::_subproblem.getVariable(0, getParam<VariableName>("u"))),
    _v(isParamValid("v") ? &UserObject::_subproblem.getVariable(0, getParam<VariableName>("v"))
                         : nullptr),
    _w(isParamValid("w") ? &UserObject::_subproblem.getVariable(0, getParam<VariableName>("w"))
                         : nullptr),
    _example(0),
    _standard_body_forces(getParam<bool>("standard_body_forces"))
{
  _var_numbers.push_back(_u.number());
  if (_v)
    _var_numbers.push_back(_v->number());
  if (_w)
    _var_numbers.push_back(_w->number());

  if (&(UserObject::_subproblem) != &(TaggingInterface::_subproblem))
    mooseError("Different subproblems in INSFVRhieChowInterpolator!");
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

VectorValue<ADReal>
INSFVRhieChowInterpolator::interpolateB(
    std::unordered_map<dof_id_type, VectorValue<ADReal>> & b_container, const FaceInfo & fi)
{
  const auto & b_elem = b_container[fi.elem().id()];
  if (!fi.neighborPtr())
    return b_elem;

  const auto & b_neighbor = b_container[fi.neighbor().id()];
  const auto elem_weighting_factor = fi.gC();
  const auto neighbor_weighting_factor = 1 - elem_weighting_factor;

  return elem_weighting_factor * b_elem + neighbor_weighting_factor * b_neighbor;
}

void
INSFVRhieChowInterpolator::computeFirstAndSecondOverBars()
{
  const auto & dof_map = _sys.dofMap();
  const auto & all_fi = _moose_mesh.allFaceInfo();
  // reserve to avoid re-allocating all the time
  _b2.reserve(_fe_problem.getEvaluableElementRange().size());

  for (const auto & fi : all_fi)
  {
    if (fi.neighborPtr() == remote_elem)
      // Let's assume that this face information object is at our ghosted elements boundary and in
      // that case we don't care about its face computation
      continue;

    const auto it = _b1.emplace(std::make_pair(&fi, interpolateB(_b, fi))).first;

    if (_standard_body_forces)
      continue;

    Real face_coord;
    coordTransformFactor(
        UserObject::_subproblem, fi.elem().subdomain_id(), fi.faceCentroid(), face_coord);
    const Point surface_vector = fi.normal() * fi.faceArea() * face_coord;
    auto product = (it->second * fi.dCF()) * surface_vector;

    // Now should we compute _b2 for this element?
    if (dof_map.is_evaluable(fi.elem(), _var_numbers[0]))
      // Second term in RHS of Mercinger equation 42
      _b2[fi.elem().id()] += product * fi.gC() / _assembly.elementVolume(&fi.elem());

    // Or for the neighbor?
    if (fi.neighborPtr() && dof_map.is_evaluable(fi.neighbor(), _var_numbers[0]))
      // Second term in RHS of Mercinger equation 42. Apply both a minus sign to the surface vector
      // and to dCF such that result is a + so we don't have to change the sign
      _b2[fi.neighbor().id()] +=
          std::move(product) * (1. - fi.gC()) / _assembly.elementVolume(fi.neighborPtr());
  }

  if (_standard_body_forces)
    for (const auto & pr : _b)
      _b2[pr.first] = pr.second;

  // We now no longer need to store _b so we can drop its memory
  _b.clear();
}

void
INSFVRhieChowInterpolator::computeThirdOverBar()
{
  const auto & local_fi = _moose_mesh.faceInfo();

  for (auto * const fi : local_fi)
    _b3.emplace(std::make_pair(fi, interpolateB(_b2, *fi)));
}

void
INSFVRhieChowInterpolator::applyBData()
{
  const auto s = _sys.number();
  for (auto * const elem : *_moose_mesh.getActiveLocalElementRange())
  {
    const auto elem_volume = _assembly.elementVolume(elem);
    for (const auto i : index_range(_var_numbers))
    {
      const auto vn = _var_numbers[i];
      // negative here because we swapped the sign in addToB and so now we need to swap it back
      const auto residual = -elem_volume * libmesh_map_find(_b2, elem->id())(i);
      const auto dof_index = elem->dof_number(s, vn, 0);

      if (_fe_problem.currentlyComputingJacobian())
        _assembly.processDerivatives(residual, dof_index, _matrix_tags);
      else
        _assembly.cacheResidual(dof_index, residual.value(), _vector_tags);
    }
  }

  if (_fe_problem.currentlyComputingJacobian())
    _assembly.addCachedJacobian();
  else
    _assembly.addCachedResiduals();
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
      _communicator, pull_requests, gather_functor, action_functor, &_example);

  // We can proceed to all the overbar operations for _b
  computeFirstAndSecondOverBars();
  computeThirdOverBar();

  // Add the b data to the residual/Jacobian
  applyBData();
}

void
INSFVRhieChowInterpolator::finalize()
{
  finalizeAData();
  finalizeBData();
}
