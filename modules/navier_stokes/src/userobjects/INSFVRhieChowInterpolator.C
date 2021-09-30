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
  params.addRequiredParam<VariableName>(NS::pressure, "The pressure variable");
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
    _p(getFunctor<ADReal>(NS::pressure)),
    _p_num(UserObject::_subproblem.getVariable(0, getParam<VariableName>(NS::pressure)).number()),
    _example(0),
    _has_rz(false)
{
  _var_numbers.push_back(_u.number());
  if (_v)
    _var_numbers.push_back(_v->number());
  if (_w)
    _var_numbers.push_back(_w->number());

  const auto & sub_ids = _moose_mesh.meshSubdomains();
  for (const auto sub_id : sub_ids)
  {
    const auto coord_type = _fe_problem.getCoordSystem(sub_id);
    switch (coord_type)
    {
      case Moose::COORD_RZ:
        _has_rz = true;
        break;

      case Moose::COORD_RSPHERICAL:
        mooseError("We don't yet support r-spherical for INSFV");
        break;

      default:
        break;
    }
  }

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
INSFVRhieChowInterpolator::interpolateB(const VectorValue<ADReal> & b_elem,
                                        const VectorValue<ADReal> & b_neighbor,
                                        const FaceInfo & fi) const
{
  if (!fi.neighborPtr())
    return b_elem;

  Real coord_elem;
  coordTransformFactor(
      UserObject::_subproblem, fi.elem().subdomain_id(), fi.elemCentroid(), coord_elem);
  Real coord_neighbor;
  coordTransformFactor(
      UserObject::_subproblem, fi.neighbor().subdomain_id(), fi.neighborCentroid(), coord_neighbor);

  const auto elem_coord_weighting_factor = coord_elem / (coord_elem + coord_neighbor);
  const auto elem_distance_weighting_factor = fi.gC();
  const auto elem_weighting_factor =
      (elem_coord_weighting_factor + elem_distance_weighting_factor) / 2;
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

    const auto elem_id = fi.elem().id();
    // if there are no forces then keys will not already exist so we take advantage of operator[]
    // zero instertion here
    const auto & b_elem = _b[elem_id];
    const auto & b_neighbor = fi.neighborPtr() ? _b[fi.neighborPtr()->id()] : b_elem;
    const auto it = _b1.emplace(std::make_pair(&fi, interpolateB(b_elem, b_neighbor, fi))).first;

    Real coord;
    coordTransformFactor(
        UserObject::_subproblem, fi.elem().subdomain_id(), fi.faceCentroid(), coord);
    const Point surface_vector = fi.normal() * fi.faceArea() * coord;
    auto product = (it->second * fi.dCF()) * surface_vector;

    // Now should we compute _b2 for this element?
    if (dof_map.is_evaluable(fi.elem(), _p_num))
    {
      coordTransformFactor(
          UserObject::_subproblem, fi.elem().subdomain_id(), fi.elemCentroid(), coord);
      // Face info volume just uses libMesh::Elem::volume which has no knowledge of the coordinate
      // system
      const auto elem_volume = coord * fi.elemVolume();
      // Second term in RHS of Mercinger equation 42
      _b2[elem_id] += product * fi.gC() / elem_volume;
      // First term in RHS of Mercinger equation 42
      _b2[elem_id] += surface_vector * _p(&fi.elem()) / elem_volume;
    }

    // Or for the neighbor?
    if (fi.neighborPtr() && dof_map.is_evaluable(fi.neighbor(), _p_num))
    {
      coordTransformFactor(
          UserObject::_subproblem, fi.neighborPtr()->subdomain_id(), fi.neighborCentroid(), coord);
      // Face info volume just uses libMesh::Elem::volume which has no knowledge of the coordinate
      // system
      const auto neighbor_volume = coord * fi.neighborVolume();
      // Second term in RHS of Mercinger equation 42. Apply both a minus sign to the surface vector
      // and to dCF such that result is a + so we don't have to change the sign
      _b2[fi.neighborPtr()->id()] += std::move(product) * (1. - fi.gC()) / neighbor_volume;
      // First term in RHS of Mercinger equation 42. Apply a minus sign to the surface vector
      _b2[fi.neighborPtr()->id()] += -surface_vector * _p(fi.neighborPtr()) / neighbor_volume;
    }
  }

  // We now no longer need to store _b so we can drop its memory
  _b.clear();

  if (!_has_rz)
    return;

  const bool displaced = &(UserObject::_subproblem) != &_fe_problem;
  for (auto * const candidate_elem : _fe_problem.getEvaluableElementRange())
  {
    auto * const elem = displaced ? _moose_mesh.elemPtr(candidate_elem->id()) : candidate_elem;

    const auto coord_system = UserObject::_subproblem.getCoordSystem(elem->subdomain_id());
    if (coord_system == Moose::CoordinateSystemType::COORD_RZ)
    {
      const auto r_coord = UserObject::_subproblem.getAxisymmetricRadialCoord();
      _b2[elem->id()](r_coord) -= _p(elem) / elem->vertex_average()(r_coord);
    }
  }
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

    _b3.emplace(std::make_pair(fi, interpolateB(b_elem, b_neighbor, *fi)));
  }
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
