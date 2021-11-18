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
#include "Reconstructions.h"

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
  params += BlockRestrictable::validParams();
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
    BlockRestrictable(this),
    _moose_mesh(UserObject::_subproblem.mesh()),
    _mesh(_moose_mesh.getMesh()),
    _sys(*getCheckedPointerParam<SystemBase *>("_sys")),
    _u(UserObject::_subproblem.getVariable(0, getParam<VariableName>("u"))),
    _v(isParamValid("v") ? &UserObject::_subproblem.getVariable(0, getParam<VariableName>("v"))
                         : nullptr),
    _w(isParamValid("w") ? &UserObject::_subproblem.getVariable(0, getParam<VariableName>("w"))
                         : nullptr),
    _example(0),
    _standard_body_forces(getParam<bool>("standard_body_forces")),
    _b(_moose_mesh, true),
    _b2(_moose_mesh, true),
    _bx(_b, 0),
    _by(_b, 1),
    _b2x(_b2, 0),
    _b2y(_b2, 1),
    _sub_ids(blockRestricted() ? blockIDs() : _moose_mesh.meshSubdomains())
{
  _var_numbers.push_back(_u.number());
  if (_v)
    _var_numbers.push_back(_v->number());
  if (_w)
    _var_numbers.push_back(_w->number());

  if (&(UserObject::_subproblem) != &(TaggingInterface::_subproblem))
    mooseError("Different subproblems in INSFVRhieChowInterpolator!");

  UserObject::_subproblem.addFunctor("bx", _bx, 0);
  UserObject::_subproblem.addFunctor("by", _by, 0);
  UserObject::_subproblem.addFunctor("b2x", _b2x, 0);
  UserObject::_subproblem.addFunctor("b2y", _b2y, 0);
}

bool
INSFVRhieChowInterpolator::isFaceGeometricallyRelevant(const FaceInfo & fi) const
{
  if (&fi.elem() == libMesh::remote_elem)
    return false;

  bool on_us = _sub_ids.count(fi.elem().subdomain_id());

  if (fi.neighborPtr())
  {
    if (&fi.neighbor() == libMesh::remote_elem)
      return false;

    on_us = on_us || _sub_ids.count(fi.neighbor().subdomain_id());
  }

  return on_us;
}

void
INSFVRhieChowInterpolator::interpolatorSetup()
{
  _elem_range =
      std::make_unique<ConstElemRange>(_mesh.active_local_subdomain_set_elements_begin(_sub_ids),
                                       _mesh.active_local_subdomain_set_elements_end(_sub_ids));

  const auto & all_fi = _moose_mesh.allFaceInfo();
  _evaluable_fi.reserve(all_fi.size());

  const auto & eq = UserObject::_subproblem.es();
  std::vector<const DofMap *> dof_maps(eq.n_systems());
  for (const auto i : make_range(eq.n_systems()))
  {
    const auto & sys = eq.get_system(i);
    dof_maps[i] = &sys.get_dof_map();
  }

  auto is_fi_evaluable = [this, &dof_maps](const FaceInfo & fi) {
    if (!isFaceGeometricallyRelevant(fi))
      return false;

    auto is_elem_evaluable = [&dof_maps](const Elem & elem) {
      for (const auto * const dof_map : dof_maps)
        if (!dof_map->is_evaluable(elem))
          return false;

      return true;
    };

    // We definitely need the face info element to be evaluable in all cases
    if (!is_elem_evaluable(fi.elem()))
      return false;

    if (fi.neighborPtr())
      // We're on an internal face and interpolation to this face will just entail linear
      // interpolation from neighboring cell values. We'll just check then if the neighbor is
      // evaluable (we've already checked the element)
      return is_elem_evaluable(fi.neighbor());

    // Else we are a boundary face. Two-term boundary face extrapolation will require a cell value
    // and gradient, which will require evaluations on all surrounding elements

    mooseAssert(fi.elem().on_boundary(),
                "If we don't have a neighbor pointer, then this element has to be on a boundary");

    for (auto * const neighbor : fi.elem().neighbor_ptr_range())
    {
      if (!neighbor)
        continue;

      if ((neighbor == libMesh::remote_elem) || !is_elem_evaluable(*neighbor))
        return false;
    }

    return true;
  };

  for (const auto & fi : all_fi)
    if (is_fi_evaluable(fi))
      _evaluable_fi.push_back(&fi);

  _evaluable_fi.shrink_to_fit();
}

void
INSFVRhieChowInterpolator::initialSetup()
{
  interpolatorSetup();
}

void
INSFVRhieChowInterpolator::meshChanged()
{
  interpolatorSetup();
}

void
INSFVRhieChowInterpolator::initialize()
{
  _elements_to_push_pull.clear();

  _a.clear();
  _b.clear();
  _b2.clear();
}

void
INSFVRhieChowInterpolator::execute()
{
  PARALLEL_TRY
  {
    GatherRCDataElementThread et(_fe_problem, _var_numbers);
    Threads::parallel_reduce(*_elem_range, et);
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
#ifdef MOOSE_GLOBAL_AD_INDEXING
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
#else
  mooseError("INSFVRhieChowInterpolator only supported for global AD indexing.");
#endif
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
  _b2.reserve(_fe_problem.getEvaluableElementRange().size());

  if (_standard_body_forces)
  {
    for (const auto & pr : _b)
      _b2[pr.first] = pr.second;

    return;
  }

  Moose::FV::reconstruct(_b2, _b, 1, false, false, _evaluable_fi, *this);
}

void
INSFVRhieChowInterpolator::applyBData()
{
  const auto s = _sys.number();
  for (auto * const elem : *_elem_range)
  {
    const auto elem_volume = _assembly.elementVolume(elem);
    for (const auto i : index_range(_var_numbers))
    {
      const auto vn = _var_numbers[i];
      // negative here because we swapped the sign in addToB and so now we need to swap it back
      const auto residual = -elem_volume * _b2[elem->id()](i);
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
#ifdef MOOSE_GLOBAL_AD_INDEXING
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

  // We can proceed to the overbar operations for _b. We only do the first and second overbars. The
  // third overbar is done on the fly when it is requested
  computeFirstAndSecondOverBars();

  // Add the b data to the residual/Jacobian
  applyBData();
#else
  mooseError("INSFVRhieChowInterpolator only supported for global AD indexing.");
#endif
}

void
INSFVRhieChowInterpolator::finalize()
{
  finalizeAData();
  finalizeBData();
}
