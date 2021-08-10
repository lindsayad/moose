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
#include "timpi/parallel_sync.h"

using namespace libMesh;

registerMooseObject("NavierStokesApp", INSFVRhieChowInterpolator);

namespace TIMPI
{
using MetaPhysicL::DynamicStdArrayWrapper;
using MetaPhysicL::SemiDynamicSparseNumberArray;

template <typename T, typename NType>
class StandardType<DynamicStdArrayWrapper<T, NType>> : public DataType
{
public:
  explicit StandardType(const DynamicStdArrayWrapper<T, NType> * example = nullptr)
  {
    // We need an example for MPI_Address to use
    static const DynamicStdArrayWrapper<T, NType> p{};
    if (!example)
      example = &p;

#ifdef TIMPI_HAVE_MPI

    // Get the sub-data-types, and make sure they live long enough
    // to construct the derived type
    StandardType<std::array<T, NType::size>> d1(&example->_data);
    StandardType<std::size_t> d2(&example->_dynamic_n);

    MPI_Datatype types[] = {(data_type)d1, (data_type)d2};
    int blocklengths[] = {1, 1};
    MPI_Aint displs[2], start;

    timpi_call_mpi(
        MPI_Get_address(const_cast<DynamicStdArrayWrapper<T, NType> *>(example), &start));
    timpi_call_mpi(
        MPI_Get_address(const_cast<std::array<T, NType::size> *>(&example->_data), &displs[0]));
    timpi_call_mpi(MPI_Get_address(const_cast<std::size_t *>(&example->_dynamic_n), &displs[1]));
    displs[0] -= start;
    displs[1] -= start;

    // create a prototype structure
    MPI_Datatype tmptype;
    timpi_call_mpi(MPI_Type_create_struct(2, blocklengths, displs, types, &tmptype));
    timpi_call_mpi(MPI_Type_commit(&tmptype));

    // resize the structure type to account for padding, if any
    timpi_call_mpi(
        MPI_Type_create_resized(tmptype, 0, sizeof(DynamicStdArrayWrapper<T, NType>), &_datatype));
    timpi_call_mpi(MPI_Type_free(&tmptype));

    this->commit();

#endif // TIMPI_HAVE_MPI
  }

  StandardType(const StandardType<DynamicStdArrayWrapper<T, NType>> & timpi_mpi_var(t))
  {
    timpi_call_mpi(MPI_Type_dup(t._datatype, &_datatype));
  }

  ~StandardType() { this->free(); }

  static const bool is_fixed_type = true;
};

template <typename T, typename I, typename N>
class StandardType<SemiDynamicSparseNumberArray<T, I, N>> : public DataType
{
public:
  explicit StandardType(const SemiDynamicSparseNumberArray<T, I, N> * example = nullptr)
  {
    // We need an example for MPI_Address to use
    static const SemiDynamicSparseNumberArray<T, I, N> p;
    if (!example)
      example = &p;

#ifdef TIMPI_HAVE_MPI

    // Get the sub-data-types, and make sure they live long enough
    // to construct the derived type
    StandardType<DynamicStdArrayWrapper<T, N>> d1(&example->nude_data());
    StandardType<DynamicStdArrayWrapper<I, N>> d2(&example->nude_indices());

    MPI_Datatype types[] = {(data_type)d1, (data_type)d2};
    int blocklengths[] = {1, 1};
    MPI_Aint displs[2], start;

    timpi_call_mpi(
        MPI_Get_address(const_cast<SemiDynamicSparseNumberArray<T, I, N> *>(example), &start));
    timpi_call_mpi(MPI_Get_address(
        const_cast<DynamicStdArrayWrapper<T, N> *>(&example->nude_data()), &displs[0]));
    timpi_call_mpi(MPI_Get_address(
        const_cast<DynamicStdArrayWrapper<I, N> *>(&example->nude_indices()), &displs[1]));
    displs[0] -= start;
    displs[1] -= start;

    // create a prototype structure
    MPI_Datatype tmptype;
    timpi_call_mpi(MPI_Type_create_struct(2, blocklengths, displs, types, &tmptype));
    timpi_call_mpi(MPI_Type_commit(&tmptype));

    // resize the structure type to account for padding, if any
    timpi_call_mpi(MPI_Type_create_resized(
        tmptype, 0, sizeof(SemiDynamicSparseNumberArray<T, I, N>), &_datatype));
    timpi_call_mpi(MPI_Type_free(&tmptype));

    this->commit();

#endif // TIMPI_HAVE_MPI
  }

  StandardType(const StandardType<SemiDynamicSparseNumberArray<T, I, N>> & timpi_mpi_var(t))
  {
    timpi_call_mpi(MPI_Type_dup(t._datatype, &_datatype));
  }

  ~StandardType() { this->free(); }

  static const bool is_fixed_type = true;
};
}

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
    _w(isParamValid("w") ? &_subproblem.getVariable(0, getParam<VariableName>("w")) : nullptr)
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

  _b[elem->id()](component) += value;
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
INSFVRhieChowInterpolator::finalize()
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

  {
    // First push data
    auto action_functor = [this](const processor_id_type libmesh_dbg_var(pid),
                                 const std::vector<Datum> & sent_data) {
      mooseAssert(pid != this->processor_id(), "We do not send messages to ourself here");
      for (const auto & pr : sent_data)
        _a[pr.first] += pr.second;
    };
    TIMPI::push_parallel_vector_data(_communicator, push_data, action_functor);
  }

  {
    // Then pull data
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
