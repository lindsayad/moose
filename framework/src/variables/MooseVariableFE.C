//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseVariableFE.h"
#include <typeinfo>
#include "TimeIntegrator.h"
#include "NonlinearSystemBase.h"
#include "DisplacedSystem.h"
#include "Assembly.h"
#include "MooseVariableData.h"

template <>
InputParameters
MooseVariableFE<Real>::validParams()
{
  auto params = MooseVariableField<Real>::validParams();
  params.addClassDescription(
      "Represents standard field variables, e.g. Lagrange, Hermite, or non-constant Monomials");
  return params;
}

template <>
InputParameters
MooseVariableFE<RealVectorValue>::validParams()
{
  auto params = MooseVariableField<RealVectorValue>::validParams();
  params.addClassDescription("Represents vector field variables, e.g. Vector Lagrange or Nedelec");
  return params;
}

template <>
InputParameters
MooseVariableFE<RealEigenVector>::validParams()
{
  auto params = MooseVariableField<RealEigenVector>::validParams();
  params.addClassDescription(
      "Used for grouping standard field variables with the same finite element family and order");
  return params;
}

template <typename RawOutputType>
MooseVariableFE<RawOutputType>::MooseVariableFE(const InputParameters & parameters)
  : MooseVariableField<RawOutputType>(parameters)
{
  _element_data = std::make_unique<MooseVariableData<RawOutputType>>(*this,
                                                                     _sys,
                                                                     _tid,
                                                                     Moose::ElementType::Element,
                                                                     this->_assembly.qRule(),
                                                                     this->_assembly.qRuleFace(),
                                                                     this->_assembly.node(),
                                                                     this->_assembly.elem());
  _neighbor_data = std::make_unique<MooseVariableData<RawOutputType>>(
      *this,
      _sys,
      _tid,
      Moose::ElementType::Neighbor,
      this->_assembly.qRuleNeighbor(), // Place holder
      this->_assembly.qRuleNeighbor(),
      this->_assembly.nodeNeighbor(),
      this->_assembly.neighbor());
  _lower_data = std::make_unique<MooseVariableData<RawOutputType>>(
      *this,
      _sys,
      _tid,
      Moose::ElementType::Lower,
      this->_assembly.qRuleFace(),
      this->_assembly.qRuleFace(), // Place holder
      this->_assembly.node(),      // Place holder
      this->_assembly.lowerDElem());
}

template <typename RawOutputType>
const std::set<SubdomainID> &
MooseVariableFE<RawOutputType>::activeSubdomains() const
{
  return this->_sys.system().variable(_var_num).active_subdomains();
}

template <typename RawOutputType>
Moose::VarFieldType
MooseVariableFE<RawOutputType>::fieldType() const
{
  if (std::is_same<RawOutputType, Real>::value)
    return Moose::VarFieldType::VAR_FIELD_STANDARD;
  else if (std::is_same<RawOutputType, RealVectorValue>::value)
    return Moose::VarFieldType::VAR_FIELD_VECTOR;
  else if (std::is_same<RawOutputType, RealEigenVector>::value)
    return Moose::VarFieldType::VAR_FIELD_ARRAY;
  else
    mooseError("Unknown variable field type");
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::activeOnSubdomain(SubdomainID subdomain) const
{
  return this->_sys.system().variable(_var_num).active_on_subdomain(subdomain);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::clearDofIndices()
{
  _element_data->clearDofIndices();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::prepare()
{
  _element_data->prepare();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::prepareNeighbor()
{
  _neighbor_data->prepare();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::prepareLowerD()
{
  _lower_data->prepare();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::prepareAux()
{
  _element_data->prepareAux();
  _neighbor_data->prepareAux();
  _lower_data->prepareAux();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::reinitNode()
{
  _element_data->reinitNode();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::reinitAux()
{
  _element_data->reinitAux();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::reinitAuxNeighbor()
{
  _neighbor_data->reinitAux();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::reinitNodes(const std::vector<dof_id_type> & nodes)
{
  _element_data->reinitNodes(nodes);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::reinitNodesNeighbor(const std::vector<dof_id_type> & nodes)
{
  _neighbor_data->reinitNodes(nodes);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::getDofIndices(const Elem * elem,
                                              std::vector<dof_id_type> & dof_indices) const
{
  _element_data->getDofIndices(elem, dof_indices);
}

template <typename RawOutputType>
typename MooseVariableFE<RawOutputType>::OutputData
MooseVariableFE<RawOutputType>::getNodalValue(const Node & node) const
{
  return _element_data->getNodalValue(node, Moose::Current);
}

template <typename RawOutputType>
typename MooseVariableFE<RawOutputType>::OutputData
MooseVariableFE<RawOutputType>::getNodalValueOld(const Node & node) const
{
  return _element_data->getNodalValue(node, Moose::Old);
}

template <typename RawOutputType>
typename MooseVariableFE<RawOutputType>::OutputData
MooseVariableFE<RawOutputType>::getNodalValueOlder(const Node & node) const
{
  return _element_data->getNodalValue(node, Moose::Older);
}

template <typename RawOutputType>
typename MooseVariableFE<RawOutputType>::OutputData
MooseVariableFE<RawOutputType>::getElementalValue(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Current, idx);
}

template <typename RawOutputType>
typename MooseVariableFE<RawOutputType>::OutputData
MooseVariableFE<RawOutputType>::getElementalValueOld(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Old, idx);
}

template <typename RawOutputType>
typename MooseVariableFE<RawOutputType>::OutputData
MooseVariableFE<RawOutputType>::getElementalValueOlder(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Older, idx);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::insert(NumericVector<Number> & residual)
{
  _element_data->insert(residual);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::add(NumericVector<Number> & residual)
{
  _element_data->add(residual);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::addSolution(const DenseVector<Number> & v)
{
  _element_data->addSolution(this->_sys.solution(), v);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::addSolutionNeighbor(const DenseVector<Number> & v)
{
  _neighbor_data->addSolution(this->_sys.solution(), v);
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValue() const
{
  mooseDeprecated("Use dofValues instead of dofValue");
  return dofValues();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValues() const
{
  return _element_data->dofValues();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesOld() const
{
  return _element_data->dofValuesOld();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesOlder() const
{
  return _element_data->dofValuesOlder();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesPreviousNL() const
{
  return _element_data->dofValuesPreviousNL();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesNeighbor() const
{
  return _neighbor_data->dofValues();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesOldNeighbor() const
{
  return _neighbor_data->dofValuesOld();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesOlderNeighbor() const
{
  return _neighbor_data->dofValuesOlder();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesPreviousNLNeighbor() const
{
  return _neighbor_data->dofValuesPreviousNL();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDot() const
{
  return _element_data->dofValuesDot();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDotDot() const
{
  return _element_data->dofValuesDotDot();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDotOld() const
{
  return _element_data->dofValuesDotOld();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDotDotOld() const
{
  return _element_data->dofValuesDotDotOld();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDotNeighbor() const
{
  return _neighbor_data->dofValuesDot();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDotDotNeighbor() const
{
  return _neighbor_data->dofValuesDotDot();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDotOldNeighbor() const
{
  return _neighbor_data->dofValuesDotOld();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::dofValuesDotDotOldNeighbor() const
{
  return _neighbor_data->dofValuesDotDotOld();
}

template <typename RawOutputType>
const VariableValue &
MooseVariableFE<RawOutputType>::dofValuesDuDotDu() const
{
  return _element_data->dofValuesDuDotDu();
}

template <typename RawOutputType>
const VariableValue &
MooseVariableFE<RawOutputType>::dofValuesDuDotDotDu() const
{
  return _element_data->dofValuesDuDotDotDu();
}

template <typename RawOutputType>
const VariableValue &
MooseVariableFE<RawOutputType>::dofValuesDuDotDuNeighbor() const
{
  return _neighbor_data->dofValuesDuDotDu();
}

template <typename RawOutputType>
const VariableValue &
MooseVariableFE<RawOutputType>::dofValuesDuDotDotDuNeighbor() const
{
  return _neighbor_data->dofValuesDuDotDotDu();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::prepareIC()
{
  _element_data->prepareIC();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeElemValues()
{
  _element_data->setGeometry(Moose::Volume);
  _element_data->computeValues();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeElemValuesFace()
{
  _element_data->setGeometry(Moose::Face);
  _element_data->computeValues();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeNeighborValuesFace()
{
  _neighbor_data->setGeometry(Moose::Face);
  _neighbor_data->computeValues();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeNeighborValues()
{
  _neighbor_data->setGeometry(Moose::Volume);
  _neighbor_data->computeValues();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeLowerDValues()
{
  _lower_data->setGeometry(Moose::Volume);
  _lower_data->computeValues();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeIncrementAtQps(const NumericVector<Number> & increment_vec)
{
  _element_data->computeIncrementAtQps(increment_vec);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeIncrementAtNode(const NumericVector<Number> & increment_vec)
{
  _element_data->computeIncrementAtNode(increment_vec);
}

template <typename RawOutputType>
RawOutputType
MooseVariableFE<RawOutputType>::getValue(const Elem * elem,
                                         const std::vector<std::vector<OutputShape>> & phi) const
{
  std::vector<dof_id_type> dof_indices;
  this->_dof_map.dof_indices(elem, dof_indices, _var_num);

  OutputType value = 0;
  if (isNodal())
  {
    mooseAssert(dof_indices.size() == phi.size(),
                "The number of shapes does not match the number of dof indices on the elem");

    for (unsigned int i = 0; i < dof_indices.size(); ++i)
    {
      // The zero index is because we only have one point that the phis are evaluated at
      value += phi[i][0] * (*this->_sys.currentSolution())(dof_indices[i]);
    }
  }
  else
  {
    mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
    value = (*this->_sys.currentSolution())(dof_indices[0]);
  }

  return value;
}

template <>
RealEigenVector
MooseVariableFE<RealEigenVector>::getValue(const Elem * elem,
                                           const std::vector<std::vector<Real>> & phi) const
{
  std::vector<dof_id_type> dof_indices;
  this->_dof_map.dof_indices(elem, dof_indices, _var_num);

  RealEigenVector value(_count);
  if (isNodal())
  {
    for (unsigned int i = 0; i < dof_indices.size(); ++i)
      for (unsigned int j = 0; j < _count; j++)
      {
        // The zero index is because we only have one point that the phis are evaluated at
        value(j) += phi[i][0] * (*this->_sys.currentSolution())(dof_indices[i] + j);
      }
  }
  else
  {
    mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
    unsigned int n = 0;
    for (unsigned int j = 0; j < _count; j++)
    {
      value(j) = (*this->_sys.currentSolution())(dof_indices[0] + n);
      n += this->_dof_indices.size();
    }
  }

  return value;
}

template <typename RawOutputType>
typename OutputTools<RawOutputType>::OutputGradient
MooseVariableFE<RawOutputType>::getGradient(
    const Elem * elem, const std::vector<std::vector<OutputShapeGradient>> & grad_phi) const
{
  std::vector<dof_id_type> dof_indices;
  this->_dof_map.dof_indices(elem, dof_indices, _var_num);

  OutputGradient value;
  if (isNodal())
  {
    for (unsigned int i = 0; i < dof_indices.size(); ++i)
    {
      // The zero index is because we only have one point that the phis are evaluated at
      value += grad_phi[i][0] * (*this->_sys.currentSolution())(dof_indices[i]);
    }
  }
  else
  {
    mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
    value = 0.0;
  }

  return value;
}

template <>
RealVectorArrayValue
MooseVariableFE<RealEigenVector>::getGradient(
    const Elem * elem, const std::vector<std::vector<RealVectorValue>> & grad_phi) const
{
  std::vector<dof_id_type> dof_indices;
  this->_dof_map.dof_indices(elem, dof_indices, _var_num);

  RealVectorArrayValue value(_count, LIBMESH_DIM);
  if (isNodal())
  {
    for (unsigned int i = 0; i < dof_indices.size(); ++i)
      for (unsigned int j = 0; j < _count; ++j)
        for (const auto k : make_range(Moose::dim))
        {
          // The zero index is because we only have one point that the phis are evaluated at
          value(j, k) += grad_phi[i][0](k) * (*this->_sys.currentSolution())(dof_indices[i] + j);
        }
  }
  else
  {
    mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
  }

  return value;
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValue() const
{
  return _element_data->nodalValue(Moose::Current);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueNeighbor() const
{
  return _neighbor_data->nodalValue(Moose::Current);
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::nodalVectorTagValue(TagID tag) const
{
  return _element_data->nodalVectorTagValue(tag);
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::DoFValue &
MooseVariableFE<RawOutputType>::nodalMatrixTagValue(TagID tag) const
{
  return _element_data->nodalMatrixTagValue(tag);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueOld() const
{
  return _element_data->nodalValue(Moose::Old);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueOldNeighbor() const
{
  return _neighbor_data->nodalValue(Moose::Old);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueOlder() const
{
  return _element_data->nodalValue(Moose::Older);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueOlderNeighbor() const
{
  return _neighbor_data->nodalValue(Moose::Older);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValuePreviousNL() const
{
  return _element_data->nodalValue(Moose::PreviousNL);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValuePreviousNLNeighbor() const
{
  return _neighbor_data->nodalValue(Moose::PreviousNL);
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueDot() const
{
  return _element_data->nodalValueDot();
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueDotDot() const
{
  return _element_data->nodalValueDotDot();
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueDotOld() const
{
  return _element_data->nodalValueDotOld();
}

template <typename RawOutputType>
const RawOutputType &
MooseVariableFE<RawOutputType>::nodalValueDotDotOld() const
{
  return _element_data->nodalValueDotDotOld();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeNodalValues()
{
  _element_data->computeNodalValues();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::computeNodalNeighborValues()
{
  _neighbor_data->computeNodalValues();
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::setNodalValue(const OutputType & value, unsigned int idx)
{
  _element_data->setNodalValue(value, idx);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::setDofValue(const OutputData & value, unsigned int index)
{
  _element_data->setDofValue(value, index);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::setDofValues(const DenseVector<OutputData> & values)
{
  _element_data->setDofValues(values);
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::insertNodalValue(NumericVector<Number> & residual,
                                                 const OutputData & v)
{
  _element_data->insertNodalValue(residual, v);
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::isArray() const
{
  return std::is_same<RawOutputType, RealEigenVector>::value;
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::isVector() const
{
  return std::is_same<RawOutputType, RealVectorValue>::value;
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiSecond &
MooseVariableFE<RawOutputType>::secondPhi() const
{
  return _element_data->secondPhi();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiCurl &
MooseVariableFE<RawOutputType>::curlPhi() const
{
  return _element_data->curlPhi();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiSecond &
MooseVariableFE<RawOutputType>::secondPhiFace() const
{
  return _element_data->secondPhiFace();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiCurl &
MooseVariableFE<RawOutputType>::curlPhiFace() const
{
  return _element_data->curlPhiFace();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiSecond &
MooseVariableFE<RawOutputType>::secondPhiNeighbor() const
{
  return _neighbor_data->secondPhi();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiCurl &
MooseVariableFE<RawOutputType>::curlPhiNeighbor() const
{
  return _neighbor_data->curlPhi();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiSecond &
MooseVariableFE<RawOutputType>::secondPhiFaceNeighbor() const
{
  return _neighbor_data->secondPhiFace();
}

template <typename RawOutputType>
const typename MooseVariableFE<RawOutputType>::FieldVariablePhiCurl &
MooseVariableFE<RawOutputType>::curlPhiFaceNeighbor() const
{
  return _neighbor_data->curlPhiFace();
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::usesSecondPhi() const
{
  return _element_data->usesSecondPhi();
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::usesSecondPhiNeighbor() const
{
  return _neighbor_data->usesSecondPhi();
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::computingCurl() const
{
  return _element_data->computingCurl();
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::isNodalDefined() const
{
  return _element_data->isNodalDefined();
}

template <typename RawOutputType>
bool
MooseVariableFE<RawOutputType>::isNodalNeighborDefined() const
{
  return _neighbor_data->isNodalDefined();
}

template <typename RawOutputType>
unsigned int
MooseVariableFE<RawOutputType>::oldestSolutionStateRequested() const
{
  unsigned int state = 0;
  state = std::max(state, _element_data->oldestSolutionStateRequested());
  state = std::max(state, _neighbor_data->oldestSolutionStateRequested());
  state = std::max(state, _lower_data->oldestSolutionStateRequested());
  return state;
}

template <typename RawOutputType>
void
MooseVariableFE<RawOutputType>::clearAllDofIndices()
{
  _element_data->clearDofIndices();
  _neighbor_data->clearDofIndices();
  _lower_data->clearDofIndices();
}

template class MooseVariableFE<Real>;
template class MooseVariableFE<RealVectorValue>;
template class MooseVariableFE<RealEigenVector>;
