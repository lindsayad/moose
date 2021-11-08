//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "TaggingInterface.h"
#include "ADReal.h"
#include "MooseTypes.h"
#include "libmesh/vector_value.h"
#include "libmesh/id_types.h"
#include <unordered_map>
#include <set>
#include <unordered_set>

class MooseMesh;
namespace libMesh
{
class Elem;
class MeshBase;
}

class INSFVRhieChowInterpolator : public GeneralUserObject, public TaggingInterface
{
public:
  static InputParameters validParams();
  INSFVRhieChowInterpolator(const InputParameters & params);

  void addToA(const libMesh::Elem * elem, unsigned int component, const ADReal & value);
  void addToB(const libMesh::Elem * elem, unsigned int component, const ADReal & value);
  const ADReal & getB2(const libMesh::Elem & elem, unsigned int component) const;
  const VectorValue<ADReal> & getB1(const FaceInfo & fi) const;
  const VectorValue<ADReal> & getB3(const FaceInfo & fi) const;
  const VectorValue<ADReal> & rcCoeff(const libMesh::Elem * elem) const;

  void initialize() override final;
  void execute() override final;
  void finalize() override final;

private:
  void finalizeAData();
  VectorValue<ADReal> interpolateB(const VectorValue<ADReal> & b_elem,
                                   const VectorValue<ADReal> & b_neighbor,
                                   const FaceInfo & fi) const;
  void computeFirstAndSecondOverBars();
  void computeThirdOverBar();
  void applyBData();
  void finalizeBData();

  std::vector<unsigned int> _var_numbers;
  std::unordered_set<const Elem *> _elements_to_push_pull;
  std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>> _a;
  std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>> _b;
  // Here the suffix on _b refers to the number of bar operations we've performed
  std::map<const FaceInfo *, libMesh::VectorValue<ADReal>> _b1;
  std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>> _b2;
  std::map<const FaceInfo *, libMesh::VectorValue<ADReal>> _b3;
  MooseMesh & _moose_mesh;
  const libMesh::MeshBase & _mesh;
  SystemBase & _sys;
  MooseVariableFieldBase & _u;
  MooseVariableFieldBase * const _v;
  MooseVariableFieldBase * const _w;
  const VectorValue<ADReal> _example;
};

inline const ADReal &
INSFVRhieChowInterpolator::getB2(const libMesh::Elem & elem, const unsigned int component) const
{
  return libmesh_map_find(_b2, elem.id())(component);
}

inline const VectorValue<ADReal> &
INSFVRhieChowInterpolator::getB1(const FaceInfo & fi) const
{
  return libmesh_map_find(_b1, &fi);
}

inline const VectorValue<ADReal> &
INSFVRhieChowInterpolator::getB3(const FaceInfo & fi) const
{
  return libmesh_map_find(_b3, &fi);
}

inline void
INSFVRhieChowInterpolator::addToA(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  if (elem->processor_id() != this->processor_id())
    _elements_to_push_pull.insert(elem);

  _a[elem->id()](component) += value;
}

inline void
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

inline const VectorValue<ADReal> &
INSFVRhieChowInterpolator::rcCoeff(const libMesh::Elem * const elem) const
{
  const auto it = _a.find(elem->id());
  mooseAssert(it != _a.end(), "Could not find the requested element with id " << elem->id());
  return it->second;
}
