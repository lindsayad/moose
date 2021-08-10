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

class INSFVRhieChowInterpolator : public GeneralUserObject
{
public:
  static InputParameters validParams();
  INSFVRhieChowInterpolator(const InputParameters & params);

  void addToA(const libMesh::Elem * elem, unsigned int component, const ADReal & value);
  void addToB(const libMesh::Elem * elem, unsigned int component, const ADReal & value);
  const VectorValue<ADReal> & rcCoeff(const libMesh::Elem * elem) const;

  void initialize() override final;
  void execute() override final;
  void finalize() override final;

private:
  std::set<unsigned int> _var_numbers;
  std::unordered_set<const Elem *> _elements_to_push_pull;
  std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>> _a;
  std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>> _b;
  MooseMesh & _moose_mesh;
  const libMesh::MeshBase & _mesh;
  MooseVariableFieldBase & _u;
  MooseVariableFieldBase * const _v;
  MooseVariableFieldBase * const _w;
};
