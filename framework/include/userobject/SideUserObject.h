//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SIDEUSEROBJECT_H
#define SIDEUSEROBJECT_H

// MOOSE includes
#include "UserObject.h"
#include "BoundaryRestrictableRequired.h"
#include "MaterialPropertyInterface.h"
#include "Coupleable.h"
#include "MooseVariableDependencyInterface.h"
#include "UserObjectInterface.h"
#include "TransientInterface.h"
#include "PostprocessorInterface.h"

// Forward Declarations
class SideUserObject;

template <>
InputParameters validParams<SideUserObject>();

class SideUserObject : public UserObject,
                       public BoundaryRestrictableRequired,
                       public MaterialPropertyInterface,
                       public Coupleable,
                       public MooseVariableDependencyInterface,
                       public UserObjectInterface,
                       public TransientInterface,
                       protected PostprocessorInterface
{
public:
  SideUserObject(const InputParameters & parameters);

protected:
  MooseMesh & _mesh;

  const MooseArray<Point> & _q_point;
  QBase *& _qrule;
  const MooseArray<ADPointReal> & _JxW;
  const MooseArray<Real> & _coord;
  const MooseArray<TypeVector<ADPointReal>> & _normals;

  const Elem *& _current_elem;
  /// current side of the current element
  unsigned int & _current_side;

  const Elem *& _current_side_elem;
  const Real & _current_side_volume;
};

#endif
