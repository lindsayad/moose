//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVSymmetryBC.h"
#include "INSFVResidualObject.h"

/**
 * A class for setting a symmetry boundary condition on the velocity. It should be
 * used in conjunction with an INSFVSymmetryPressureBC.
 */
class INSFVSymmetryVelocityBC : public INSFVSymmetryBC, public INSFVResidualObject
{
public:
  static InputParameters validParams();
  INSFVSymmetryVelocityBC(const InputParameters & params);

  void gatherRCData(const Elem &) override {}
  void gatherRCData(const FaceInfo & fi) override;

protected:
  ADReal computeQpResidual() override;

  /// x-velocity on the FaceInfo elem
  const ADVariableValue & _u_elem;
  /// y-velocity on the FaceInfo elem
  const ADVariableValue & _v_elem;
  /// z-velocity on the FaceInfo elem
  const ADVariableValue & _w_elem;

  /// x-velocity on the FaceInfo neighbor
  const ADVariableValue & _u_neighbor;
  /// y-velocity on the FaceInfo neighbor
  const ADVariableValue & _v_neighbor;
  /// z-velocity on the FaceInfo neighbor
  const ADVariableValue & _w_neighbor;

  /// What component of the velocity this object is acting on
  const unsigned int _comp;

  /// The dynamic viscosity
  const Moose::Functor<ADReal> & _mu;

  /// The mesh dimension
  const unsigned int _dim;
};
