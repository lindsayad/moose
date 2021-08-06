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
#include "INSFVMomentumResidualObject.h"

/**
 * A class for setting a symmetry boundary condition on the velocity. It should be
 * used in conjunction with an INSFVSymmetryPressureBC.
 */
class INSFVSymmetryVelocityBC : public INSFVSymmetryBC, public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  INSFVSymmetryVelocityBC(const InputParameters & params);

  void gatherRCData(const Elem &) override {}
  void gatherRCData(const FaceInfo & fi) override;

protected:
  ADReal computeQpResidual() override;

  /// x-velocity
  const Moose::Functor<ADReal> & _u_functor;
  /// y-velocity
  const Moose::Functor<ADReal> & _v_functor;
  /// z-velocity
  const Moose::Functor<ADReal> & _w_functor;

  /// The dynamic viscosity
  const Moose::Functor<ADReal> & _mu;

  /// The mesh dimension
  const unsigned int _dim;

  /// Whether we are computing RC data
  bool _computing_rc_data;

  ADReal _a = 0;
};
