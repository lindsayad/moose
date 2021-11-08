//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFluxKernel.h"
#include "INSFVMomentumResidualObject.h"

// Forward declare variable class
class INSFVVelocityVariable;

class INSFVMixingLengthReynoldsStress : public FVFluxKernel, public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();

  INSFVMixingLengthReynoldsStress(const InputParameters & params);

  void gatherRCData(const Elem &) override final {}
  void gatherRCData(const FaceInfo &) override final;

protected:
  ADReal computeQpResidual() override;

  /// the dimension of the simulation
  const unsigned int _dim;

  /// index x|y|z
  const unsigned int _axis_index;

  /// x-velocity
  const INSFVVelocityVariable * const _u_var;
  /// y-velocity
  const INSFVVelocityVariable * const _v_var;
  /// z-velocity
  const INSFVVelocityVariable * const _w_var;

  /// Density
  const Moose::Functor<ADReal> & _rho;

  /// Turbulent eddy mixing length
  const Moose::Functor<ADReal> & _mixing_len;

  /// Whether we are currently computing Rhie-Chow data
  bool _computing_rc_data = false;

  /// Rhie-Chow element coefficient
  ADReal _ae = 0;

  /// Rhie-Chow neighbor coefficient
  ADReal _an = 0;
};
