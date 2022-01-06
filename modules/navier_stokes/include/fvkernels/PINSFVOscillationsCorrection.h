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
#include "INSFVFluxKernelInterface.h"
#include "INSFVMomentumResidualObject.h"

class PINSFVOscillationsCorrection : public FVFluxKernel,
                                     public INSFVFluxKernelInterface,
                                     public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  PINSFVOscillationsCorrection(const InputParameters & params);
  void gatherRCData(const Elem &) override final {}
  void gatherRCData(const FaceInfo & fi) override final;
  void initialSetup() override { INSFVFluxKernelInterface::initialSetup(*this); }

protected:
  ADReal computeQpResidual() override;

  /// Darcy coefficient
  const Moose::Functor<ADRealVectorValue> * const _cL;
  /// Forchheimer coefficient
  const Moose::Functor<ADRealVectorValue> * const _cQ;

  /// Booleans to select the right models
  const bool _use_Darcy_friction_model;
  const bool _use_Forchheimer_friction_model;

  /// Porosity to compute the intersitial velocity from the superficial velocity
  const Moose::Functor<ADReal> & _eps;
  /// Density as a functor
  const Moose::Functor<ADReal> & _rho;

  /// Whether we are computing RhieChow data at the time of our residual evaluation
  bool _computing_rc_data;

  /// The a coefficient for the element
  ADReal _ae = 0;

  /// The a coefficient for the neighbor
  ADReal _an = 0;

  /// Parameter for scaling the consistent pressure interpolation
  Real _consistent_scaling;
};
