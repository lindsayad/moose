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

class INSFVMomentumDiffusion : public FVFluxKernel,
                               public INSFVFluxKernelInterface,
                               public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  INSFVMomentumDiffusion(const InputParameters & params);
  void gatherRCData(const Elem &) override final {}
  void gatherRCData(const FaceInfo & fi) override final;
  void initialSetup() override { INSFVFluxKernelInterface::initialSetup(*this); }

protected:
  ADReal computeQpResidual() override;

  /// The dynamic viscosity
  const Moose::Functor<ADReal> & _mu;

  /// Whether we are computing RhieChow data at the time of our residual evaluation
  bool _computing_rc_data;

  /// The a coefficient for the element
  ADReal _ae = 0;

  /// The a coefficient for the neighbor
  ADReal _an = 0;
};
