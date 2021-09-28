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
  void gatherRCData(const Elem &) override {}
  void gatherRCData(const FaceInfo & fi) override;
  void initialSetup() override { INSFVFluxKernelInterface::initialSetup(*this); }

protected:
  ADReal computeQpResidual() override;

  /// The dynamic viscosity
  const Moose::Functor<ADReal> & _mu;
};
