//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVTimeKernel.h"
#include "INSFVMomentumResidualObject.h"

class INSFVMomentumTimeDerivative : public FVTimeKernel, public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  INSFVMomentumTimeDerivative(const InputParameters & params);

  void gatherRCData(const Elem &) override;
  void gatherRCData(const FaceInfo &) override {}

protected:
  ADReal computeQpResidual(const Elem & elem);
  ADReal computeQpResidual() override final;

  const Real & _rho;
};
