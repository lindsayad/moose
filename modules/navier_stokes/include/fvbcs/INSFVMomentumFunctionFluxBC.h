//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFluxBC.h"
#include "INSFVMomentumResidualObject.h"

class InputParameters;
class Function;

/**
 * A parent class for slip/no-slip wall boundary conditions
 */
class INSFVMomentumFunctionFluxBC : public FVFluxBC, public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  INSFVMomentumFunctionFluxBC(const InputParameters & params);

  void gatherRCData(const FaceInfo &) override {}
  void gatherRCData(const Elem &) override {}

protected:
  ADReal computeQpResidual() override;

  const Function & _function;
};
