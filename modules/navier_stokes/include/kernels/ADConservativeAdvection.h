//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
class ADConservativeAdvection : public ADKernel
{
public:
  static InputParameters validParams();

  ADConservativeAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// advection velocity
  const ADMaterialProperty<RealVectorValue> & _velocity;

  /// advected quantity
  const MooseArray<ADReal> & _adv_quant;
};
