//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADIntegratedBC.h"

class Function;

class ADDGConvectionBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADDGConvectionBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// The velocity as a material property
  const ADMaterialProperty<RealVectorValue> * const _velocity_mat_prop;

  /// The velocity as a function
  const Function * const _velocity_function;

  /// The advected quantity
  const MooseArray<ADReal> & _adv_quant;

  /// Dirichlet value for the primal variable
  const Function * const _primal_dirichlet;

  /// The density
  const ADMaterialProperty<Real> * const _rho;

  /// The specific heat capacity
  const ADMaterialProperty<Real> * const _cp;
};
