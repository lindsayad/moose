//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeDerivative.h"
#include "SinglePhaseFluidProperties.h"

class FluidVelocityTimeDerivative : public TimeDerivative
{
public:
  static InputParameters validParams();

  FluidVelocityTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  // Coupled variables
  const VariableValue & _rho;
  const VariableValue & _pressure;
  const VariableValue & _pressure_dot;
  const VariableValue & _d_pressuredot_du;
  unsigned _pressure_var_number;
  const VariableValue & _temperature;
  const VariableValue & _temperature_dot;
  const VariableValue & _d_temperaturedot_du;
  unsigned _temperature_var_number;
  const SinglePhaseFluidProperties & _eos;
};
