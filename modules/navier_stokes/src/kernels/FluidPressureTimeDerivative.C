//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluidPressureTimeDerivative.h"

registerMooseObject("NavierStokesApp", FluidPressureTimeDerivative);

InputParameters
FluidPressureTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addRequiredCoupledVar("temperature", "coupled temperature");
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

FluidPressureTimeDerivative::FluidPressureTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
    _temperature(coupledValue("temperature")),
    _temperature_dot(coupledDot("temperature")),
    _d_temperaturedot_du(coupledDotDu("temperature")),
    _temperature_var_number(coupled("temperature")),
    _eos(getUserObject<SinglePhaseFluidProperties>("eos"))
{
}

Real
FluidPressureTimeDerivative::computeQpResidual()
{
  Real rho, drho_dp, drho_dT;
  _eos.rho_from_p_T(_u[_qp], _temperature[_qp], rho, drho_dp, drho_dT);
  return (drho_dT * _temperature_dot[_qp] + drho_dp * _u_dot[_qp]) * _test[_i][_qp];
}

Real
FluidPressureTimeDerivative::computeQpJacobian()
{
  Real rho, drho_dp, drho_dT;
  _eos.rho_from_p_T(_u[_qp], _temperature[_qp], rho, drho_dp, drho_dT);
  return drho_dp * _du_dot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}

Real
FluidPressureTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _temperature_var_number)
  {
    Real rho, drho_dp, drho_dT;
    _eos.rho_from_p_T(_u[_qp], _temperature[_qp], rho, drho_dp, drho_dT);
    return drho_dT * _d_temperaturedot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  }
  else
    return 0.;
}
