//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluidVelocityTimeDerivative.h"

registerMooseObject("NavierStokesApp", FluidVelocityTimeDerivative);

InputParameters
FluidVelocityTimeDerivative::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addRequiredCoupledVar("rho", "coupled density");
  params.addRequiredCoupledVar("pressure", "coupled pressure");
  params.addRequiredCoupledVar("temperature", "coupled temperature");
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

FluidVelocityTimeDerivative::FluidVelocityTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _rho(coupledValue("rho")),
    _pressure(coupledValue("pressure")),
    _pressure_dot(coupledDot("pressure")),
    _d_pressuredot_du(coupledDotDu("pressure")),
    _pressure_var_number(coupled("pressure")),
    _temperature(coupledValue("temperature")),
    _temperature_dot(coupledDot("temperature")),
    _d_temperaturedot_du(coupledDotDu("temperature")),
    _temperature_var_number(coupled("temperature")),
    _eos(getUserObject<SinglePhaseFluidProperties>("eos"))
{
}

Real
FluidVelocityTimeDerivative::computeQpResidual()
{
  Real rho, drho_dp, drho_dT;
  _eos.rho_from_p_T(_pressure[_qp], _temperature[_qp], rho, drho_dp, drho_dT);
  Real rho_dot = drho_dT * _temperature_dot[_qp] + drho_dp * _pressure_dot[_qp];

  return _rho[_qp] * TimeDerivative::computeQpResidual() + _u[_qp] * rho_dot * _test[_i][_qp];
}

Real
FluidVelocityTimeDerivative::computeQpJacobian()
{
  Real rho, drho_dp, drho_dT;
  _eos.rho_from_p_T(_pressure[_qp], _temperature[_qp], rho, drho_dp, drho_dT);
  Real rho_dot = drho_dT * _temperature_dot[_qp] + drho_dp * _pressure_dot[_qp];

  return _rho[_qp] * TimeDerivative::computeQpJacobian() + _phi[_j][_qp] * rho_dot * _test[_i][_qp];
}

Real
FluidVelocityTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _temperature_var_number)
  {
    Real rho, drho_dp, drho_dT;
    _eos.rho_from_p_T(_pressure[_qp], _temperature[_qp], rho, drho_dp, drho_dT);
    return drho_dT * _phi[_j][_qp] * TimeDerivative::computeQpResidual() +
           drho_dT * _u[_qp] * _d_temperaturedot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  }
  else if (jvar == _pressure_var_number)
  {
    Real rho, drho_dp, drho_dT;
    _eos.rho_from_p_T(_pressure[_qp], _temperature[_qp], rho, drho_dp, drho_dT);
    return drho_dp * _phi[_j][_qp] * TimeDerivative::computeQpResidual() +
           drho_dp * _u[_qp] * _d_pressuredot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  }
  else
    return 0;
}
