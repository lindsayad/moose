//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluidTemperatureTimeDerivative.h"

registerMooseObject("NavierStokesApp", FluidTemperatureTimeDerivative);

InputParameters
FluidTemperatureTimeDerivative::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addRequiredCoupledVar("rho", "coupled density");
  params.addCoupledVar("pressure", "coupled pressure");
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

FluidTemperatureTimeDerivative::FluidTemperatureTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _rho(coupledValue("rho")),
    _pressure(coupledValue("pressure")),
    _pressure_dot(coupledDot("pressure")),
    _d_pressuredot_du(coupledDotDu("pressure")),
    _pressure_var_number(coupled("pressure")),
    _eos(getUserObject<SinglePhaseFluidProperties>("eos"))
{
}

Real
FluidTemperatureTimeDerivative::computeQpResidual()
{
  Real Cp = _eos.cp_from_p_T(_pressure[_qp], _u[_qp]);
  Real rho, drho_dp, drho_dT;
  _eos.rho_from_p_T(_pressure[_qp], _u[_qp], rho, drho_dp, drho_dT);
  Real enthalpy = _eos.h_from_p_T(_pressure[_qp], _u[_qp]);
  Real rho_dot = drho_dT * _u_dot[_qp] + drho_dp * _pressure_dot[_qp];

  return _rho[_qp] * Cp * TimeDerivative::computeQpResidual() + enthalpy * rho_dot * _test[_i][_qp];
}

Real
FluidTemperatureTimeDerivative::computeQpJacobian()
{
  Real Cp = _eos.cp_from_p_T(_pressure[_qp], _u[_qp]);
  Real rho, drho_dp, drho_dT;
  _eos.rho_from_p_T(_pressure[_qp], _u[_qp], rho, drho_dp, drho_dT);
  Real enthalpy = _eos.h_from_p_T(_pressure[_qp], _u[_qp]);

  return (_rho[_qp] * Cp + enthalpy * drho_dT) * TimeDerivative::computeQpJacobian();
}

Real
FluidTemperatureTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _pressure_var_number)
  {
    Real Cp = _eos.cp_from_p_T(_pressure[_qp], _u[_qp]);
    Real rho, drho_dp, drho_dT;
    _eos.rho_from_p_T(_pressure[_qp], _u[_qp], rho, drho_dp, drho_dT);
    Real enthalpy = _eos.h_from_p_T(_pressure[_qp], _u[_qp]);

    return drho_dp * Cp * TimeDerivative::computeQpResidual() +
           enthalpy * drho_dp * _d_pressuredot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  }
  else
    return 0;
}
