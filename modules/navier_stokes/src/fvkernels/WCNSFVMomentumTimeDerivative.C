//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFVMomentumTimeDerivative.h"
#include "SystemBase.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFVMomentumTimeDerivative);

InputParameters
WCNSFVMomentumTimeDerivative::validParams()
{
  InputParameters params = FVTimeKernel::validParams();
  params += INSFVMomentumResidualObject::validParams();
  params.addClassDescription(
      "Adds the time derivative term to the incompressible Navier-Stokes momentum equation.");
  params.addRequiredParam<MaterialPropertyName>(NS::density, "The density material property");
  params.addRequiredParam<MaterialPropertyName>(
      NS::time_deriv(NS::density), "The time derivative of the density material property");
  return params;
}

WCNSFVMomentumTimeDerivative::WCNSFVMomentumTimeDerivative(const InputParameters & params)
  : FVTimeKernel(params),
    INSFVMomentumResidualObject(*this),
    _rho(getFunctor<ADReal>(NS::density)),
    _rho_dot(getFunctor<ADReal>(NS::time_deriv(NS::density)))
{
}

ADReal
WCNSFVMomentumTimeDerivative::computeQpResidual(const Elem & elem)
{
  const auto elem_arg = makeElemArg(&elem);
  const auto rho_dot = _rho_dot(elem_arg);
  const auto var_dot = _var.dot(elem_arg);
  const auto rho = _rho(elem_arg);
  const auto var = _var(elem_arg);

  if (_computing_rc_data)
  {
    const auto dof_number = elem.dof_number(_sys.number(), _var.number(), 0);
    mooseAssert(var.derivatives()[dof_number] == 1.,
                "This is an implicit assumption in our coefficient calculation.");
    _a = rho_dot;
    _a += rho * var_dot.derivatives()[dof_number];
  }

  return rho_dot * var + rho * var_dot;
}

ADReal
WCNSFVMomentumTimeDerivative::computeQpResidual()
{
  return computeQpResidual(*_current_elem);
}

void
WCNSFVMomentumTimeDerivative::gatherRCData(const Elem & elem)
{
  // _rho and _rho_dot could ultimately be functions of the nonlinear variables making our residual
  // nonlinear so we cannot do the simple treatment that is done in
  // INSFVMomentumTimeDerivative::gatherRCData

  _computing_rc_data = true;
  // Fill-in the coefficient _a (but without multiplication by V)
  computeQpResidual(elem);
  _computing_rc_data = false;

  _rc_uo.addToA(&elem, _index, _a * _assembly.elementVolume(&elem));
}
