//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGMomentumPressure.h"
#include "Function.h"

registerMooseObject("NavierStokesApp", DGMomentumPressure);

InputParameters
DGMomentumPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredCoupledVar(NS::pressure, "The pressure variable");
  params.addRequiredParam<unsigned short>("component", "The velocity component.");
  return params;
}

DGMomentumPressure::DGMomentumPressure(const InputParameters & parameters)
  : ADKernel(parameters),
    _grad_pressure(adCoupledGradient(NS::pressure)),
    _component(getParam<unsigned short>("component"))
{
}

ADReal
DGMomentumPressure::computeQpResidual()
{
  return _grad_pressure[_qp](_component) * _test[_i][_qp];
}
