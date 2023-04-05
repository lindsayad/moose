//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CGMass.h"
#include "Function.h"

registerMooseObject("NavierStokesApp", CGMass);

InputParameters
CGMass::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("velocity", "The velocity");
  return params;
}

CGMass::CGMass(const InputParameters & parameters)
  : ADKernel(parameters), _velocity(getADMaterialProperty<RealVectorValue>("velocity"))
{
}

ADReal
CGMass::computeQpResidual()
{
  // (div u) * q
  // Note: we (arbitrarily) multiply this term by -1 so that it matches the -p(div v)
  // term in the momentum equation.  Not sure if that is really important?
  return _velocity[_qp] * _grad_test[_i][_qp];
}
