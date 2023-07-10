//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VectorMass.h"

registerMooseObject("MooseApp", VectorMass);

InputParameters
VectorMass::validParams()
{
  InputParameters params = VectorKernel::validParams();
  params.addClassDescription(
      "Implements a simple mass term with weak form $(\\psi_i, \\rho u_h)$.");
  params.addParam<Real>("density", 1.0, "The $(\\rho)$ multiplier");
  params.set<MultiMooseEnum>("vector_tags") = "";
  params.set<MultiMooseEnum>("matrix_tags") = "";
  params.suppressParameter<MultiMooseEnum>("vector_tags");
  return params;
}

VectorMass::VectorMass(const InputParameters & parameters)
  : VectorKernel(parameters), _density(getParam<Real>("density"))
{
}

Real
VectorMass::computeQpResidual()
{
  return _test[_i][_qp] * _density * _u[_qp];
}

Real
VectorMass::computeQpJacobian()
{
  return _test[_i][_qp] * _density * _phi[_j][_qp];
}
