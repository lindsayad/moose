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
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  return params;
}

CGMass::CGMass(const InputParameters & parameters)
  : ADKernel(parameters),
    _grad_u_vel(adCoupledGradient("u")),
    _grad_v_vel(isCoupled("v") ? adCoupledGradient("v") : _ad_grad_zero),
    _grad_w_vel(isCoupled("w") ? adCoupledGradient("w") : _ad_grad_zero)
{
}

ADReal
CGMass::computeQpResidual()
{
  // (div u) * q
  // Note: we (arbitrarily) multiply this term by -1 so that it matches the -p(div v)
  // term in the momentum equation.  Not sure if that is really important?
  return -(_grad_u_vel[_qp](0) + _grad_v_vel[_qp](1) + _grad_w_vel[_qp](2)) * _test[_i][_qp];
}
