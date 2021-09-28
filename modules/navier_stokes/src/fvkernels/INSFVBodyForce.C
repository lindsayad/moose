//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVBodyForce.h"

registerMooseObject("NavierStokesApp", INSFVBodyForce);

InputParameters
INSFVBodyForce::validParams()
{
  auto params = FVBodyForce::validParams();
  params += INSFVMomentumResidualObject::validParams();
  params.addClassDescription("Body force that contributes to the Rhie-Chow interpolation");
  return params;
}

INSFVBodyForce::INSFVBodyForce(const InputParameters & params)
  : FVBodyForce(params), INSFVMomentumResidualObject(*this)
{
}
