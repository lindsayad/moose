//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumFunctionFluxBC.h"
#include "InputParameters.h"
#include "Function.h"

registerMooseObject("NavierStokesApp", INSFVMomentumFunctionFluxBC);

InputParameters
INSFVMomentumFunctionFluxBC::validParams()
{
  auto params = FVFluxBC::validParams();
  params += INSFVMomentumResidualObject::validParams();
  params.addClassDescription("Allows implementation of a generic function on a boundary for a flux "
                             "BC. This is useful in MMS studies.");
  params.addRequiredParam<FunctionName>("function", "The value of the flux crossing the boundary.");
  return params;
}

INSFVMomentumFunctionFluxBC::INSFVMomentumFunctionFluxBC(const InputParameters & params)
  : FVFluxBC(params), INSFVMomentumResidualObject(*this), _function(getFunction("function"))
{
}

ADReal
INSFVMomentumFunctionFluxBC::computeQpResidual()
{
  return _function.value(_t, _face_info->faceCentroid());
}
