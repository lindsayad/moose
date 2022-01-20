//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMaterial.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", INSFVMaterial);

InputParameters
INSFVMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("This is the material class used to compute advected quantities for "
                             "the finite-volume implementation of the Navier-Stokes equations.");
  params.addRequiredParam<MaterialPropertyName>("rho", "The value for the density");
  params.addCoupledVar("temperature", "the temperature");
  params.addParam<MaterialPropertyName>("cp_name", "cp", "the name of the specific heat capacity");
  return params;
}

INSFVMaterial::INSFVMaterial(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _rho(getFunctor<ADReal>("rho")),
    _has_temperature(isParamValid("temperature")),
    _temperature(_has_temperature ? getVarHelper<MooseVariableFVReal>("temperature", 0) : nullptr),
    _cp(_has_temperature ? &getFunctor<ADReal>("cp_name") : nullptr)
{
  if (_has_temperature)
    addFunctorProperty<ADReal>("rho_cp_temp",
                               [this](const auto & r, const auto & t) -> ADReal
                               { return _rho(r, t) * (*_cp)(r, t) * (*_temperature)(r, t); });
}
