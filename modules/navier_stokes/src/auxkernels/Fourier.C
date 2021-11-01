//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Fourier.h"
#include "MooseMesh.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", Fourier);

InputParameters
Fourier::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addParam<MaterialPropertyName>(NS::mu, NS::mu, "The name of the dynamic viscosity");
  params.addParam<MaterialPropertyName>(NS::density, NS::density, "The name of the density");
  return params;
}

Fourier::Fourier(const InputParameters & parameters)
  : AuxKernel(parameters),
    // Material properties
    _mu(getMaterialProperty<Real>(NS::mu)),
    _rho(getMaterialProperty<Real>(NS::density))
{
}

Real
Fourier::computeValue()
{
  const auto nu = _mu[_qp] / _rho[_qp];
  const auto h = _current_elem->hmin();
  return _dt * nu / (h * h);
}
