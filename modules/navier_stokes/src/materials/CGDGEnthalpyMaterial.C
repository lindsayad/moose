//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CGDGEnthalpyMaterial.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", CGDGEnthalpyMaterial);

InputParameters
CGDGEnthalpyMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "This is the material class used to compute enthalpy for "
      "the incompressible/weakly-compressible hybrid finite-elememnt implementation "
      "of the Navier-Stokes equations.");
  params.addRequiredParam<MaterialPropertyName>(NS::density, "The value for the density");
  params.addRequiredCoupledVar(NS::temperature, "The temperature");
  params.addRequiredParam<MaterialPropertyName>(NS::cp, "The name of the specific heat capacity");
  return params;
}

CGDGEnthalpyMaterial::CGDGEnthalpyMaterial(const InputParameters & parameters)
  : Material(parameters),
    _rho(getADMaterialProperty<Real>(NS::density)),
    _temperature(adCoupledValue(NS::temperature)),
    _cp(getADMaterialProperty<Real>(NS::cp)),
    _rho_cp_temp(declareADProperty<Real>("rho_cp_temp"))
{
}

void
CGDGEnthalpyMaterial::computeQpProperties()
{
  _rho_cp_temp[_qp] = _rho[_qp] * _cp[_qp] * _temperature[_qp];
}
