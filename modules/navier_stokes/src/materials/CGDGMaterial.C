//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CGDGMaterial.h"
#include "Function.h"
#include "Assembly.h"
#include "INSADObjectTracker.h"
#include "FEProblemBase.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", CGDGMaterial);

InputParameters
CGDGMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredParam<MaterialPropertyName>(NS::density, "The density");
  return params;
}

CGDGMaterial::CGDGMaterial(const InputParameters & parameters)
  : Material(parameters),
    _velocity(declareADProperty<RealVectorValue>(NS::velocity)),
    _vel_x(adCoupledValue("u")),
    _vel_y(isCoupled("v") ? adCoupledValue("v") : _ad_zero),
    _vel_z(isCoupled("w") ? adCoupledValue("w") : _ad_zero),
    _mom_x(declareADProperty<Real>(NS::momentum_x)),
    _mom_y(declareADProperty<Real>(NS::momentum_y)),
    _mom_z(declareADProperty<Real>(NS::momentum_z)),
    _rho(getADMaterialProperty<Real>(NS::density))
{
}

void
CGDGMaterial::computeQpProperties()
{
  _velocity[_qp] = ADRealVectorValue{_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]};
  _mom_x[_qp] = _rho[_qp] * _vel_x[_qp];
  _mom_y[_qp] = _rho[_qp] * _vel_y[_qp];
  _mom_z[_qp] = _rho[_qp] * _vel_z[_qp];
}
