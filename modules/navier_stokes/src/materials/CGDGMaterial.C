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
  return params;
}

CGDGMaterial::CGDGMaterial(const InputParameters & parameters)
  : Material(parameters),
    _velocity(declareADProperty<RealVectorValue>("velocity")),
    _u_vel(adCoupledValue("u")),
    _v_vel(isCoupled("v") ? adCoupledValue("v") : _ad_zero),
    _w_vel(isCoupled("w") ? adCoupledValue("w") : _ad_zero)
{
}

void
CGDGMaterial::computeQpProperties()
{
  _velocity[_qp] = ADRealVectorValue{_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]};
}
