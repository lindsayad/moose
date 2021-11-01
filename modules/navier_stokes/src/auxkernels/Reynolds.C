//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Reynolds.h"
#include "MooseMesh.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", Reynolds);

InputParameters
Reynolds::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription("Computes |u| dt / h_min.");
  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D
  params.addParam<MaterialPropertyName>(NS::mu, NS::mu, "The name of the dynamic viscosity");
  params.addParam<MaterialPropertyName>(NS::density, NS::density, "The name of the density");

  return params;
}

Reynolds::Reynolds(const InputParameters & parameters)
  : AuxKernel(parameters),
    _u_vel(coupledValue("u")),
    _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
    _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero),
    // Material properties
    _mu(getMaterialProperty<Real>(NS::mu)),
    _rho(getMaterialProperty<Real>(NS::density))
{
}

Real
Reynolds::computeValue()
{
  const RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real vel_mag = U.norm();

  return _rho[_qp] * vel_mag * _current_elem->hmax() / _mu[_qp];
}
