//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGConvection.h"
#include "NSFVUtils.h"

registerMooseObject("NavierStokesApp", ADDGConvection);

using namespace Moose::FV;

InputParameters
ADDGConvection::validParams()
{
  InputParameters params = ADDGKernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("velocity", "Velocity vector");
  params.addClassDescription("DG for convection");
  params.addParam<MooseEnum>(
      "advected_interp_method",
      interpolationMethods(),
      "The interpolation to use for the advected quantity. Options are "
      "'upwind', 'average', 'sou' (for second-order upwind), 'min_mod', 'vanLeer', 'quick', and "
      "'skewness-corrected' with the default being 'upwind'.");
  MooseEnum velocity_interp_method("average rc", "rc");
  return params;
}

ADDGConvection::ADDGConvection(const InputParameters & parameters)
  : ADDGKernel(parameters),
    _velocity(getADMaterialProperty<RealVectorValue>("velocity")),
    _velocity_neighbor(getNeighborADMaterialProperty<RealVectorValue>("velocity"))
{
  setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");
  if ((_advected_interp_method != InterpMethod::Average) &&
      (_advected_interp_method != InterpMethod::Upwind))
    paramError("advected_interp_method", "Currently only support 'upwind' or 'average'");
}

ADReal
ADDGConvection::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;

  auto average = [](const auto & elem_value, const auto & neighbor_value)
  { return (elem_value + neighbor_value) / 2; };

  const auto vdotn = average(_velocity[_qp], _velocity_neighbor[_qp]) * _normals[_qp];

  const auto face_u = [&]()
  {
    if (_advected_interp_method == InterpMethod::Average)
      return average(_u[_qp], _u_neighbor[_qp]);
    else
    {
      if (vdotn >= 0)
        return _u[_qp];
      else
        return _u_neighbor[_qp];
    }
    mooseError("We should never get here");
  }();

  switch (type)
  {
    case Moose::Element:
      r += vdotn * face_u * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r -= vdotn * face_u * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
