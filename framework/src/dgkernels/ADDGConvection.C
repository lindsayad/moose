//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGConvection.h"

registerMooseObject("MooseApp", ADDGConvection);

InputParameters
ADDGConvection::validParams()
{
  InputParameters params = ADDGKernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("velocity", "Velocity vector");
  params.addClassDescription("DG for convection");
  return params;
}

ADDGConvection::ADDGConvection(const InputParameters & parameters)
  : ADDGKernel(parameters),
    _velocity(getADMaterialProperty<RealVectorValue>("velocity")),
    _velocity_neighbor(getNeighborADMaterialProperty<RealVectorValue>("velocity"))
{
}

ADReal
ADDGConvection::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;

  auto average = [](const auto & elem_value, const auto & neighbor_value)
  { return (elem_value + neighbor_value) / 2; };

  const auto vdotn = average(_velocity[_qp], _velocity_neighbor[_qp]) * _normals[_qp];

  switch (type)
  {
    case Moose::Element:
      r += vdotn * average(_u[_qp], _u_neighbor[_qp]) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r -= vdotn * average(_u[_qp], _u_neighbor[_qp]) * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
