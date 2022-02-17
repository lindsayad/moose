//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseCoordTransform.h"
#include "InputParameters.h"

#include "libmesh/tensor_value.h"
#include "libmesh/point.h"

using namespace libMesh;

MooseCoordTransform::MooseCoordTransform(const InputParameters & params)
{
  // Rotation
  const bool has_alpha = params.isParamValid("alpha_rotation");
  const bool has_beta = params.isParamValid("beta_rotation");
  const bool has_gamma = params.isParamValid("gamma_rotation");
  if (has_alpha || has_beta || has_gamma)
  {
    const auto alpha = (has_alpha ? params.get<Real>("alpha_rotation") : Real(0));
    const auto beta = (has_beta ? params.get<Real>("beta_rotation") : Real(0));
    const auto gamma = (has_gamma ? params.get<Real>("gamma_rotation") : Real(0));

    // We actually want the inverse of the rotation matrix since the user is telling us how much
    // they are rotated with respect to the canonical coordinate system
    _rotate = std::make_unique<RealTensorValue>(
        RealTensorValue::inverse_rotation_matrix(alpha, beta, gamma));
  }

  // Scale
  if (params.isParamValid("length_units_per_meter"))
  {
    const auto scale = 1 / params.get<Real>("length_units_per_meter");

    _scale =
        std::make_unique<RealTensorValue>(RealTensorValue(scale, 0, 0, 0, scale, 0, 0, 0, scale));
  }
}

Point
MooseCoordTransform::operator()(const Point & point) const
{
  Point ret;
  if (_rotate)
    ret = (*_rotate) * point;
  if (_scale)
    ret = (*_scale) * ret;

  // If this shows up in profiling we can make _translation a pointer
  ret += _translation;

  return ret;
}

void
MooseCoordTransform::setTranslationVector(const Point & translation)
{
  // Recall that we are defining how to *get back* to the reference frame and the input describes
  // how much we are translated *from* the reference frame
  _translation = -translation;
}
