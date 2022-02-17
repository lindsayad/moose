//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <memory>
#include "libmesh/point.h"

class InputParameters;
namespace libMesh
{
template <typename>
class TensorValue;
typedef TensorValue<Real> RealTensorValue;
}

class MooseCoordTransform
{
public:
  MooseCoordTransform(const InputParameters & params);

  libMesh::Point operator()(const libMesh::Point & point) const;

  /**
   * Set how much this coordinate system origin is translated from the canonical/reference
   * coordinate system origin. The translation vector itself should be in reference frame
   * coordinates, e.g. if this API is being used to set the translation of a multi-app based on the
   * multi-app parameter \p positions, then a point from \p positions should be passed to the main
   * application's \p MooseCoordTransform::operator() in order to get the translation in the
   * reference frame
   */
  void setTranslationVector(const libMesh::Point & translation);

private:
  std::unique_ptr<libMesh::RealTensorValue> _scale;
  std::unique_ptr<libMesh::RealTensorValue> _rotate;
  libMesh::Point _translation;
};
