//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * A class that computes material properties relevant to a CGDG hybrid discretization. In fact the
 * material properties are very generic and could also be applied to a pure CG scheme, however, the
 * CG (AD) implementation is currently covered by \p INSADMaterial which is constructed in a way
 * naturally geared towards a Petrov-Galerkin stabilization
 */
class CGDGMaterial : public Material
{
public:
  static InputParameters validParams();

  CGDGMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// The velocity
  ADMaterialProperty<RealVectorValue> & _velocity;
  /// The x-component of velocity
  const ADVariableValue & _vel_x;
  /// The y-component of velocity
  const ADVariableValue & _vel_y;
  /// The z-component of velocity
  const ADVariableValue & _vel_z;
  /// The x-component of momentum
  ADMaterialProperty<Real> & _mom_x;
  /// The y-component of momentum
  ADMaterialProperty<Real> & _mom_y;
  /// The z-component of momentum
  ADMaterialProperty<Real> & _mom_z;
  /// The density
  const ADMaterialProperty<Real> & _rho;
};
