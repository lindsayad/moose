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

class INSADObjectTracker;

class CGDGMaterial : public Material
{
public:
  static InputParameters validParams();

  CGDGMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  ADMaterialProperty<RealVectorValue> & _velocity;
  const ADVariableValue & _u_vel;
  const ADVariableValue & _v_vel;
  const ADVariableValue & _w_vel;
};
