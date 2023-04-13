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
 * This is the material class used to compute enthalpy for the incompressible/weakly-compressible
 * hybrid finite-elememnt implementation of the Navier-Stokes equations
 */
class CGDGEnthalpyMaterial : public Material
{
public:
  static InputParameters validParams();

  CGDGEnthalpyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// density
  const ADMaterialProperty<Real> & _rho;

  /// the temperature
  const ADVariableValue & _temperature;

  /// the specific heat capacity
  const ADMaterialProperty<Real> & _cp;

  /// The enthalpy
  ADMaterialProperty<Real> & _rho_cp_temp;
};
