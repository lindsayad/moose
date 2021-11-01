//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations

/**
 * Computes |u| dt / h_min
 */
class Reynolds : public AuxKernel
{
public:
  static InputParameters validParams();

  Reynolds(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  // Velocity
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;

  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _rho;
};
