//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADDGConvectionBC.h"

/**
 * This class computes the mass equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes mass
 * equation at inflow (and possibly outflow) boundaries
 */
class CGMassBC : public ADDGConvectionBC
{
public:
  static InputParameters validParams();

  CGMassBC(const InputParameters & parameters);
};
