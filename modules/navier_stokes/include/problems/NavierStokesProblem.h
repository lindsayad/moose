//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FEProblem.h"

class NonlinearSystem;

/**
 * Specialization of SubProblem for solving nonlinear equations plus auxiliary equations
 *
 */
class NavierStokesProblem : public FEProblem
{
public:
  static InputParameters validParams();

  NavierStokesProblem(const InputParameters & parameters);

protected:
  /**
   * Reinitialize PETSc output for proper linear/nonlinear iteration display
   */
  virtual void initPetscOutput() override;

private:
  const TagName & _velocity_mass_matrix;
};
