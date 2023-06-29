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

  TagID massMatrixTagID() const { return getMatrixTagID(_velocity_mass_matrix); }
  TagID BMatrixTagID() const { return getMatrixTagID(_B_matrix); }
  TagID CMatrixTagID() const { return getMatrixTagID(_C_matrix); }
  Mat getL() { return _L; }
  const std::string & velocitySplitName() const { return _velocity_split_name; }

  virtual ~NavierStokesProblem();

protected:
  /**
   * Reinitialize PETSc output for proper linear/nonlinear iteration display
   */
  virtual void initPetscOutput() override;

private:
  const TagName & _velocity_mass_matrix;
  const TagName & _B_matrix;
  const TagName & _C_matrix;
  const std::string & _velocity_split_name;

  Mat _L = nullptr;
};
