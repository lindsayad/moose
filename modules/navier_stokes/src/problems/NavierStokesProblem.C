//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NavierStokesProblem.h"
#include "NonlinearSystemBase.h"
#include "libmesh/petsc_matrix.h"
#include <petscsnes.h>

registerMooseObject("NavierStokesApp", NavierStokesProblem);

InputParameters
NavierStokesProblem::validParams()
{
  InputParameters params = FEProblem::validParams();
  params.addRequiredParam<TagName>(
      "velocity_mass_matrix", "The matrix tag name corresponding to the velocity mass matrix.");
  return params;
}

NavierStokesProblem::NavierStokesProblem(const InputParameters & parameters)
  : FEProblem(parameters), _velocity_mass_matrix(getParam<TagName>("velocity_mass_matrix"))
{
}

void
NavierStokesProblem::initPetscOutput()
{
  PetscErrorCode ierr = 0;
  KSP ksp, schur_ksp;
  KSP * subksp;
  PC fs_pc, lsc_pc;
  PetscInt num_splits;
  Mat lsc_pc_pmat;
  auto snes = getNonlinearSystemBase(0).getSNES();
  ierr = SNESGetKSP(snes, &ksp);
  LIBMESH_CHKERR2(this->comm(), ierr);
  ierr = KSPGetPC(ksp, &fs_pc);
  LIBMESH_CHKERR2(this->comm(), ierr);
  // Need to call this before getting the sub ksps
  ierr = PCSetUp(fs_pc);
  LIBMESH_CHKERR2(this->comm(), ierr);
  ierr = PCFieldSplitGetSubKSP(fs_pc, &num_splits, &subksp);
  LIBMESH_CHKERR2(this->comm(), ierr);
  schur_ksp = subksp[1];
  ierr = KSPGetPC(schur_ksp, &lsc_pc);
  LIBMESH_CHKERR2(this->comm(), ierr);
  ierr = PCGetOperators(lsc_pc, NULL, &lsc_pc_pmat);
  LIBMESH_CHKERR2(this->comm(), ierr);

  const auto mass_matrix_tag_id = getMatrixTagID(_velocity_mass_matrix);
  auto & libmesh_sparse_mass_matrix = getNonlinearSystemBase(0).getMatrix(mass_matrix_tag_id);
  auto & libmesh_petsc_mass_matrix = static_cast<PetscMatrix<Number> &>(libmesh_sparse_mass_matrix);
  auto petsc_mass_matrix = libmesh_petsc_mass_matrix.mat();

  ierr = PetscObjectCompose((PetscObject)lsc_pc_pmat, "Qv", (PetscObject)petsc_mass_matrix);
  LIBMESH_CHKERR2(this->comm(), ierr);
}
