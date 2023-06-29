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
  params.addRequiredParam<TagName>(
      "B_matrix",
      "The matrix tag name corresponding to the velocity-pressure portion of the system matrix.");
  params.addRequiredParam<TagName>(
      "C_matrix",
      "The matrix tag name corresponding to the pressure-velocity portion of the system matrix.");
  params.addRequiredParam<std::string>("velocity_split_name",
                                       "The name of the velocity field split");
  return params;
}

NavierStokesProblem::NavierStokesProblem(const InputParameters & parameters)
  : FEProblem(parameters),
    _velocity_mass_matrix(getParam<TagName>("velocity_mass_matrix")),
    _B_matrix(getParam<TagName>("B_matrix")),
    _C_matrix(getParam<TagName>("C_matrix")),
    _velocity_split_name(getParam<std::string>("velocity_split_name"))
{
}

NavierStokesProblem::~NavierStokesProblem()
{
  if (_L)
    // We're destructing so don't check for errors which can throw
    MatDestroy(&_L);
}

PetscErrorCode
navierStokesKSPPreSolve(KSP ksp, Vec /*rhs*/, Vec /*x*/, void * context)
{
  KSP * subksp;
  KSP schur_ksp;
  PC fs_pc, lsc_pc;
  PetscInt num_splits;
  Mat lsc_pc_pmat, B, C, Qv, CQvdiaginv;
  Vec Qvdiaginv;
  IS velocity_is, pressure_is = NULL;
  PetscInt rstart, rend;

  PetscFunctionBegin;
  PetscCall(KSPGetPC(ksp, &fs_pc));
  // Need to call this before getting the sub ksps
  PetscCall(PCSetUp(fs_pc));
  PetscCall(PCFieldSplitGetSubKSP(fs_pc, &num_splits, &subksp));
  schur_ksp = subksp[1];
  PetscCall(KSPGetPC(schur_ksp, &lsc_pc));
  PetscCall(PCGetOperators(lsc_pc, NULL, &lsc_pc_pmat));

  auto * ns_problem = static_cast<NavierStokesProblem *>(context);
  auto Q = static_cast<PetscMatrix<Number> &>(
               ns_problem->getNonlinearSystemBase(0).getMatrix(ns_problem->massMatrixTagID()))
               .mat();
  auto B_full = static_cast<PetscMatrix<Number> &>(
                    ns_problem->getNonlinearSystemBase(0).getMatrix(ns_problem->BMatrixTagID()))
                    .mat();
  auto C_full = static_cast<PetscMatrix<Number> &>(
                    ns_problem->getNonlinearSystemBase(0).getMatrix(ns_problem->CMatrixTagID()))
                    .mat();
  auto L = ns_problem->getL();

  PetscCall(MatGetOwnershipRange(Q, &rstart, &rend));
  PetscCall(PCFieldSplitGetIS(fs_pc, ns_problem->velocitySplitName().c_str(), &velocity_is));
  PetscCall(ISComplement(velocity_is, rstart, rend, &pressure_is));
  PetscCall(MatCreateSubMatrix(Q, velocity_is, velocity_is, MAT_INITIAL_MATRIX, &Qv));
  PetscCall(MatCreateSubMatrix(B_full, velocity_is, pressure_is, MAT_INITIAL_MATRIX, &B));
  PetscCall(MatCreateSubMatrix(C_full, pressure_is, velocity_is, MAT_INITIAL_MATRIX, &C));
  PetscCall(ISDestroy(&pressure_is));

  // We'll be right-multiplying C so need right compatible
  PetscCall(MatCreateVecs(C, &Qvdiaginv, NULL));
  PetscCall(MatGetDiagonal(Qv, Qvdiaginv));
  PetscCall(VecReciprocal(Qvdiaginv));
  PetscCall(MatConvert(C, MATSAME, MAT_INITIAL_MATRIX, &CQvdiaginv));
  PetscCall(MatDiagonalScale(CQvdiaginv, NULL, Qvdiaginv));
  if (!L)
    PetscCall(MatMatMult(CQvdiaginv, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &L));
  else
    PetscCall(MatMatMult(CQvdiaginv, B, MAT_REUSE_MATRIX, PETSC_DEFAULT, &L));

  PetscCall(PetscObjectCompose((PetscObject)lsc_pc_pmat, "LSC_L", (PetscObject)L));
  PetscCall(PetscObjectCompose((PetscObject)lsc_pc_pmat, "LSC_Lp", (PetscObject)L));
  PetscCall(PetscObjectCompose((PetscObject)lsc_pc_pmat, "Q", (PetscObject)Q));

  PetscCall(VecDestroy(&Qvdiaginv));
  PetscCall(MatDestroy(&Qv));
  PetscCall(MatDestroy(&B));
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&CQvdiaginv));

  PetscFunctionReturn(PETSC_SUCCESS);
}

void
NavierStokesProblem::initPetscOutput()
{
  FEProblem::initPetscOutput();

  PetscErrorCode ierr = 0;
  KSP ksp;
  auto snes = getNonlinearSystemBase(0).getSNES();
  ierr = SNESGetKSP(snes, &ksp);
  LIBMESH_CHKERR2(this->comm(), ierr);
  ierr = KSPSetPreSolve(ksp, &navierStokesKSPPreSolve, this);
  LIBMESH_CHKERR2(this->comm(), ierr);
}
