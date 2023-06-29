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
  params.addRequiredParam<TagName>("mass_matrix",
                                   "The matrix tag name corresponding to the mass matrix.");
  params.addRequiredParam<TagName>(
      "physics_matrix",
      "The matrix tag name corresponding to just the physics portion of the system matrix.");
  params.addRequiredParam<std::string>("velocity_split_name",
                                       "The name of the velocity field split");
  return params;
}

NavierStokesProblem::NavierStokesProblem(const InputParameters & parameters)
  : FEProblem(parameters),
    _mass_matrix(getParam<TagName>("mass_matrix")),
    _physics_matrix(getParam<TagName>("physics_matrix")),
    _velocity_split_name(getParam<std::string>("velocity_split_name"))
{
}

NavierStokesProblem::~NavierStokesProblem()
{
  auto destroy_mat = [](auto mat)
  {
    if (mat)
      // We're destructing so don't check for errors which can throw
      MatDestroy(&mat);
  };
  destroy_mat(_L);
  destroy_mat(_A);
  destroy_mat(_B);
  destroy_mat(_C);
}

PetscErrorCode
navierStokesKSPPreSolve(KSP ksp, Vec /*rhs*/, Vec /*x*/, void * context)
{
  KSP * subksp;
  KSP schur_ksp;
  PC fs_pc, lsc_pc;
  PetscInt num_splits;
  Mat lsc_pc_pmat, Qv, CQvdiaginv;
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
  auto physics = static_cast<PetscMatrix<Number> &>(ns_problem->getNonlinearSystemBase(0).getMatrix(
                                                        ns_problem->physicsMatrixTagID()))
                     .mat();
  auto L = ns_problem->getL();
  auto A = ns_problem->getA();
  auto B = ns_problem->getB();
  auto C = ns_problem->getC();

  PetscCall(MatGetOwnershipRange(Q, &rstart, &rend));
  PetscCall(PCFieldSplitGetIS(fs_pc, ns_problem->velocitySplitName().c_str(), &velocity_is));
  PetscCall(ISComplement(velocity_is, rstart, rend, &pressure_is));
  PetscCall(MatCreateSubMatrix(Q, velocity_is, velocity_is, MAT_INITIAL_MATRIX, &Qv));

  if (!A)
    PetscCall(MatCreateSubMatrix(physics, velocity_is, velocity_is, MAT_INITIAL_MATRIX, &A));
  else
    PetscCall(MatCreateSubMatrix(physics, velocity_is, velocity_is, MAT_REUSE_MATRIX, &A));

  if (!B)
    PetscCall(MatCreateSubMatrix(physics, velocity_is, pressure_is, MAT_INITIAL_MATRIX, &B));
  else
    PetscCall(MatCreateSubMatrix(physics, velocity_is, pressure_is, MAT_REUSE_MATRIX, &B));

  if (!C)
    PetscCall(MatCreateSubMatrix(physics, pressure_is, velocity_is, MAT_INITIAL_MATRIX, &C));
  else
    PetscCall(MatCreateSubMatrix(physics, pressure_is, velocity_is, MAT_REUSE_MATRIX, &C));

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
  PetscCall(PetscObjectCompose((PetscObject)lsc_pc_pmat, "A", (PetscObject)A));
  PetscCall(PetscObjectCompose((PetscObject)lsc_pc_pmat, "B", (PetscObject)B));
  PetscCall(PetscObjectCompose((PetscObject)lsc_pc_pmat, "C", (PetscObject)C));

  PetscCall(VecDestroy(&Qvdiaginv));
  PetscCall(MatDestroy(&Qv));
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
