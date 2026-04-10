//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MultiAppMFEMCopyTransfer.h"
#include "FEProblemBase.h"
#include "MultiApp.h"
#include "SystemBase.h"
#include "MFEMProblem.h"
#include "MFEMMesh.h"

registerMooseObject("MooseApp", MultiAppMFEMCopyTransfer);

InputParameters
MultiAppMFEMCopyTransfer::validParams()
{
  InputParameters params = MultiAppTransfer::validParams();
  params.addRequiredParam<std::vector<AuxVariableName>>(
      "variable", "AuxVariable to store transferred value in.");
  params.addRequiredParam<std::vector<VariableName>>("source_variable",
                                                     "Variable to transfer from");
  params.addClassDescription("Copies variable values from one MFEM application to another");
  return params;
}

MultiAppMFEMCopyTransfer::MultiAppMFEMCopyTransfer(InputParameters const & params)
  : MultiAppTransfer(params),
    _from_var_names(getParam<std::vector<VariableName>>("source_variable")),
    _to_var_names(getParam<std::vector<AuxVariableName>>("variable"))
{
  auto bad_problem = [this]()
  {
    mooseError(type(),
               " only works with MFEMProblem based applications. Check that all your inputs "
               "involved in this transfer are MFEMProblem based");
  };
  if (hasToMultiApp())
  {
    if (!dynamic_cast<MFEMProblem *>(&getToMultiApp()->problemBase()))
      bad_problem();
    for (const auto i : make_range(getToMultiApp()->numGlobalApps()))
      if (getToMultiApp()->hasLocalApp(i) &&
          !dynamic_cast<MFEMProblem *>(&getToMultiApp()->appProblemBase(i)))
        bad_problem();
  }
  if (hasFromMultiApp())
  {
    if (!dynamic_cast<MFEMProblem *>(&getFromMultiApp()->problemBase()))
      bad_problem();
    for (const auto i : make_range(getFromMultiApp()->numGlobalApps()))
      if (getFromMultiApp()->hasLocalApp(i) &&
          !dynamic_cast<MFEMProblem *>(&getFromMultiApp()->appProblemBase(i)))
        bad_problem();
  }
}

void
MultiAppMFEMCopyTransfer::transfer(MFEMProblem & to_problem, MFEMProblem & from_problem)
{
  if (numToVar() != numFromVar())
    mooseError("Number of variables transferred must be same in both systems.");

  // This transfer copies local MFEM vectors directly and requires both apps to have
  // identical local DoF layouts. That is only guaranteed when both apps build their
  // mfem::ParMesh via the same path: either both via MFEMMesh (native MFEM reader) or
  // both via LibMesh mesh conversion. Mixing the two paths can produce different METIS
  // partitions and therefore different local DoF counts.
  const bool from_is_mfem_mesh = (from_problem.mesh().type() == "MFEMMesh");
  const bool to_is_mfem_mesh = (to_problem.mesh().type() == "MFEMMesh");
  if (from_is_mfem_mesh != to_is_mfem_mesh)
    mooseError(type(),
               " requires both apps to build their MFEM mesh via the same path. "
               "Either both should use MFEMMesh or both should use a LibMesh-based mesh "
               "(e.g. FileMesh). Mixing the two gives different METIS partitions and "
               "mismatched local DoF counts.");

  auto getGF = [&](MFEMProblem & problem, const std::string & name) -> mfem::Vector &
  {
    if (problem.getProblemData().gridfunctions.Has(name))
      return *problem.getGridFunction(name);
    if (problem.getProblemData().cmplx_gridfunctions.Has(name))
      return *problem.getComplexGridFunction(name);
    mooseError("No real or complex variable named '", name, "' found.");
  };

  for (const auto v : make_range(numToVar()))
  {
    mfem::Vector & from_var = getGF(from_problem, getFromVarName(v));
    mfem::Vector & to_var = getGF(to_problem, getToVarName(v));

    if (from_var.Size() != to_var.Size())
      mooseError("'", getFromVarName(v), "' and '", getToVarName(v), "' differ in no. of DoFs.");

    to_var = from_var;
  }
}

void
MultiAppMFEMCopyTransfer::execute()
{
  TIME_SECTION("MultiAppMFEMCopyTransfer::execute", 5, "Copies variables");
  if (_current_direction == TO_MULTIAPP)
  {
    for (unsigned int i = 0; i < getToMultiApp()->numGlobalApps(); i++)
    {
      if (getToMultiApp()->hasLocalApp(i))
      {
        transfer(static_cast<MFEMProblem &>(getToMultiApp()->appProblemBase(i)),
                 static_cast<MFEMProblem &>(getToMultiApp()->problemBase()));
      }
    }
  }
  else if (_current_direction == FROM_MULTIAPP)
  {
    for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); i++)
    {
      if (getFromMultiApp()->hasLocalApp(i))
      {
        transfer(static_cast<MFEMProblem &>(getFromMultiApp()->problemBase()),
                 static_cast<MFEMProblem &>(getFromMultiApp()->appProblemBase(i)));
      }
    }
  }
  else if (_current_direction == BETWEEN_MULTIAPP)
  {
    int transfers_done = 0;
    for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); i++)
    {
      if (getFromMultiApp()->hasLocalApp(i) && getToMultiApp()->hasLocalApp(i))
      {
        transfer(static_cast<MFEMProblem &>(getToMultiApp()->appProblemBase(i)),
                 static_cast<MFEMProblem &>(getFromMultiApp()->appProblemBase(i)));
        ++transfers_done;
      }
    }
    if (!transfers_done)
      mooseError("BETWEEN_MULTIAPP transfer not supported if there is not at least one subapp "
                 "per multiapp involved on each rank");
  }
}

void
MultiAppMFEMCopyTransfer::checkSiblingsTransferSupported() const
{
  // Check that we are in the supported configuration: same number of source and target apps
  // The allocation of the child apps on the processors must be the same
  if (getFromMultiApp()->numGlobalApps() == getToMultiApp()->numGlobalApps())
  {
    for (const auto i : make_range(getToMultiApp()->numGlobalApps()))
      if (getFromMultiApp()->hasLocalApp(i) + getToMultiApp()->hasLocalApp(i) == 1)
        mooseError("Child application allocation on parallel processes must be the same to support "
                   "siblings variable field copy transfer");
  }
  else
    mooseError("Number of source and target child apps must match for siblings transfer");
}

#endif
