/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "NeighborCoupleableMooseVariableDependencyIntermediateInterface.h"
#include "MooseObject.h"
#include "InputParameters.h"
#include "SubProblem.h"

NeighborCoupleableMooseVariableDependencyIntermediateInterface::
    NeighborCoupleableMooseVariableDependencyIntermediateInterface(const MooseObject * moose_object,
                                                                   bool nodal,
                                                                   bool neighbor_nodal)
  : NeighborCoupleable(moose_object, nodal, neighbor_nodal), ScalarCoupleable(moose_object)
{
  const std::vector<MooseVariableFE *> & coupled_vars = getCoupledMooseVars();
  for (unsigned int i = 0; i < coupled_vars.size(); i++)
    addMooseVariableDependency(coupled_vars[i]);

  const InputParameters & parameters = moose_object->parameters();

  SubProblem & problem = *parameters.get<SubProblem *>("_subproblem");

  THREAD_ID tid = parameters.get<THREAD_ID>("_tid");

  // Try the scalar version first
  std::string variable_name = parameters.getMooseType("variable");
  if (variable_name == "")
    // When using vector variables, we are only going to use the first one in the list at the
    // interface level...
    variable_name = parameters.getVecMooseType("variable")[0];

  addMooseVariableDependency(&problem.getVariable(tid, variable_name));
}
