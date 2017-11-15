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

#include "NeighborCoupleable.h"

#include "FEProblem.h"
#include "MooseError.h" // mooseDeprecated
#include "MooseVariableField.h"
#include "Problem.h"
#include "SubProblem.h"

NeighborCoupleable::NeighborCoupleable(const MooseObject * moose_object,
                                       bool nodal,
                                       bool neighbor_nodal)
  : Coupleable(moose_object, nodal), _neighbor_nodal(neighbor_nodal)
{
}

NeighborCoupleable::~NeighborCoupleable() {}

const VariableValue &
NeighborCoupleable::coupledNeighborValue(const std::string & var_name, unsigned int comp)
{
  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  if (_neighbor_nodal)
    return (_c_is_implicit) ? var->nodalValueNeighbor() : var->nodalValueOldNeighbor();
  else
    return (_c_is_implicit) ? var->slnNeighbor() : var->slnOldNeighbor();
}

const VariableValue &
NeighborCoupleable::coupledNeighborValueOld(const std::string & var_name, unsigned int comp)
{
  validateExecutionerType(var_name);

  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  if (_neighbor_nodal)
    return (_c_is_implicit) ? var->nodalValueOldNeighbor() : var->nodalValueOlderNeighbor();
  else
    return (_c_is_implicit) ? var->slnOldNeighbor() : var->slnOlderNeighbor();
}

const VariableValue &
NeighborCoupleable::coupledNeighborValueOlder(const std::string & var_name, unsigned int comp)
{
  validateExecutionerType(var_name);

  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  if (_neighbor_nodal)
  {
    if (_c_is_implicit)
      return var->nodalValueOlderNeighbor();
    else
      mooseError("Older values not available for explicit schemes");
  }
  else
  {
    if (_c_is_implicit)
      return var->slnOlderNeighbor();
    else
      mooseError("Older values not available for explicit schemes");
  }
}

const VariableGradient &
NeighborCoupleable::coupledNeighborGradient(const std::string & var_name, unsigned int comp)
{
  if (_neighbor_nodal)
    mooseError("Nodal variables do not have gradients");

  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  return (_c_is_implicit) ? var->gradSlnNeighbor() : var->gradSlnOldNeighbor();
}

const VariableGradient &
NeighborCoupleable::coupledNeighborGradientOld(const std::string & var_name, unsigned int comp)
{
  if (_neighbor_nodal)
    mooseError("Nodal variables do not have gradients");

  validateExecutionerType(var_name);
  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  return (_c_is_implicit) ? var->gradSlnOldNeighbor() : var->gradSlnOlderNeighbor();
}

const VariableGradient &
NeighborCoupleable::coupledNeighborGradientOlder(const std::string & var_name, unsigned int comp)
{
  if (_neighbor_nodal)
    mooseError("Nodal variables do not have gradients");

  validateExecutionerType(var_name);
  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  if (_c_is_implicit)
    return var->gradSlnOlderNeighbor();
  else
    mooseError("Older values not available for explicit schemes");
}

const VariableSecond &
NeighborCoupleable::coupledNeighborSecond(const std::string & var_name, unsigned int comp)
{
  if (_neighbor_nodal)
    mooseError("Nodal variables do not have second derivatives");

  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  return (_c_is_implicit) ? var->secondSlnNeighbor() : var->secondSlnOldNeighbor();
}

const DenseVector<Number> &
NeighborCoupleable::coupledNeighborSolutionDoFs(const std::string & var_name, unsigned int comp)
{
  if (_neighbor_nodal)
    mooseError("nodal objects should not call coupledSolutionDoFs");

  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  return (_c_is_implicit) ? var->solutionDoFsNeighbor() : var->solutionDoFsOldNeighbor();
}

const DenseVector<Number> &
NeighborCoupleable::coupledNeighborSolutionDoFsOld(const std::string & var_name, unsigned int comp)
{
  if (_neighbor_nodal)
    mooseError("nodal objects should not call coupledSolutionDoFsOld");

  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  return (_c_is_implicit) ? var->solutionDoFsOldNeighbor() : var->solutionDoFsOlderNeighbor();
}

const DenseVector<Number> &
NeighborCoupleable::coupledNeighborSolutionDoFsOlder(const std::string & var_name,
                                                     unsigned int comp)
{
  if (_neighbor_nodal)
    mooseError("nodal objects should not call coupledSolutionDoFsOlder");

  MooseVariable * var = dynamic_cast<MooseVariable *>(getVar(var_name, comp));
  if (_c_is_implicit)
    return var->solutionDoFsOlderNeighbor();
  else
    mooseError("Older values not available for explicit schemes");
}
