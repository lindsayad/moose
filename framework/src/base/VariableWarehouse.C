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

#include "VariableWarehouse.h"
#include "MooseVariableField.h"
#include "MooseVariableScalar.h"
#include "MooseTypes.h"

VariableWarehouse::VariableWarehouse() {}

VariableWarehouse::~VariableWarehouse()
{
  for (auto & var : _all_objects)
    delete var;
}

void
VariableWarehouse::add(const std::string & var_name, MooseVariableBase * var)
{
  _names.push_back(var_name);
  _var_name[var_name] = var;
  _all_objects.push_back(var);

  if (dynamic_cast<MooseVariableFE *>(var) != NULL)
  {
    _vars.push_back(dynamic_cast<MooseVariableFE *>(var));
    if (dynamic_cast<MooseVariable *>(var) != NULL)
      _regular_vars.push_back(dynamic_cast<MooseVariable *>(var));
    else if (dynamic_cast<MooseVariableVector *>(var) != NULL)
      _vector_vars.push_back(dynamic_cast<MooseVariableVector *>(var));
    else
      mooseError("Unknown variable class passed into VariableWarehouse. Attempt to hack us?");
  }
  else if (dynamic_cast<MooseVariableScalar *>(var) != NULL)
    _scalar_vars.push_back(dynamic_cast<MooseVariableScalar *>(var));
  else
    mooseError("Unknown variable class passed into VariableWarehouse. Attempt to hack us?");
}

void
VariableWarehouse::addBoundaryVar(BoundaryID bnd, MooseVariableFE * var)
{
  _boundary_vars[bnd].insert(var);
}

void
VariableWarehouse::addBoundaryVar(const std::set<BoundaryID> & boundary_ids, MooseVariableFE * var)
{
  for (const auto & bid : boundary_ids)
    addBoundaryVar(bid, var);
}

void
VariableWarehouse::addBoundaryVars(
    const std::set<BoundaryID> & boundary_ids,
    const std::map<std::string, std::vector<MooseVariableFE *>> & vars)
{
  for (const auto & bid : boundary_ids)
    for (const auto & it : vars)
      for (const auto & var : it.second)
        addBoundaryVar(bid, var);
}

MooseVariableBase *
VariableWarehouse::getVariable(const std::string & var_name)
{
  return _var_name[var_name];
}

MooseVariableBase *
VariableWarehouse::getVariable(unsigned int var_number)
{
  if (var_number < _all_objects.size())
    return _all_objects[var_number];
  else
    return NULL;
}

const std::vector<VariableName> &
VariableWarehouse::names() const
{
  return _names;
}

const std::vector<MooseVariableFE *> &
VariableWarehouse::variables()
{
  return _vars;
}

const std::vector<MooseVariable *> &
VariableWarehouse::regularVariables()
{
  return _regular_vars;
}

const std::vector<MooseVariableVector *> &
VariableWarehouse::vectorVariables()
{
  return _vector_vars;
}

const std::vector<MooseVariableScalar *> &
VariableWarehouse::scalars()
{
  return _scalar_vars;
}

const std::set<MooseVariableFE *> &
VariableWarehouse::boundaryVars(BoundaryID bnd)
{
  return _boundary_vars[bnd];
}
