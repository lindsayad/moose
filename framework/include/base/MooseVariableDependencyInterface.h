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

#ifndef MOOSEVARIABLEDEPENDENCYINTERFACE_H
#define MOOSEVARIABLEDEPENDENCYINTERFACE_H

#include <set>

// Forward declarations
class MooseVariableFE;

class MooseVariableDependencyInterface
{
public:
  MooseVariableDependencyInterface() {}

  /**
   * Retrieve the set of MooseVariableFEs that _this_ object depends on.
   * @return The MooseVariableFEs that MUST be reinited before evaluating this object
   */
  const std::set<MooseVariableFE *> & getMooseVariableDependencies() const
  {
    return _moose_variable_dependencies;
  }

protected:
  /**
   * Call this function to add the passed in MooseVariableFE as a variable that _this_ object
   * depends
   * on.
   */
  void addMooseVariableDependency(MooseVariableFE * var)
  {
    _moose_variable_dependencies.insert(var);
  }
  void addMooseVariableDependency(std::vector<MooseVariableFE *> vars)
  {
    _moose_variable_dependencies.insert(vars.begin(), vars.end());
  }

private:
  std::set<MooseVariableFE *> _moose_variable_dependencies;
};

#endif // MOOSEVARIABLEDEPENDENCYINTERFACE_H
