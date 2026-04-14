//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "MFEMObject.h"
#include "SetupInterface.h"
#include "DependencyResolverInterface.h"

#include <algorithm>
#include <optional>
#include <set>
#include <type_traits>

/**
 * Base class for MFEM objects that participate in execution ordering but are not UserObjects.
 */
class MFEMExecutedObject : public MFEMObject,
                           public SetupInterface,
                           public DependencyResolverInterface
{
public:
  /**
   * Metadata describing one dependency-bearing input parameter on this object.
   */
  struct DependencyParam
  {
    /// The resource category requested through this parameter.
    std::string kind;
    /// The parameter name storing the dependency value(s).
    std::string param_name;
  };

  /**
   * Declare the common parameters used by MFEM executed objects.
   */
  static InputParameters validParams();

  /**
   * Construct an executed MFEM object and materialize its dependency metadata.
   */
  MFEMExecutedObject(const InputParameters & parameters);

  /**
   * Perform any pre-execution setup for this object.
   */
  virtual void initialize() {}
  /**
   * Perform the main work for this object.
   */
  virtual void execute() {}
  /**
   * Perform any post-execution finalization for this object.
   */
  virtual void finalize() {}

  /**
   * Return the variable names produced by this object.
   */
  virtual std::set<std::string> producedVariableNames() const;
  /**
   * Return the postprocessor names produced by this object.
   */
  virtual std::set<std::string> producedPostprocessorNames() const;
  /**
   * Return the vector postprocessor names produced by this object.
   */
  virtual std::set<std::string> producedVectorPostprocessorNames() const;

  virtual const std::set<std::string> & getRequestedItems() override;
  virtual const std::set<std::string> & getSuppliedItems() override;

  /**
   * Add an optional dependency-bearing parameter and register it with the MFEM scheduler.
   */
  template <typename T>
  static void addDependencyParam(InputParameters & params,
                                 const std::string & param_name,
                                 const std::string & doc_string);

  /**
   * Add a required dependency-bearing parameter and register it with the MFEM scheduler.
   */
  template <typename T>
  static void addRequiredDependencyParam(InputParameters & params,
                                         const std::string & param_name,
                                         const std::string & doc_string);

protected:
  /**
   * Build the dependency key used for a produced/requested variable.
   */
  static std::string variableDependencyKey(const std::string & name);
  /**
   * Build the dependency key used for a produced/requested postprocessor.
   */
  static std::string postprocessorDependencyKey(const std::string & name);
  /**
   * Build the dependency key used for a produced/requested vector postprocessor.
   */
  static std::string vectorPostprocessorDependencyKey(const std::string & name);

  /**
   * Record one dependency-bearing parameter in the private parameter metadata.
   */
  static void appendDependencyParam(InputParameters & params,
                                    const std::string & param_name,
                                    const std::string & kind);

private:
  /// Lazily constructed requested dependency keys for this object's registered dependencies.
  std::optional<std::set<std::string>> _requested_items;
  /// Lazily constructed supplied dependency keys for this object's produced resources.
  std::optional<std::set<std::string>> _supplied_items;
  /// Compact typed metadata describing dependency-bearing parameters on this object.
  std::vector<DependencyParam> _dependency_params;
};

template <typename T>
void
MFEMExecutedObject::addDependencyParam(InputParameters & params,
                                       const std::string & param_name,
                                       const std::string & doc_string)
{
  params.addParam<T>(param_name, doc_string);
  appendDependencyParam(params, param_name, typeid(T).name());
}

template <typename T>
void
MFEMExecutedObject::addRequiredDependencyParam(InputParameters & params,
                                               const std::string & param_name,
                                               const std::string & doc_string)
{
  params.addRequiredParam<T>(param_name, doc_string);
  appendDependencyParam(params, param_name, typeid(T).name());
}

#endif
