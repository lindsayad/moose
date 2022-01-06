//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * FunctorMaterials compute functor material properties
 */
class FunctorMaterial : public Material
{
public:
  static InputParameters validParams();

  FunctorMaterial(const InputParameters & parameters);
  void computeProperties() override final {}

protected:
  void computeQpProperties() override final {}

  /**
   * Declare a functor material property
   */
  template <typename T, typename PolymorphicLambda>
  const Moose::Functor<T> &
  addFunctorProperty(const std::string & name,
                     PolymorphicLambda my_lammy,
                     const std::set<ExecFlagType> & clearance_schedule = {EXEC_ALWAYS},
                     bool nonorthogonal_correction = true);

  /**
   * Declare a functor material property with specified subdomain ids
   */
  template <typename T, typename PolymorphicLambda>
  const Moose::Functor<T> &
  addFunctorPropertyByBlocks(const std::string & name,
                             PolymorphicLambda my_lammy,
                             const std::set<SubdomainID> & sub_ids,
                             const std::set<ExecFlagType> & clearance_schedule = {EXEC_ALWAYS},
                             bool nonorthogonal_correction = true);
};

template <typename T, typename PolymorphicLambda>
const Moose::Functor<T> &
FunctorMaterial::addFunctorProperty(const std::string & name,
                                    PolymorphicLambda my_lammy,
                                    const std::set<ExecFlagType> & clearance_schedule,
                                    const bool nonorthogonal_correction)
{
  return addFunctorPropertyByBlocks<T>(
      name, my_lammy, blockIDs(), clearance_schedule, nonorthogonal_correction);
}

template <typename T, typename PolymorphicLambda>
const Moose::Functor<T> &
FunctorMaterial::addFunctorPropertyByBlocks(const std::string & name,
                                            PolymorphicLambda my_lammy,
                                            const std::set<SubdomainID> & sub_ids,
                                            const std::set<ExecFlagType> & clearance_schedule,
                                            const bool nonorthogonal_correction)
{
  // Check if the supplied parameter is a valid input parameter key
  std::string prop_name = name;
  if (_pars.have_parameter<MaterialPropertyName>(name))
    prop_name = _pars.get<MaterialPropertyName>(name);

  return _subproblem.addPiecewiseByBlockLambdaFunctor<T>(
      prop_name, my_lammy, clearance_schedule, nonorthogonal_correction, _mesh, sub_ids, _tid);
}
