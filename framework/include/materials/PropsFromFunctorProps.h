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
 * This material automatically declares as material properties whatever is passed to it
 * through the parameters 'prop_names' and uses the values from 'prop_values' as the values
 * for those properties.
 */
template <bool is_ad>
class PropsFromFunctorPropsTempl : public Material
{
public:
  static InputParameters validParams();

  PropsFromFunctorPropsTempl(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const std::vector<std::string> & _prop_names;

  unsigned int _num_props;

  std::vector<GenericMaterialProperty<Real, is_ad> *> _properties;
  std::vector<const FunctorMaterialProperty<GenericReal<is_ad>> *> _functor_properties;
};

typedef PropsFromFunctorPropsTempl<false> PropsFromFunctorProps;
typedef PropsFromFunctorPropsTempl<true> ADPropsFromFunctorProps;
