//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PropsFromFunctorProps.h"

registerMooseObject("MooseApp", PropsFromFunctorProps);
registerMooseObject("MooseApp", ADPropsFromFunctorProps);

template <bool is_ad>
InputParameters
PropsFromFunctorPropsTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Declares material properties based on names and values prescribed by input parameters.");
  params.addRequiredParam<std::vector<std::string>>(
      "prop_names", "The names of the properties this material will have");
  return params;
}

template <bool is_ad>
PropsFromFunctorPropsTempl<is_ad>::PropsFromFunctorPropsTempl(const InputParameters & parameters)
  : Material(parameters), _prop_names(getParam<std::vector<std::string>>("prop_names"))
{
  unsigned int num_names = _prop_names.size();

  _num_props = num_names;

  _properties.resize(num_names);
  _functor_properties.resize(num_names);

  for (unsigned int i = 0; i < _num_props; i++)
    for (const auto i : make_range(_num_props))
    {
      _properties[i] = &declareGenericProperty<Real, is_ad>(_prop_names[i]);
      _functor_properties[i] = &getFunctorMaterialProperty<GenericReal<is_ad>>(_prop_names[i]);
    }
}

template <bool is_ad>
void
PropsFromFunctorPropsTempl<is_ad>::computeQpProperties()
{
  for (const auto i : make_range(_num_props))
    (*_properties[i])[_qp] = (*_functor_properties[i])(_current_elem);
}

template class PropsFromFunctorPropsTempl<false>;
template class PropsFromFunctorPropsTempl<true>;
