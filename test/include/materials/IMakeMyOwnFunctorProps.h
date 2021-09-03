//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

template <typename>
class CustomFunctorProp;

/**
 * A material that optionally couples both a finite element and finite volume variable (strictly
 * speaking they don't have to be one or the either). This class also uses a specialization of
 * FunctorMaterialProperty that overrides the \p evaluate methods and doesn't use lambdas
 */
class IMakeMyOwnFunctorProps : public FunctorMaterial
{
public:
  static InputParameters validParams();

  IMakeMyOwnFunctorProps(const InputParameters & parameters);

protected:
  const FunctorInterface<ADReal> & _fe_var;
  const FunctorInterface<ADReal> & _fv_var;
  const FunctorInterface<ADReal> * const _retrieved_prop;
  CustomFunctorProp<ADReal> * const _fe_prop;
  CustomFunctorProp<ADReal> * const _fv_prop;
};
