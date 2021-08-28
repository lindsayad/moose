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

/**
 * A material that optionally couples both a finite element and finite volume variable (strictly
 * speaking they don't have to be one or the either)
 */
class FEFVCouplingMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  FEFVCouplingMaterial(const InputParameters & parameters);

protected:
  const FunctorInterface<ADReal> & _fe_var;
  const FunctorInterface<ADReal> & _fv_var;
  FunctorMaterialProperty<ADReal> * const _fe_prop;
  FunctorMaterialProperty<ADReal> * const _fv_prop;
  FunctorMaterialProperty<ADReal> * const _declared_prop;
  const FunctorInterface<ADReal> * const _retrieved_prop;
};
