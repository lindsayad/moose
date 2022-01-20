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

class INSFVMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  INSFVMaterial(const InputParameters & parameters);

protected:
  /// density
  const Moose::Functor<ADReal> & _rho;

  const bool _has_temperature;

  const MooseVariableFVReal * const _temperature;
  const Moose::Functor<ADReal> * const _cp;
};
