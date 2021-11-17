//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVElementalKernel.h"

/**
 * An elemental kernel to add the inverse porosity gradient term to the momentum equation
 */
class PINSFVMomentumAdvectionPorosityGradient : public INSFVElementalKernel
{
public:
  static InputParameters validParams();
  PINSFVMomentumAdvectionPorosityGradient(const InputParameters & params);

  using INSFVElementalKernel::gatherRCData;
  void gatherRCData(const Elem &) override;

protected:
  /// porosity functor to compute gradients
  const Moose::Functor<ADReal> & _eps;
  /// superficial velocity x-component
  const Moose::Functor<ADReal> & _u;
  /// superficial velocity y-component
  const Moose::Functor<ADReal> & _v;
  /// superficial velocity z-component
  const Moose::Functor<ADReal> & _w;
  /// density
  const Real & _rho;
  /// whether the porosity has discontinuities, in which case this kernel should not be used
  const bool _smooth_porosity;
};
