//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVMomentumAdvectionPorosityGradient.h"
#include "PINSFVSuperficialVelocityVariable.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", PINSFVMomentumAdvectionPorosityGradient);

InputParameters
PINSFVMomentumAdvectionPorosityGradient::validParams()
{
  auto params = INSFVElementalKernel::validParams();
  params.addClassDescription(
      "Porosity gradient spun from the advection term for the porous media Navier Stokes "
      "momentum equation.");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity");

  params.addRequiredParam<MooseFunctorName>("u", "The superficial velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", 0, "The superficial velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", 0, "The superficial velocity in the z direction.");

  params.addRequiredParam<Real>(NS::density, "The value for the density");
  params.declareControllable(NS::density);

  params.addRequiredParam<bool>("smooth_porosity", "Whether the porosity has no discontinuities");
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

PINSFVMomentumAdvectionPorosityGradient::PINSFVMomentumAdvectionPorosityGradient(
    const InputParameters & params)
  : INSFVElementalKernel(params),
    _eps(getFunctor<ADReal>(NS::porosity)),
    _u(getFunctor<ADReal>("u")),
    _v(getFunctor<ADReal>("v")),
    _w(getFunctor<ADReal>("w")),
    _rho(getParam<Real>(NS::density)),
    _smooth_porosity(getParam<bool>("smooth_porosity"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("PINSFV is not supported by local AD indexing. In order to use PINSFV, please run "
             "the configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif
  if (!dynamic_cast<PINSFVSuperficialVelocityVariable *>(&_var))
    mooseError("PINSFVMomentumAdvectionPorosityGradient may only be used with a superficial "
               "velocity variable, of variable type PINSFVSuperficialVelocityVariable.");

  if (!_smooth_porosity)
    paramError(
        "smooth_porosity",
        "The MomentumAdvectionContinuousPorosity may only be used with a continuous porosity.");
}

void
PINSFVMomentumAdvectionPorosityGradient::gatherRCData(const Elem & elem)
{
  const Real one_over_eps = 1 / MetaPhysicL::raw_value(_eps(&elem));
  ADRealVectorValue V = {_u(&elem), _v(&elem), _w(&elem)};

  _rc_uo.addToB(&elem,
                _index,
                _rho * V(_index) * (-one_over_eps * one_over_eps) *
                    (V * MetaPhysicL::raw_value(_eps.gradient(&elem))));
}
