//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorIsotropicDragCoefficients.h"
#include "NS.h"

/**
 * Abstract base class to compute isotropic drag coefficients in a
 * pebble bed. These correlations generally share common proportionalities to
 * porosity, hydraulic diameter, etc. that are placed into this
 * base class. The Darcy coefficient in this case scales a term with
 * proportionality \f$\frac{\mu_f}{\rho_f}\frac{\epsilon}{D_h^2}\f$,
 * while the Forchheimer coefficient scales a term with proportionality
 * \f$\epsilon\|\vec{V}\|/D_h\f$.
 */
template <typename Derived>
class FunctorPebbleBedDragCoefficients : public FunctorIsotropicDragCoefficients<Derived>
{
  friend class FunctorIsotropicDragCoefficients<Derived>;

public:
  FunctorPebbleBedDragCoefficients(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  template <typename Space, typename Time>
  ADReal computeDarcyPrefactor(const Space & r, const Time & t);

  template <typename Space, typename Time>
  ADReal computeForchheimerPrefactor(const Space & r, const Time & t);

  /// Compute hydraulic diameter in the bed
  template <typename Space, typename Time>
  ADReal computeHydraulicDiameter(const Space & r, const Time & t);

  /// porosity
  const Moose::Functor<ADReal> & _eps;

  /// fluid density
  const Moose::Functor<ADReal> & _rho;

  /// fluid dynamic viscosity
  const Moose::Functor<ADReal> & _mu;

  /// pebble diameter
  const Real & _d_pebble;
};

template <typename Derived>
InputParameters
FunctorPebbleBedDragCoefficients<Derived>::validParams()
{
  InputParameters params = FunctorIsotropicDragCoefficients<Derived>::validParams();
  params.addRequiredRangeCheckedParam<Real>(
      NS::pebble_diameter, NS::pebble_diameter + " > 0.0", "Pebble diameter");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity");

  params.addParam<MooseFunctorName>(NS::density, NS::density, "Density");
  params.addParam<MooseFunctorName>(NS::mu, NS::mu, "Dynamic viscosity");

  return params;
}

template <typename Derived>
FunctorPebbleBedDragCoefficients<Derived>::FunctorPebbleBedDragCoefficients(
    const InputParameters & parameters)
  : FunctorIsotropicDragCoefficients<Derived>(parameters),
    _eps(this->template getFunctor<ADReal>(NS::porosity)),
    _rho(this->template getFunctor<ADReal>(NS::density)),
    _mu(this->template getFunctor<ADReal>(NS::mu)),
    _d_pebble(this->template getParam<Real>(NS::pebble_diameter))
{
}

template <typename Derived>
template <typename Space, typename Time>
ADReal
FunctorPebbleBedDragCoefficients<Derived>::computeDarcyPrefactor(const Space & r, const Time & t)
{
  ADReal Dh = computeHydraulicDiameter(r, t);
  return (_mu(r, t) / _rho(r, t)) / (Dh * Dh);
}

template <typename Derived>
template <typename Space, typename Time>
ADReal
FunctorPebbleBedDragCoefficients<Derived>::computeForchheimerPrefactor(const Space & r,
                                                                       const Time & t)
{
  return 1 / computeHydraulicDiameter(r, t);
}

template <typename Derived>
template <typename Space, typename Time>
ADReal
FunctorPebbleBedDragCoefficients<Derived>::computeHydraulicDiameter(const Space & r, const Time & t)
{
  mooseAssert(1.0 - _eps(r, t).value() > 1e-8,
              "Bed hydraulic diameter is ill-defined at porosity of 1.");
  return _eps(r, t) * _d_pebble / (1.0 - _eps(r, t));
}
