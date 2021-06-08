//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Navier-Stokes includes
#include "GasMixPorousMixedVarMaterial.h"
#include "NS.h"

// FluidProperties includes
#include "SinglePhaseFluidProperties.h"

registerMooseObject("NavierStokesApp", GasMixPorousMixedVarMaterial);

InputParameters
GasMixPorousMixedVarMaterial::validParams()
{
  auto params = Material::validParams();
  params.addRequiredParam<UserObjectName>(NS::fluid, "fluid userobject");
  params.addRequiredCoupledVar("secondary_fraction", "Secondary mass fractions for gas species");
  params.addRequiredCoupledVar(NS::pressure, "The pressure");
  params.addRequiredCoupledVar(NS::T_fluid, "The fluid temperature");
  params.addRequiredCoupledVar(NS::superficial_momentum_x, "The x-momentum times the porosity");
  params.addCoupledVar(NS::superficial_momentum_y, "The y-momentum times the porosity");
  params.addCoupledVar(NS::superficial_momentum_z, "The z-momentum times the porosity");
  params.addClassDescription("Provides access to variables for a primitive variable set "
                             "of pressure, temperature, and superficial velocity");
  params.addRequiredParam<MaterialPropertyName>(NS::porosity, "the porosity");
  return params;
}

GasMixPorousMixedVarMaterial::GasMixPorousMixedVarMaterial(const InputParameters & params)
  : Material(params),
    _fluid(UserObjectInterface::getUserObject<SinglePhaseFluidProperties>(NS::fluid)),
    // mixed variables
    _var_fraction(adCoupledValue("secondary_fraction")),
    _grad_var_fraction(adCoupledGradient("secondary_fraction")),
    _fraction_dot(_is_transient ? adCoupledDot("secondary_fraction") : _ad_zero),
    _var_pressure(adCoupledValue(NS::pressure)),
    _grad_var_pressure(adCoupledGradient(NS::pressure)),
    _pressure_dot(_is_transient ? adCoupledDot(NS::pressure) : _ad_zero),
    _var_T_fluid(adCoupledValue(NS::T_fluid)),
    _grad_var_T_fluid(adCoupledGradient(NS::T_fluid)),
    _T_fluid_dot(_is_transient ? adCoupledDot(NS::T_fluid) : _ad_zero),
    _var_sup_mom_x(adCoupledValue(NS::superficial_momentum_x)),
    _grad_var_sup_mom_x(adCoupledGradient(NS::superficial_momentum_x)),
    _var_sup_mom_y(isCoupled(NS::superficial_momentum_y)
                       ? adCoupledValue(NS::superficial_momentum_y)
                       : _ad_zero),
    _grad_var_sup_mom_y(isCoupled(NS::superficial_momentum_y)
                            ? adCoupledGradient(NS::superficial_momentum_y)
                            : _ad_grad_zero),
    _var_sup_mom_z(isCoupled(NS::superficial_momentum_z)
                       ? adCoupledValue(NS::superficial_momentum_z)
                       : _ad_zero),
    _grad_var_sup_mom_z(isCoupled(NS::superficial_momentum_z)
                            ? adCoupledGradient(NS::superficial_momentum_z)
                            : _ad_grad_zero),
    _var_sup_mom_x_dot(_is_transient ? adCoupledDot(NS::superficial_momentum_x) : _ad_zero),
    _var_sup_mom_y_dot((isCoupled(NS::superficial_momentum_y) && _is_transient)
                           ? adCoupledDot(NS::superficial_momentum_y)
                           : _ad_zero),
    _var_sup_mom_z_dot((isCoupled(NS::superficial_momentum_z) && _is_transient)
                           ? adCoupledDot(NS::superficial_momentum_z)
                           : _ad_zero),
    // porosity
    _epsilon(getMaterialProperty<Real>(NS::porosity)),
    // properties: primitives
    _fraction(declareADProperty<Real>("fraction")),
    _grad_fraction(declareADProperty<RealVectorValue>(NS::grad("fraction"))),
    _pressure(declareADProperty<Real>(NS::pressure)),
    _grad_pressure(declareADProperty<RealVectorValue>(NS::grad(NS::pressure))),
    _T_fluid(declareADProperty<Real>(NS::T_fluid)),
    _grad_T_fluid(declareADProperty<RealVectorValue>(NS::grad(NS::T_fluid))),
    _sup_vel(declareADProperty<RealVectorValue>(NS::superficial_velocity)),
    _sup_vel_x(declareADProperty<Real>(NS::superficial_velocity_x)),
    _grad_sup_vel_x(declareADProperty<RealVectorValue>(NS::grad(NS::superficial_velocity_x))),
    _sup_vel_y(declareADProperty<Real>(NS::superficial_velocity_y)),
    _grad_sup_vel_y(declareADProperty<RealVectorValue>(NS::grad(NS::superficial_velocity_y))),
    _sup_vel_z(declareADProperty<Real>(NS::superficial_velocity_z)),
    _grad_sup_vel_z(declareADProperty<RealVectorValue>(NS::grad(NS::superficial_velocity_z))),
    // properties: for viz
    _rho(declareADProperty<Real>(NS::density)),
    _sup_rho_dot(declareADProperty<Real>(NS::time_deriv(NS::superficial_density))),
    _vel_x(declareADProperty<Real>(NS::velocity_x)),
    _vel_y(declareADProperty<Real>(NS::velocity_y)),
    _vel_z(declareADProperty<Real>(NS::velocity_z)),
    _sup_mom_x(declareADProperty<Real>(NS::superficial_momentum_x)),
    _sup_mom_y(declareADProperty<Real>(NS::superficial_momentum_y)),
    _sup_mom_z(declareADProperty<Real>(NS::superficial_momentum_z)),
    _sup_mom_x_dot(declareADProperty<Real>(NS::time_deriv(NS::superficial_momentum_x))),
    _sup_mom_y_dot(declareADProperty<Real>(NS::time_deriv(NS::superficial_momentum_y))),
    _sup_mom_z_dot(declareADProperty<Real>(NS::time_deriv(NS::superficial_momentum_z))),
    _sup_rho_et_dot(declareADProperty<Real>(NS::time_deriv(NS::superficial_total_energy_density))),
    _mom(declareADProperty<RealVectorValue>(NS::momentum)),
    _mom_x(declareADProperty<Real>(NS::momentum_x)),
    _mom_y(declareADProperty<Real>(NS::momentum_y)),
    _mom_z(declareADProperty<Real>(NS::momentum_z)),
    _speed(declareADProperty<Real>(NS::speed)),
    _rho_et(declareADProperty<Real>(NS::total_energy_density)),
    _sup_rho_f_dot(declareADProperty<Real>(NS::time_deriv("rho_f")))
{
  if (_mesh.dimension() >= 2 && !isCoupled(NS::superficial_momentum_y))
    mooseError("You must couple in a superficial y-momentum when solving 2D or 3D problems.");

  if (_mesh.dimension() >= 3 && !isCoupled(NS::superficial_momentum_z))
    mooseError("You must couple in a superficial z-momentum when solving 3D problems.");
}

void
GasMixPorousMixedVarMaterial::computeQpProperties()
{
  // Our primitive variable set
  _fraction[_qp] = _var_fraction[_qp];
  _grad_fraction[_qp] = _grad_var_fraction[_qp];
  _pressure[_qp] = _var_pressure[_qp];
  _grad_pressure[_qp] = _grad_var_pressure[_qp];
  _T_fluid[_qp] = _var_T_fluid[_qp];
  _grad_T_fluid[_qp] = _grad_var_T_fluid[_qp];
  const VectorValue<ADReal> superficial_momentum = {
      _var_sup_mom_x[_qp], _var_sup_mom_y[_qp], _var_sup_mom_z[_qp]};
  _sup_mom_x[_qp] = superficial_momentum(0);
  _sup_mom_y[_qp] = superficial_momentum(1);
  _sup_mom_z[_qp] = superficial_momentum(2);
  _sup_mom_x_dot[_qp] = _var_sup_mom_x_dot[_qp];
  _sup_mom_y_dot[_qp] = _var_sup_mom_y_dot[_qp];
  _sup_mom_z_dot[_qp] = _var_sup_mom_z_dot[_qp];

  // ----------------------------------------------------------------------------------------------------------------------------------------
  std::vector<ADReal> mass_fractions(1);
  mass_fractions[0] = _fraction[_qp];
  ADReal drho_dp, drho_dT, drho_dx;
  _fluid.rho_from_p_T_X(_pressure[_qp], _T_fluid[_qp], mass_fractions, _rho[_qp], drho_dp, drho_dT, drho_dx);
  const auto rho_dot = drho_dp * _pressure_dot[_qp] + drho_dT * _T_fluid_dot[_qp] + drho_dx * _fraction_dot[_qp];
  const auto grad_rho = drho_dp * _grad_pressure[_qp] + drho_dT * _grad_T_fluid[_qp] + drho_dx * _grad_fraction[_qp] ;
  _sup_rho_f_dot[_qp] = _epsilon[_qp] * (rho_dot * _fraction[_qp]  + _rho[_qp] * _fraction_dot[_qp]);
  // ----------------------------------------------------------------------------------------------------------------------------------------
  _sup_rho_dot[_qp] = _epsilon[_qp] * rho_dot;
  _sup_vel[_qp] = superficial_momentum / _rho[_qp];
  _sup_vel_x[_qp] = _sup_vel[_qp](0);
  _sup_vel_y[_qp] = _sup_vel[_qp](1);
  _sup_vel_z[_qp] = _sup_vel[_qp](2);
  _grad_sup_vel_x[_qp] = _grad_var_sup_mom_x[_qp] / _rho[_qp] -
                         superficial_momentum(0) / (_rho[_qp] * _rho[_qp]) * grad_rho;
  _grad_sup_vel_y[_qp] = _grad_var_sup_mom_y[_qp] / _rho[_qp] -
                         superficial_momentum(1) / (_rho[_qp] * _rho[_qp]) * grad_rho;
  _grad_sup_vel_z[_qp] = _grad_var_sup_mom_z[_qp] / _rho[_qp] -
                         superficial_momentum(2) / (_rho[_qp] * _rho[_qp]) * grad_rho;
  const auto sup_vel_x_dot = _var_sup_mom_x_dot[_qp] / _rho[_qp] -
                             superficial_momentum(0) / (_rho[_qp] * _rho[_qp]) * rho_dot;
  const auto sup_vel_y_dot = _var_sup_mom_y_dot[_qp] / _rho[_qp] -
                             superficial_momentum(1) / (_rho[_qp] * _rho[_qp]) * rho_dot;
  const auto sup_vel_z_dot = _var_sup_mom_z_dot[_qp] / _rho[_qp] -
                             superficial_momentum(2) / (_rho[_qp] * _rho[_qp]) * rho_dot;

  const auto velocity = _sup_vel[_qp] / _epsilon[_qp];
  _vel_x[_qp] = velocity(0);
  _vel_y[_qp] = velocity(1);
  _vel_z[_qp] = velocity(2);

  const auto v = 1. / _rho[_qp];
  const auto v_dot = -rho_dot / (_rho[_qp] * _rho[_qp]);
  // ----------------------------------------------------------------------------------------------------------------------------------------
  ADReal e, de_dp, de_dT, de_dx;
  _fluid.e_from_p_T_X(_pressure[_qp], _T_fluid[_qp], mass_fractions, e, de_dp, de_dT, de_dx);
  const auto e_dot = de_dp * _pressure_dot[_qp] + de_dT * _T_fluid_dot[_qp] + de_dx * _fraction_dot[_qp];
  // ----------------------------------------------------------------------------------------------------------------------------------------
  const auto et = e + velocity * velocity / 2.;
  const auto velocity_dot =
      VectorValue<ADReal>(sup_vel_x_dot, sup_vel_y_dot, sup_vel_z_dot) / _epsilon[_qp];
  const auto et_dot = e_dot + velocity * velocity_dot;
  _sup_rho_et_dot[_qp] = _epsilon[_qp] * (rho_dot * et + et_dot * _rho[_qp]);

  _mom_x[_qp] = _sup_mom_x[_qp] / _epsilon[_qp];
  _mom_y[_qp] = _sup_mom_y[_qp] / _epsilon[_qp];
  _mom_z[_qp] = _sup_mom_z[_qp] / _epsilon[_qp];
  _mom[_qp] = {_mom_x[_qp], _mom_y[_qp], _mom_z[_qp]};

  // if the velocity is zero, then the norm function call fails because AD tries to calculate the
  // derivatives which causes a divide by zero - because d/dx(sqrt(f(x))) = 1/2/sqrt(f(x))*df/dx.
  // So add a bit of noise to avoid this failure mode.
  if ((MooseUtils::absoluteFuzzyEqual(velocity(0), 0)) &&
      (MooseUtils::absoluteFuzzyEqual(velocity(1), 0)) &&
      (MooseUtils::absoluteFuzzyEqual(velocity(2), 0)))
    _speed[_qp] = 1e-42;
  else
    _speed[_qp] = velocity.norm();

  _rho_et[_qp] = _rho[_qp] * et;
}
