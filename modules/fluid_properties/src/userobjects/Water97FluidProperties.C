//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Water97FluidProperties.h"
#include "MathUtils.h"
#include "libmesh/utility.h"

registerMooseObject("FluidPropertiesApp", Water97FluidProperties);

InputParameters
Water97FluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addClassDescription("Fluid properties for water and steam (H2O) using IAPWS-IF97");
  return params;
}

Water97FluidProperties::Water97FluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    _Mh2o(18.015e-3),
    _Rw(461.526),
    _p_critical(22.064e6),
    _T_critical(647.096),
    _rho_critical(322.0),
    _p_triple(611.657),
    _T_triple(273.16)
{
}

Water97FluidProperties::~Water97FluidProperties() {}

std::string
Water97FluidProperties::fluidName() const
{
  return "water";
}

Real
Water97FluidProperties::molarMass() const
{
  return _Mh2o;
}

Real
Water97FluidProperties::criticalPressure() const
{
  return _p_critical;
}

Real
Water97FluidProperties::criticalTemperature() const
{
  return _T_critical;
}

Real
Water97FluidProperties::criticalDensity() const
{
  return _rho_critical;
}

Real
Water97FluidProperties::triplePointPressure() const
{
  return _p_triple;
}

Real
Water97FluidProperties::triplePointTemperature() const
{
  return _T_triple;
}

template <typename T>
T
Water97FluidProperties::dgamma1_dpi(const T & pi, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 0; i < _n1.size(); ++i)
    sum += -_n1[i] * _I1[i] * MathUtils::pow(7.1 - pi, _I1[i] - 1) *
           MathUtils::pow(tau - 1.222, _J1[i]);

  return sum;
}

template <typename T>
T
Water97FluidProperties::dgamma2_dpi(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 1.0 / pi;

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n2.size(); ++i)
    dgr += _n2[i] * _I2[i] * MathUtils::pow(pi, _I2[i] - 1) * MathUtils::pow(tau - 0.5, _J2[i]);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::gamma1(const T & pi, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 0; i < _n1.size(); ++i)
    sum += _n1[i] * MathUtils::pow(7.1 - pi, _I1[i]) * MathUtils::pow(tau - 1.222, _J1[i]);

  return sum;
}

template <typename T>
T
Water97FluidProperties::d2gamma1_dpi2(const T & pi, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 0; i < _n1.size(); ++i)
    sum += _n1[i] * _I1[i] * (_I1[i] - 1) * MathUtils::pow(7.1 - pi, _I1[i] - 2) *
           MathUtils::pow(tau - 1.222, _J1[i]);

  return sum;
}

template <typename T>
T
Water97FluidProperties::dgamma1_dtau(const T & pi, const T & tau) const
{
  T g = 0.0;
  for (std::size_t i = 0; i < _n1.size(); ++i)
    g += _n1[i] * _J1[i] * MathUtils::pow(7.1 - pi, _I1[i]) *
         MathUtils::pow(tau - 1.222, _J1[i] - 1);

  return g;
}

template <typename T>
T
Water97FluidProperties::d2gamma1_dtau2(const T & pi, const T & tau) const
{
  T dg = 0.0;
  for (std::size_t i = 0; i < _n1.size(); ++i)
    dg += _n1[i] * _J1[i] * (_J1[i] - 1) * MathUtils::pow(7.1 - pi, _I1[i]) *
          MathUtils::pow(tau - 1.222, _J1[i] - 2);

  return dg;
}

template <typename T>
T
Water97FluidProperties::d2gamma1_dpitau(const T & pi, const T & tau) const
{
  T dg = 0.0;
  for (std::size_t i = 0; i < _n1.size(); ++i)
    dg += -_n1[i] * _I1[i] * _J1[i] * MathUtils::pow(7.1 - pi, _I1[i] - 1) *
          MathUtils::pow(tau - 1.222, _J1[i] - 1);

  return dg;
}

template <typename T>
T
Water97FluidProperties::gamma2(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T sum0 = 0.0;
  for (std::size_t i = 0; i < _n02.size(); ++i)
    sum0 += _n02[i] * MathUtils::pow(tau, _J02[i]);

  T g0 = std::log(pi) + sum0;

  // Residual part of the Gibbs free energy
  T gr = 0.0;
  for (std::size_t i = 0; i < _n2.size(); ++i)
    gr += _n2[i] * MathUtils::pow(pi, _I2[i]) * MathUtils::pow(tau - 0.5, _J2[i]);

  return g0 + gr;
}

template <typename T>
T
Water97FluidProperties::d2gamma2_dpi2(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = -1.0 / pi / pi;

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n2.size(); ++i)
    dgr += _n2[i] * _I2[i] * (_I2[i] - 1) * MathUtils::pow(pi, _I2[i] - 2) *
           MathUtils::pow(tau - 0.5, _J2[i]);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::dgamma2_dtau(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 0.0;
  for (std::size_t i = 0; i < _n02.size(); ++i)
    dg0 += _n02[i] * _J02[i] * MathUtils::pow(tau, _J02[i] - 1);

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n2.size(); ++i)
    dgr += _n2[i] * _J2[i] * MathUtils::pow(pi, _I2[i]) * MathUtils::pow(tau - 0.5, _J2[i] - 1);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::d2gamma2_dtau2(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 0.0;
  for (std::size_t i = 0; i < _n02.size(); ++i)
    dg0 += _n02[i] * _J02[i] * (_J02[i] - 1) * MathUtils::pow(tau, _J02[i] - 2);

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n2.size(); ++i)
    dgr += _n2[i] * _J2[i] * (_J2[i] - 1) * MathUtils::pow(pi, _I2[i]) *
           MathUtils::pow(tau - 0.5, _J2[i] - 2);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::d2gamma2_dpitau(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 0.0;

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n2.size(); ++i)
    dgr += _n2[i] * _I2[i] * _J2[i] * MathUtils::pow(pi, _I2[i] - 1) *
           MathUtils::pow(tau - 0.5, _J2[i] - 1);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::phi3(const T & delta, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 1; i < _n3.size(); ++i)
    sum += _n3[i] * MathUtils::pow(delta, _I3[i]) * MathUtils::pow(tau, _J3[i]);

  return _n3[0] * std::log(delta) + sum;
}

template <typename T>
T
Water97FluidProperties::dphi3_ddelta(const T & delta, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 1; i < _n3.size(); ++i)
    sum += _n3[i] * _I3[i] * MathUtils::pow(delta, _I3[i] - 1) * MathUtils::pow(tau, _J3[i]);

  return _n3[0] / delta + sum;
}

template <typename T>
T
Water97FluidProperties::d2phi3_ddelta2(const T & delta, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 1; i < _n3.size(); ++i)
    sum += _n3[i] * _I3[i] * (_I3[i] - 1) * MathUtils::pow(delta, _I3[i] - 2) *
           MathUtils::pow(tau, _J3[i]);

  return -_n3[0] / delta / delta + sum;
}

template <typename T>
T
Water97FluidProperties::dphi3_dtau(const T & delta, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 1; i < _n3.size(); ++i)
    sum += _n3[i] * _J3[i] * MathUtils::pow(delta, _I3[i]) * MathUtils::pow(tau, _J3[i] - 1);

  return sum;
}

template <typename T>
T
Water97FluidProperties::d2phi3_dtau2(const T & delta, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 1; i < _n3.size(); ++i)
    sum += _n3[i] * _J3[i] * (_J3[i] - 1) * MathUtils::pow(delta, _I3[i]) *
           MathUtils::pow(tau, _J3[i] - 2);

  return sum;
}

template <typename T>
T
Water97FluidProperties::d2phi3_ddeltatau(const T & delta, const T & tau) const
{
  T sum = 0.0;
  for (std::size_t i = 1; i < _n3.size(); ++i)
    sum += _n3[i] * _I3[i] * _J3[i] * MathUtils::pow(delta, _I3[i] - 1) *
           MathUtils::pow(tau, _J3[i] - 1);

  return sum;
}

template <typename T>
T
Water97FluidProperties::gamma5(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T sum0 = 0.0;
  for (std::size_t i = 0; i < _n05.size(); ++i)
    sum0 += _n05[i] * MathUtils::pow(tau, _J05[i]);

  T g0 = std::log(pi) + sum0;

  // Residual part of the Gibbs free energy
  T gr = 0.0;
  for (std::size_t i = 0; i < _n5.size(); ++i)
    gr += _n5[i] * MathUtils::pow(pi, _I5[i]) * MathUtils::pow(tau, _J5[i]);

  return g0 + gr;
}

template <typename T>
T
Water97FluidProperties::dgamma5_dpi(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 1.0 / pi;

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n5.size(); ++i)
    dgr += _n5[i] * _I5[i] * MathUtils::pow(pi, _I5[i] - 1) * MathUtils::pow(tau, _J5[i]);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::d2gamma5_dpi2(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = -1.0 / pi / pi;

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n5.size(); ++i)
    dgr += _n5[i] * _I5[i] * (_I5[i] - 1) * MathUtils::pow(pi, _I5[i] - 2) *
           MathUtils::pow(tau, _J5[i]);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::dgamma5_dtau(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 0.0;
  for (std::size_t i = 0; i < _n05.size(); ++i)
    dg0 += _n05[i] * _J05[i] * MathUtils::pow(tau, _J05[i] - 1);

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n5.size(); ++i)
    dgr += _n5[i] * _J5[i] * MathUtils::pow(pi, _I5[i]) * MathUtils::pow(tau, _J5[i] - 1);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::d2gamma5_dtau2(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 0.0;
  for (std::size_t i = 0; i < _n05.size(); ++i)
    dg0 += _n05[i] * _J05[i] * (_J05[i] - 1) * MathUtils::pow(tau, _J05[i] - 2);

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n5.size(); ++i)
    dgr += _n5[i] * _J5[i] * (_J5[i] - 1) * MathUtils::pow(pi, _I5[i]) *
           MathUtils::pow(tau, _J5[i] - 2);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::d2gamma5_dpitau(const T & pi, const T & tau) const
{
  // Ideal gas part of the Gibbs free energy
  T dg0 = 0.0;

  // Residual part of the Gibbs free energy
  T dgr = 0.0;
  for (std::size_t i = 0; i < _n5.size(); ++i)
    dgr +=
        _n5[i] * _I5[i] * _J5[i] * MathUtils::pow(pi, _I5[i] - 1) * MathUtils::pow(tau, _J5[i] - 1);

  return dg0 + dgr;
}

template <typename T>
T
Water97FluidProperties::tempXY(const T & pressure, const subregionEnum xy) const
{
  T pi = pressure / 1.0e6;

  // Choose the constants based on the string xy
  unsigned int row;

  switch (xy)
  {
    case AB:
      row = 0;
      break;
    case CD:
      row = 1;
      break;
    case GH:
      row = 2;
      break;
    case IJ:
      row = 3;
      break;
    case JK:
      row = 4;
      break;
    case MN:
      row = 5;
      break;
    case OP:
      row = 6;
      break;
    case QU:
      row = 7;
      break;
    case RX:
      row = 8;
      break;
    case UV:
      row = 9;
      break;
    case WX:
      row = 10;
      break;
    default:
      row = 0;
  }

  T sum = 0.0;

  if (xy == AB || xy == OP || xy == WX)
    for (std::size_t i = 0; i < _tempXY_n[row].size(); ++i)
      sum += _tempXY_n[row][i] * MathUtils::pow(std::log(pi), _tempXY_I[row][i]);
  else if (xy == EF)
    sum += 3.727888004 * (pi - _p_critical / 1.0e6) + _T_critical;
  else
    for (std::size_t i = 0; i < _tempXY_n[row].size(); ++i)
      sum += _tempXY_n[row][i] * MathUtils::pow(pi, _tempXY_I[row][i]);

  return sum;
}

template <typename T>
void
Water97FluidProperties::vaporPressureTemplate(const T & temperature, T & psat, T & dpsat_dT) const
{
  // Check whether the input temperature is within the region of validity of this equation.
  // Valid for 273.15 K <= t <= 647.096 K
  if (temperature < 273.15 || temperature > _T_critical)
    mooseException(name(),
                   ": vaporPressure(): Temperature is outside range 273.15 K <= T <= 647.096 K");

  T theta, dtheta_dT, theta2, a, b, c, da_dtheta, db_dtheta, dc_dtheta;
  theta = temperature + _n4[8] / (temperature - _n4[9]);
  dtheta_dT = 1.0 - _n4[8] / (temperature - _n4[9]) / (temperature - _n4[9]);
  theta2 = theta * theta;

  a = theta2 + _n4[0] * theta + _n4[1];
  b = _n4[2] * theta2 + _n4[3] * theta + _n4[4];
  c = _n4[5] * theta2 + _n4[6] * theta + _n4[7];

  da_dtheta = 2.0 * theta + _n4[0];
  db_dtheta = 2.0 * _n4[2] * theta + _n4[3];
  dc_dtheta = 2.0 * _n4[5] * theta + _n4[6];

  T denominator = -b + std::sqrt(b * b - 4.0 * a * c);

  psat = Utility::pow<4>(2.0 * c / denominator) * 1.0e6;

  // The derivative wrt temperature is given by the chain rule
  T dpsat = 4.0 * Utility::pow<3>(2.0 * c / denominator);
  dpsat *= (2.0 * dc_dtheta / denominator -
            2.0 * c / denominator / denominator *
                (-db_dtheta + std::pow(b * b - 4.0 * a * c, -0.5) *
                                  (b * db_dtheta - 2.0 * da_dtheta * c - 2.0 * a * dc_dtheta)));
  dpsat_dT = dpsat * dtheta_dT * 1.0e6;
}

inline void
Water97FluidProperties::vaporPressure(const Real temperature, Real & psat, Real & dpsat_dT) const
{
  vaporPressureTemplate(temperature, psat, dpsat_dT);
}

template <typename T>
T
Water97FluidProperties::v_from_p_T_template(const T & pressure, const T & temperature) const
{
  return 1 / rho_from_p_T_template(pressure, temperature);
}

template <typename T>
void
Water97FluidProperties::v_from_p_T_template(
    const T & pressure, const T & temperature, T & v, T & dv_dp, T & dv_dT) const
{
  auto functor = [this](const auto & pressure, const auto & temperature)
  { return this->v_from_p_T_template(pressure, temperature); };

  xyDerivatives(pressure, temperature, v, dv_dp, dv_dT, functor);
}

template <typename T>
T
Water97FluidProperties::e_from_p_T_template(const T & pressure, const T & temperature) const
{
  T internal_energy, pi, tau;

  // Determine which region the point is in
  unsigned int region =
      inRegion(MetaPhysicL::raw_value(pressure), MetaPhysicL::raw_value(temperature));
  switch (region)
  {
    case 1:
      pi = pressure / _p_star[0];
      tau = _T_star[0] / temperature;
      internal_energy =
          _Rw * temperature * (tau * dgamma1_dtau(pi, tau) - pi * dgamma1_dpi(pi, tau));
      break;

    case 2:
      pi = pressure / _p_star[1];
      tau = _T_star[1] / temperature;
      internal_energy =
          _Rw * temperature * (tau * dgamma2_dtau(pi, tau) - pi * dgamma2_dpi(pi, tau));
      break;

    case 3:
    {
      // Calculate density first, then use that in Helmholtz free energy
      T density3 = densityRegion3(pressure, temperature);
      T delta = density3 / _rho_critical;
      tau = _T_star[2] / temperature;
      internal_energy = _Rw * temperature * tau * dphi3_dtau(delta, tau);
      break;
    }

    case 5:
      pi = pressure / _p_star[4];
      tau = _T_star[4] / temperature;
      internal_energy =
          _Rw * temperature * (tau * dgamma5_dtau(pi, tau) - pi * dgamma5_dpi(pi, tau));
      break;

    default:
      mooseError("inRegion() has given an incorrect region");
  }
  // Output in J/kg
  return internal_energy;
}

template <typename T>
void
Water97FluidProperties::e_from_p_T_template(
    const T & pressure, const T & temperature, T & e, T & de_dp, T & de_dT) const
{
  T pi, tau, dinternal_energy_dp, dinternal_energy_dT;

  // Determine which region the point is in
  unsigned int region =
      inRegion(MetaPhysicL::raw_value(pressure), MetaPhysicL::raw_value(temperature));
  switch (region)
  {
    case 1:
    {
      pi = pressure / _p_star[0];
      tau = _T_star[0] / temperature;
      T dgdp = dgamma1_dpi(pi, tau);
      T d2gdpt = d2gamma1_dpitau(pi, tau);
      dinternal_energy_dp =
          _Rw * temperature * (tau * d2gdpt - dgdp - pi * d2gamma1_dpi2(pi, tau)) / _p_star[0];
      dinternal_energy_dT =
          _Rw * (pi * tau * d2gdpt - tau * tau * d2gamma1_dtau2(pi, tau) - pi * dgdp);
      break;
    }

    case 2:
    {
      pi = pressure / _p_star[1];
      tau = _T_star[1] / temperature;
      T dgdp = dgamma2_dpi(pi, tau);
      T d2gdpt = d2gamma2_dpitau(pi, tau);
      dinternal_energy_dp =
          _Rw * temperature * (tau * d2gdpt - dgdp - pi * d2gamma2_dpi2(pi, tau)) / _p_star[1];
      dinternal_energy_dT =
          _Rw * (pi * tau * d2gdpt - tau * tau * d2gamma2_dtau2(pi, tau) - pi * dgdp);
      break;
    }

    case 3:
    {
      // Calculate density first, then use that in Helmholtz free energy
      T density3 = densityRegion3(pressure, temperature);
      T delta = density3 / _rho_critical;
      tau = _T_star[2] / temperature;
      T dpdd = dphi3_ddelta(delta, tau);
      T d2pddt = d2phi3_ddeltatau(delta, tau);
      T d2pdd2 = d2phi3_ddelta2(delta, tau);
      dinternal_energy_dp =
          _T_star[2] * d2pddt / _rho_critical /
          (2.0 * temperature * delta * dpdd + temperature * delta * delta * d2pdd2);
      dinternal_energy_dT =
          -_Rw * (delta * tau * d2pddt * (dpdd - tau * d2pddt) / (2.0 * dpdd + delta * d2pdd2) +
                  tau * tau * d2phi3_dtau2(delta, tau));
      break;
    }

    case 5:
    {
      pi = pressure / _p_star[4];
      tau = _T_star[4] / temperature;
      T dgdp = dgamma5_dpi(pi, tau);
      T d2gdpt = d2gamma5_dpitau(pi, tau);
      dinternal_energy_dp =
          _Rw * temperature * (tau * d2gdpt - dgdp - pi * d2gamma5_dpi2(pi, tau)) / _p_star[4];
      dinternal_energy_dT =
          _Rw * (pi * tau * d2gdpt - tau * tau * d2gamma5_dtau2(pi, tau) - pi * dgdp);
      break;
    }

    default:
      mooseError("inRegion has given an incorrect region");
  }

  e = this->e_from_p_T(pressure, temperature);
  de_dp = dinternal_energy_dp;
  de_dT = dinternal_energy_dT;
}

template <typename T>
T
Water97FluidProperties::c_from_p_T_template(const T & pressure, const T & temperature) const
{
  T speed2, pi, tau, delta;

  // Determine which region the point is in
  unsigned int region =
      inRegion(MetaPhysicL::raw_value(pressure), MetaPhysicL::raw_value(temperature));
  switch (region)
  {
    case 1:
      pi = pressure / _p_star[0];
      tau = _T_star[0] / temperature;
      speed2 = _Rw * temperature * Utility::pow<2>(dgamma1_dpi(pi, tau)) /
               (Utility::pow<2>(dgamma1_dpi(pi, tau) - tau * d2gamma1_dpitau(pi, tau)) /
                    (tau * tau * d2gamma1_dtau2(pi, tau)) -
                d2gamma1_dpi2(pi, tau));
      break;

    case 2:
      pi = pressure / _p_star[1];
      tau = _T_star[1] / temperature;
      speed2 = _Rw * temperature * Utility::pow<2>(pi * dgamma2_dpi(pi, tau)) /
               ((-pi * pi * d2gamma2_dpi2(pi, tau)) +
                Utility::pow<2>(pi * dgamma2_dpi(pi, tau) - tau * pi * d2gamma2_dpitau(pi, tau)) /
                    (tau * tau * d2gamma2_dtau2(pi, tau)));
      break;

    case 3:
    {
      // Calculate density first, then use that in Helmholtz free energy
      T density3 = densityRegion3(pressure, temperature);
      delta = density3 / _rho_critical;
      tau = _T_star[2] / temperature;
      speed2 =
          _Rw * temperature *
          (2.0 * delta * dphi3_ddelta(delta, tau) + delta * delta * d2phi3_ddelta2(delta, tau) -
           Utility::pow<2>(delta * dphi3_ddelta(delta, tau) -
                           delta * tau * d2phi3_ddeltatau(delta, tau)) /
               (tau * tau * d2phi3_dtau2(delta, tau)));
      break;
    }

    case 5:
      pi = pressure / _p_star[4];
      tau = _T_star[4] / temperature;
      speed2 = _Rw * temperature * Utility::pow<2>(pi * dgamma5_dpi(pi, tau)) /
               ((-pi * pi * d2gamma5_dpi2(pi, tau)) +
                Utility::pow<2>(pi * dgamma5_dpi(pi, tau) - tau * pi * d2gamma5_dpitau(pi, tau)) /
                    (tau * tau * d2gamma5_dtau2(pi, tau)));
      break;

    default:
      mooseError("inRegion() has given an incorrect region");
  }

  return std::sqrt(speed2);
}

template <typename T>
T
Water97FluidProperties::cp_from_p_T_template(const T & pressure, const T & temperature) const
{
  T specific_heat, pi, tau, delta;

  // Determine which region the point is in
  unsigned int region =
      inRegion(MetaPhysicL::raw_value(pressure), MetaPhysicL::raw_value(temperature));
  switch (region)
  {
    case 1:
      pi = pressure / _p_star[0];
      tau = _T_star[0] / temperature;
      specific_heat = -_Rw * tau * tau * d2gamma1_dtau2(pi, tau);
      break;

    case 2:
      pi = pressure / _p_star[1];
      tau = _T_star[1] / temperature;
      specific_heat = -_Rw * tau * tau * d2gamma2_dtau2(pi, tau);
      break;

    case 3:
    {
      // Calculate density first, then use that in Helmholtz free energy
      T density3 = densityRegion3(pressure, temperature);
      delta = density3 / _rho_critical;
      tau = _T_star[2] / temperature;
      specific_heat =
          _Rw *
          (-tau * tau * d2phi3_dtau2(delta, tau) +
           (delta * dphi3_ddelta(delta, tau) - delta * tau * d2phi3_ddeltatau(delta, tau)) *
               (delta * dphi3_ddelta(delta, tau) - delta * tau * d2phi3_ddeltatau(delta, tau)) /
               (2.0 * delta * dphi3_ddelta(delta, tau) +
                delta * delta * d2phi3_ddelta2(delta, tau)));
      break;
    }

    case 5:
      pi = pressure / _p_star[4];
      tau = _T_star[4] / temperature;
      specific_heat = -_Rw * tau * tau * d2gamma5_dtau2(pi, tau);
      break;

    default:
      mooseError("inRegion() has given an incorrect region");
  }
  return specific_heat;
}

template <typename T>
T
Water97FluidProperties::cv_from_p_T_template(const T & pressure, const T & temperature) const
{
  T specific_heat, pi, tau, delta;

  // Determine which region the point is in
  unsigned int region =
      inRegion(MetaPhysicL::raw_value(pressure), MetaPhysicL::raw_value(temperature));
  switch (region)
  {
    case 1:
      pi = pressure / _p_star[0];
      tau = _T_star[0] / temperature;
      specific_heat =
          _Rw * (-tau * tau * d2gamma1_dtau2(pi, tau) +
                 Utility::pow<2>(dgamma1_dpi(pi, tau) - tau * d2gamma1_dpitau(pi, tau)) /
                     d2gamma1_dpi2(pi, tau));
      break;

    case 2:
      pi = pressure / _p_star[1];
      tau = _T_star[1] / temperature;
      specific_heat =
          _Rw * (-tau * tau * d2gamma2_dtau2(pi, tau) +
                 Utility::pow<2>(dgamma2_dpi(pi, tau) - tau * d2gamma2_dpitau(pi, tau)) /
                     d2gamma2_dpi2(pi, tau));
      break;

    case 3:
    {
      // Calculate density first, then use that in Helmholtz free energy
      T density3 = densityRegion3(pressure, temperature);
      delta = density3 / _rho_critical;
      tau = _T_star[2] / temperature;
      specific_heat = _Rw * (-tau * tau * d2phi3_dtau2(delta, tau));
      break;
    }

    case 5:
      pi = pressure / _p_star[4];
      tau = _T_star[4] / temperature;
      specific_heat =
          _Rw * (-tau * tau * d2gamma5_dtau2(pi, tau) +
                 Utility::pow<2>(dgamma5_dpi(pi, tau) - tau * d2gamma5_dpitau(pi, tau)) /
                     d2gamma5_dpi2(pi, tau));
      break;

    default:
      mooseError("inRegion() has given an incorrect region");
  }
  return specific_heat;
}

template <typename T>
T
Water97FluidProperties::k_from_rho_T_template(const T & density, const T & temperature) const
{
  // Scale the density and temperature. Note that the scales are slightly
  // different to the critical values used in IAPWS-IF97
  T Tbar = temperature / 647.26;
  T rhobar = density / 317.7;

  // Ideal gas component
  T sum0 = 0.0;

  for (std::size_t i = 0; i < _k_a.size(); ++i)
    sum0 += _k_a[i] * MathUtils::pow(Tbar, i);

  T lambda0 = std::sqrt(Tbar) * sum0;

  // The contribution due to finite density
  T lambda1 = -0.39707 + 0.400302 * rhobar +
              1.06 * std::exp(-0.171587 * Utility::pow<2>(rhobar + 2.392190));

  // Critical enhancement
  T DeltaT = std::abs(Tbar - 1.0) + 0.00308976;
  T Q = 2.0 + 0.0822994 / std::pow(DeltaT, 0.6);
  T S = (Tbar >= 1.0 ? 1.0 / DeltaT : 10.0932 / std::pow(DeltaT, 0.6));

  T lambda2 = (0.0701309 / Utility::pow<10>(Tbar) + 0.011852) * std::pow(rhobar, 1.8) *
                  std::exp(0.642857 * (1.0 - std::pow(rhobar, 2.8))) +
              0.00169937 * S * std::pow(rhobar, Q) *
                  std::exp((Q / (1.0 + Q)) * (1.0 - std::pow(rhobar, 1.0 + Q))) -
              1.02 * std::exp(-4.11717 * std::pow(Tbar, 1.5) - 6.17937 / Utility::pow<5>(rhobar));

  return lambda0 + lambda1 + lambda2;
}

template <typename T>
T
Water97FluidProperties::k_from_p_T_template(const T & pressure, const T & temperature) const
{
  T rho = this->rho_from_p_T(pressure, temperature);
  return this->k_from_rho_T_template(rho, temperature);
}

template <typename T>
T
Water97FluidProperties::k_from_v_e_template(const T & v, const T & e) const
{
  const auto [p, temperature] = p_T_from_v_e(v, e);
  return k_from_p_T_template(p, temperature);
}

template <typename T>
T
Water97FluidProperties::p_from_v_e_template(const T & v, const T & e) const
{
  const auto rho = 1 / v;
  // For NewtonSolve
  // y = e
  // x = rho
  // z = p
  auto lambda = [this](const T & rho, const T & p, T & e, T & de_drho, T & de_dp)
  { e_from_p_rho(p, rho, e, de_dp, de_drho); };
  return FluidPropertiesUtils::NewtonSolve(rho, e, _p_initial_guess, _tolerance, lambda).first;
}

template <typename T>
T
Water97FluidProperties::e_from_p_rho_template(const T & p, const T & rho) const
{
  const auto temperature = T_drhodT_from_p_rho(p, rho).first;
  return e_from_p_T_template(p, temperature);
}

template <typename T>
void
Water97FluidProperties::e_from_p_rho_template(
    const T & p, const T & rho, T & e, T & de_dp, T & de_drho) const
{
  auto functor = [this](const auto & p, const auto & rho)
  { return this->e_from_p_rho_template(p, rho); };

  xyDerivatives(p, rho, e, de_dp, de_drho, functor);
}

template <typename T>
T
Water97FluidProperties::mu_from_rho_T_template(const T & density, const T & temperature) const
{
  const T mu_star = 1.e-6;
  const T rhobar = density / _rho_critical;
  const T Tbar = temperature / _T_critical;

  // Viscosity in limit of zero density
  T sum0 = 0.0;
  for (std::size_t i = 0; i < _mu_H0.size(); ++i)
    sum0 += _mu_H0[i] / MathUtils::pow(Tbar, i);

  const T mu0 = 100.0 * std::sqrt(Tbar) / sum0;

  // Residual component due to finite density
  T sum1 = 0.0;
  for (unsigned int i = 0; i < 6; ++i)
  {
    const T fact = MathUtils::pow(1.0 / Tbar - 1.0, i);
    for (unsigned int j = 0; j < 7; ++j)
      sum1 += fact * _mu_Hij[i][j] * MathUtils::pow(rhobar - 1.0, j);
  }

  const T mu1 = std::exp(rhobar * sum1);

  // The water viscosity (in Pa.s) is then given by
  return mu_star * mu0 * mu1;
}

template <typename T>
T
Water97FluidProperties::h_from_p_T_template(const T & pressure, const T & temperature) const
{
  T enthalpy, pi, tau, delta;

  // Determine which region the point is in
  unsigned int region =
      inRegion(MetaPhysicL::raw_value(pressure), MetaPhysicL::raw_value(temperature));
  switch (region)
  {
    case 1:
      pi = pressure / _p_star[0];
      tau = _T_star[0] / temperature;
      enthalpy = _Rw * _T_star[0] * dgamma1_dtau(pi, tau);
      break;

    case 2:
      pi = pressure / _p_star[1];
      tau = _T_star[1] / temperature;
      enthalpy = _Rw * _T_star[1] * dgamma2_dtau(pi, tau);
      break;

    case 3:
    {
      // Calculate density first, then use that in Helmholtz free energy
      T density3 = densityRegion3(pressure, temperature);
      delta = density3 / _rho_critical;
      tau = _T_star[2] / temperature;
      enthalpy =
          _Rw * temperature * (tau * dphi3_dtau(delta, tau) + delta * dphi3_ddelta(delta, tau));
      break;
    }

    case 5:
      pi = pressure / _p_star[4];
      tau = _T_star[4] / temperature;
      enthalpy = _Rw * _T_star[4] * dgamma5_dtau(pi, tau);
      break;

    default:
      mooseError("Water97FluidProperties::inRegion has given an incorrect region");
  }
  return enthalpy;
}

template <typename T>
void
Water97FluidProperties::h_from_p_T_template(
    const T & pressure, const T & temperature, T & h, T & dh_dp, T & dh_dT) const
{
  auto functor = [this](const auto & pressure, const auto & temperature)
  { return this->h_from_p_T_template(pressure, temperature); };

  xyDerivatives(pressure, temperature, h, dh_dp, dh_dT, functor);
}

template <typename T>
T
Water97FluidProperties::mu_from_p_T_template(const T & pressure, const T & temperature) const
{
  T rho = this->rho_from_p_T_template(pressure, temperature);
  return this->mu_from_rho_T_template(rho, temperature);
}

template <typename T>
T
Water97FluidProperties::s_from_p_T_template(const T & pressure, const T & temperature) const
{
  T entropy, pi, tau, density3, delta;

  // Determine which region the point is in
  unsigned int region =
      inRegion(MetaPhysicL::raw_value(pressure), MetaPhysicL::raw_value(temperature));
  switch (region)
  {
    case 1:
      pi = pressure / _p_star[0];
      tau = _T_star[0] / temperature;
      entropy = _Rw * (tau * dgamma1_dtau(pi, tau) - gamma1(pi, tau));
      break;

    case 2:
      pi = pressure / _p_star[1];
      tau = _T_star[1] / temperature;
      entropy = _Rw * (tau * dgamma2_dtau(pi, tau) - gamma2(pi, tau));
      break;

    case 3:
      // Calculate density first, then use that in Helmholtz free energy
      density3 = densityRegion3(pressure, temperature);
      delta = density3 / _rho_critical;
      tau = _T_star[2] / temperature;
      entropy = _Rw * (tau * dphi3_dtau(delta, tau) - phi3(delta, tau));
      break;

    case 5:
      pi = pressure / _p_star[4];
      tau = _T_star[4] / temperature;
      entropy = _Rw * (tau * dgamma5_dtau(pi, tau) - gamma5(pi, tau));
      break;

    default:
      mooseError("inRegion() has given an incorrect region");
  }
  return entropy;
}

template <typename T>
void
Water97FluidProperties::s_from_p_T_template(
    const T & pressure, const T & temperature, T & s, T & ds_dp, T & ds_dT) const
{
  auto functor = [this](const auto & pressure, const auto & temperature)
  { return this->s_from_p_T_template(pressure, temperature); };

  xyDerivatives(pressure, temperature, s, ds_dp, ds_dT, functor);
}

Real
Water97FluidProperties::rho_from_p_T(Real pressure, Real temperature) const
{
  return rho_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::rho_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return rho_from_p_T_template(pressure, temperature);
}

void
Water97FluidProperties::rho_from_p_T(
    Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho_from_p_T_template(pressure, temperature, rho, drho_dp, drho_dT);
}

void
Water97FluidProperties::rho_from_p_T(const ADReal & pressure,
                                     const ADReal & temperature,
                                     ADReal & rho,
                                     ADReal & drho_dp,
                                     ADReal & drho_dT) const
{
  rho_from_p_T_template(pressure, temperature, rho, drho_dp, drho_dT);
}

Real
Water97FluidProperties::v_from_p_T(Real pressure, Real temperature) const
{
  return v_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::v_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return v_from_p_T_template(pressure, temperature);
}

void
Water97FluidProperties::v_from_p_T(
    Real pressure, Real temperature, Real & v, Real & dv_dp, Real & dv_dT) const
{
  v_from_p_T_template(pressure, temperature, v, dv_dp, dv_dT);
}

void
Water97FluidProperties::v_from_p_T(const ADReal & pressure,
                                   const ADReal & temperature,
                                   ADReal & v,
                                   ADReal & dv_dp,
                                   ADReal & dv_dT) const
{
  v_from_p_T_template(pressure, temperature, v, dv_dp, dv_dT);
}

Real
Water97FluidProperties::p_from_v_e(const Real v, const Real e) const
{
  return p_from_v_e_template(v, e);
}

ADReal
Water97FluidProperties::p_from_v_e(const ADReal & v, const ADReal & e) const
{
  return p_from_v_e_template(v, e);
}

Real
Water97FluidProperties::e_from_p_rho(Real p, Real rho) const
{
  return e_from_p_rho_template(p, rho);
}

ADReal
Water97FluidProperties::e_from_p_rho(const ADReal & p, const ADReal & rho) const
{
  return e_from_p_rho_template(p, rho);
}

void
Water97FluidProperties::e_from_p_rho(
    const Real p, const Real rho, Real & e, Real & de_dp, Real & de_drho) const
{
  e_from_p_rho_template(p, rho, e, de_dp, de_drho);
}

void
Water97FluidProperties::e_from_p_rho(
    const ADReal & p, const ADReal & rho, ADReal & e, ADReal & de_dp, ADReal & de_drho) const
{
  e_from_p_rho_template(p, rho, e, de_dp, de_drho);
}

Real
Water97FluidProperties::e_from_p_T(Real pressure, Real temperature) const
{
  return e_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::e_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return e_from_p_T_template(pressure, temperature);
}

void
Water97FluidProperties::e_from_p_T(
    Real pressure, Real temperature, Real & e, Real & de_dp, Real & de_dT) const
{
  e_from_p_T_template(pressure, temperature, e, de_dp, de_dT);
}

void
Water97FluidProperties::e_from_p_T(const ADReal & pressure,
                                   const ADReal & temperature,
                                   ADReal & e,
                                   ADReal & de_dp,
                                   ADReal & de_dT) const
{
  e_from_p_T_template(pressure, temperature, e, de_dp, de_dT);
}

ADReal
Water97FluidProperties::e_from_v_h(const ADReal & v, const ADReal & h) const
{
  const auto [p, T] = p_T_from_v_h(v, h);
  return e_from_p_T(p, T);
}

Real
Water97FluidProperties::T_from_v_e(const Real v, const Real e) const
{
  return p_T_from_v_e(v, e).second;
}

ADReal
Water97FluidProperties::T_from_v_e(const ADReal & v, const ADReal & e) const
{
  return p_T_from_v_e(v, e).second;
}

ADReal
Water97FluidProperties::c_from_v_e(const ADReal & v, const ADReal & e) const
{
  const auto [p, T] = p_T_from_v_e(v, e);
  return c_from_p_T(p, T);
}

Real
Water97FluidProperties::c_from_p_T(const Real p, const Real T) const
{
  return c_from_p_T_template(p, T);
}

ADReal
Water97FluidProperties::c_from_p_T(const ADReal & p, const ADReal & T) const
{
  return c_from_p_T_template(p, T);
}

Real
Water97FluidProperties::cp_from_p_T(Real pressure, Real temperature) const
{
  return cp_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::cp_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return cp_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::cp_from_v_e(const ADReal & v, const ADReal & e) const
{
  const auto [p, T] = p_T_from_v_e(v, e);
  return cp_from_p_T(p, T);
}

Real
Water97FluidProperties::cv_from_p_T(Real pressure, Real temperature) const
{
  return cv_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::cv_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return cv_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::cv_from_v_e(const ADReal & v, const ADReal & e) const
{
  const auto [p, T] = p_T_from_v_e(v, e);
  return cv_from_p_T(p, T);
}

Real
Water97FluidProperties::mu_from_p_T(Real pressure, Real temperature) const
{
  return mu_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::mu_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return mu_from_p_T_template(pressure, temperature);
}

void
Water97FluidProperties::mu_from_p_T(
    Real pressure, Real temperature, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  Real rho, drho_dp, drho_dT;
  this->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
  Real dmu_drho;
  this->mu_from_rho_T(rho, temperature, drho_dT, mu, dmu_drho, dmu_dT);
  dmu_dp = dmu_drho * drho_dp;
}

Real
Water97FluidProperties::mu_from_rho_T(Real density, Real temperature) const
{
  return mu_from_rho_T_template(density, temperature);
}

void
Water97FluidProperties::mu_from_rho_T(Real density,
                                      Real temperature,
                                      Real ddensity_dT,
                                      Real & mu,
                                      Real & dmu_drho,
                                      Real & dmu_dT) const
{
  const Real mu_star = 1.0e-6;
  const Real rhobar = density / _rho_critical;
  const Real Tbar = temperature / _T_critical;
  const Real drhobar_drho = 1.0 / _rho_critical;
  const Real dTbar_dT = 1.0 / _T_critical;

  // Limit of zero density. Derivative wrt rho is 0
  Real sum0 = 0.0, dsum0_dTbar = 0.0;
  for (std::size_t i = 0; i < _mu_H0.size(); ++i)
  {
    sum0 += _mu_H0[i] / MathUtils::pow(Tbar, i);
    dsum0_dTbar -= i * _mu_H0[i] / MathUtils::pow(Tbar, i + 1);
  }

  const Real mu0 = 100.0 * std::sqrt(Tbar) / sum0;
  const Real dmu0_dTbar =
      50.0 / std::sqrt(Tbar) / sum0 - 100.0 * std::sqrt(Tbar) * dsum0_dTbar / sum0 / sum0;

  // Residual component due to finite density
  Real sum1 = 0.0, dsum1_drhobar = 0.0, dsum1_dTbar = 0.0;
  for (unsigned int i = 0; i < 6; ++i)
  {
    const Real fact = MathUtils::pow(1.0 / Tbar - 1.0, i);
    const Real dfact_dTbar = i * MathUtils::pow(1.0 / Tbar - 1.0, i - 1) / Tbar / Tbar;

    for (unsigned int j = 0; j < 7; ++j)
    {
      sum1 += fact * _mu_Hij[i][j] * MathUtils::pow(rhobar - 1.0, j);
      dsum1_dTbar -= dfact_dTbar * _mu_Hij[i][j] * MathUtils::pow(rhobar - 1.0, j);
      dsum1_drhobar += j * fact * _mu_Hij[i][j] * MathUtils::pow(rhobar - 1.0, j - 1);
    }
  }

  const Real mu1 = std::exp(rhobar * sum1);
  const Real dmu1_drhobar = (sum1 + rhobar * dsum1_drhobar) * mu1;
  const Real dmu1_dTbar = (rhobar * dsum1_dTbar) * mu1;

  // Viscosity and its derivatives are then
  mu = mu_star * mu0 * mu1;
  dmu_drho = mu_star * mu0 * dmu1_drhobar * drhobar_drho;
  dmu_dT = mu_star * (dmu0_dTbar * mu1 + mu0 * dmu1_dTbar) * dTbar_dT + dmu_drho * ddensity_dT;
}

ADReal
Water97FluidProperties::mu_from_v_e(const ADReal & v, const ADReal & e) const
{
  const auto [rho, T] = rho_T_from_v_e(v, e);
  return mu_from_rho_T_template(rho, T);
}

void
Water97FluidProperties::rho_mu_from_p_T(Real pressure,
                                        Real temperature,
                                        Real & rho,
                                        Real & mu) const
{
  rho = this->rho_from_p_T(pressure, temperature);
  mu = this->mu_from_rho_T(rho, temperature);
}

void
Water97FluidProperties::rho_mu_from_p_T(Real pressure,
                                        Real temperature,
                                        Real & rho,
                                        Real & drho_dp,
                                        Real & drho_dT,
                                        Real & mu,
                                        Real & dmu_dp,
                                        Real & dmu_dT) const
{
  this->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
  Real dmu_drho;
  this->mu_from_rho_T(rho, temperature, drho_dT, mu, dmu_drho, dmu_dT);
  dmu_dp = dmu_drho * drho_dp;
}

Real
Water97FluidProperties::k_from_p_T(Real pressure, Real temperature) const
{
  return k_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::k_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return k_from_p_T_template(pressure, temperature);
}

void
Water97FluidProperties::k_from_p_T(Real, Real, Real &, Real &, Real &) const
{
  mooseError("k_from_p_T() is not implemented.");
}

Real
Water97FluidProperties::k_from_rho_T(Real density, Real temperature) const
{
  return k_from_rho_T_template(density, temperature);
}

Real
Water97FluidProperties::k_from_v_e(Real v, Real e) const
{
  return k_from_v_e_template(v, e);
}

ADReal
Water97FluidProperties::k_from_v_e(const ADReal & v, const ADReal & e) const
{
  return k_from_v_e_template(v, e);
}

Real
Water97FluidProperties::s_from_p_T(Real pressure, Real temperature) const
{
  return s_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::s_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return s_from_p_T_template(pressure, temperature);
}

void
Water97FluidProperties::s_from_p_T(
    const Real pressure, const Real temperature, Real & s, Real & ds_dp, Real & ds_dT) const
{
  s_from_p_T_template(pressure, temperature, s, ds_dp, ds_dT);
}

void
Water97FluidProperties::s_from_p_T(const ADReal & pressure,
                                   const ADReal & temperature,
                                   ADReal & s,
                                   ADReal & ds_dp,
                                   ADReal & ds_dT) const
{
  s_from_p_T_template(pressure, temperature, s, ds_dp, ds_dT);
}

Real
Water97FluidProperties::h_from_p_T(Real pressure, Real temperature) const
{
  return h_from_p_T_template(pressure, temperature);
}

ADReal
Water97FluidProperties::h_from_p_T(const ADReal & pressure, const ADReal & temperature) const
{
  return h_from_p_T_template(pressure, temperature);
}

void
Water97FluidProperties::h_from_p_T(
    Real pressure, Real temperature, Real & h, Real & dh_dp, Real & dh_dT) const
{
  h_from_p_T_template(pressure, temperature, h, dh_dp, dh_dT);
}

void
Water97FluidProperties::h_from_p_T(const ADReal & pressure,
                                   const ADReal & temperature,
                                   ADReal & h,
                                   ADReal & dh_dp,
                                   ADReal & dh_dT) const
{
  h_from_p_T_template(pressure, temperature, h, dh_dp, dh_dT);
}

Real
Water97FluidProperties::vaporPressure(Real temperature) const
{
  // Check whether the input temperature is within the region of validity of this equation.
  // Valid for 273.15 K <= t <= 647.096 K
  if (temperature < 273.15 || temperature > _T_critical)
    mooseException(name(),
                   ": vaporPressure(): Temperature ",
                   temperature,
                   " is outside range 273.15 K <= T "
                   "<= 647.096 K");

  Real theta, theta2, a, b, c, p;
  theta = temperature + _n4[8] / (temperature - _n4[9]);
  theta2 = theta * theta;
  a = theta2 + _n4[0] * theta + _n4[1];
  b = _n4[2] * theta2 + _n4[3] * theta + _n4[4];
  c = _n4[5] * theta2 + _n4[6] * theta + _n4[7];
  p = Utility::pow<4>(2.0 * c / (-b + std::sqrt(b * b - 4.0 * a * c)));

  return p * 1.e6;
}

FPDualReal
Water97FluidProperties::vaporTemperature_ad(const FPDualReal & pressure) const
{
  // Check whether the input pressure is within the region of validity of this equation.
  // Valid for 611.213 Pa <= p <= 22.064 MPa
  if (pressure.value() < 611.23 || pressure.value() > _p_critical)
    mooseException(name() + ": vaporTemperature(): Pressure ",
                   pressure.value(),
                   " is outside range 611.213 Pa <= p <= 22.064 MPa");

  const FPDualReal beta = std::pow(pressure / 1.e6, 0.25);
  const FPDualReal beta2 = beta * beta;
  const FPDualReal e = beta2 + _n4[2] * beta + _n4[5];
  const FPDualReal f = _n4[0] * beta2 + _n4[3] * beta + _n4[6];
  const FPDualReal g = _n4[1] * beta2 + _n4[4] * beta + _n4[7];
  const FPDualReal d = 2.0 * g / (-f - std::sqrt(f * f - 4.0 * e * g));

  return (_n4[9] + d - std::sqrt((_n4[9] + d) * (_n4[9] + d) - 4.0 * (_n4[8] + _n4[9] * d))) / 2.0;
}

Real
Water97FluidProperties::vaporTemperature(Real pressure) const
{
  const FPDualReal p = pressure;

  return vaporTemperature_ad(p).value();
}

void
Water97FluidProperties::vaporTemperature(Real pressure, Real & Tsat, Real & dTsat_dp) const
{
  FPDualReal p = pressure;
  Moose::derivInsert(p.derivatives(), 0, 1.0);

  const FPDualReal T = vaporTemperature_ad(p);

  Tsat = T.value();
  dTsat_dp = T.derivatives()[0];
}

Real
Water97FluidProperties::b23p(Real temperature) const
{
  // Check whether the input temperature is within the region of validity of this equation.
  // Valid for 623.15 K <= t <= 863.15 K
  if (temperature < 623.15 || temperature > 863.15)
    mooseException(name(),
                   ": b23p(): Temperature ",
                   temperature,
                   " is outside range of 623.15 K <= T <= 863.15 K");

  return (_n23[0] + _n23[1] * temperature + _n23[2] * temperature * temperature) * 1.e6;
}

Real
Water97FluidProperties::b23T(Real pressure) const
{
  // Check whether the input pressure is within the region of validity of this equation.
  // Valid for 16.529 MPa <= p <= 100 MPa
  if (pressure < 16.529e6 || pressure > 100.0e6)
    mooseException(
        name(), ": b23T(): Pressure ", pressure, "is outside range 16.529 MPa <= p <= 100 MPa");

  return _n23[3] + std::sqrt((pressure / 1.e6 - _n23[4]) / _n23[2]);
}

unsigned int
Water97FluidProperties::inRegionPH(Real pressure, Real enthalpy) const
{
  unsigned int region;

  // Need to calculate enthalpies at the boundaries to delineate regions
  Real p273 = vaporPressure(273.15);
  Real p623 = vaporPressure(623.15);

  if (pressure >= p273 && pressure <= p623)
  {
    if (enthalpy >= h_from_p_T(pressure, 273.15) &&
        enthalpy <= h_from_p_T(pressure, vaporTemperature(pressure)))
      region = 1;
    else if (enthalpy > h_from_p_T(pressure, vaporTemperature(pressure)) &&
             enthalpy <= h_from_p_T(pressure, 1073.15))
      region = 2;
    else if (enthalpy > h_from_p_T(pressure, 1073.15) && enthalpy <= h_from_p_T(pressure, 2273.15))
      region = 5;
    else
      mooseException("Enthalpy ", enthalpy, " is out of range in ", name(), ": inRegionPH()");
  }
  else if (pressure > p623 && pressure <= 50.0e6)
  {
    if (enthalpy >= h_from_p_T(pressure, 273.15) && enthalpy <= h_from_p_T(pressure, 623.15))
      region = 1;
    else if (enthalpy > h_from_p_T(pressure, 623.15) &&
             enthalpy <= h_from_p_T(pressure, b23T(pressure)))
      region = 3;
    else if (enthalpy > h_from_p_T(pressure, b23T(pressure)) &&
             enthalpy <= h_from_p_T(pressure, 1073.15))
      region = 2;
    else if (enthalpy > h_from_p_T(pressure, 1073.15) && enthalpy <= h_from_p_T(pressure, 2273.15))
      region = 5;
    else
      mooseException("Enthalpy ", enthalpy, " is out of range in ", name(), ": inRegionPH()");
  }
  else if (pressure > 50.0e6 && pressure <= 100.0e6)
  {
    if (enthalpy >= h_from_p_T(pressure, 273.15) && enthalpy <= h_from_p_T(pressure, 623.15))
      region = 1;
    else if (enthalpy > h_from_p_T(pressure, 623.15) &&
             enthalpy <= h_from_p_T(pressure, b23T(pressure)))
      region = 3;
    else if (enthalpy > h_from_p_T(pressure, b23T(pressure)) &&
             enthalpy <= h_from_p_T(pressure, 1073.15))
      region = 2;
    else
      mooseException("Enthalpy ", enthalpy, " is out of range in ", name(), ": inRegionPH()");
  }
  else
    mooseException("Pressure ", pressure, " is out of range in ", name(), ": inRegionPH()");

  return region;
}

unsigned int
Water97FluidProperties::subregion2ph(Real pressure, Real enthalpy) const
{
  unsigned int subregion;

  if (pressure <= 4.0e6)
    subregion = 1;
  else if (pressure > 4.0e6 && pressure < 6.5467e6)
    subregion = 2;
  else
  {
    if (enthalpy >= b2bc(pressure))
      subregion = 2;
    else
      subregion = 3;
  }

  return subregion;
}

unsigned int
Water97FluidProperties::subregion3ph(Real pressure, Real enthalpy) const
{
  unsigned int subregion;

  if (enthalpy <= b3ab(pressure))
    subregion = 1;
  else
    subregion = 2;

  return subregion;
}

Real
Water97FluidProperties::b2bc(Real pressure) const
{
  // Check whether the input pressure is within the region of validity of this equation.
  if (pressure < 6.5467e6 || pressure > 100.0e6)
    mooseException(
        name(), ": b2bc(): Pressure ", pressure, " is outside range of 6.5467 MPa <= p <= 100 MPa");

  Real pi = pressure / 1.0e6;

  return (0.26526571908428e4 + std::sqrt((pi - 0.45257578905948e1) / 0.12809002730136e-3)) * 1.0e3;
}

Real
Water97FluidProperties::T_from_p_h(Real pressure, Real enthalpy) const
{
  const FPDualReal p = pressure;
  const FPDualReal h = enthalpy;

  return T_from_p_h_ad(p, h).value();
}

void
Water97FluidProperties::T_from_p_h(
    Real pressure, Real enthalpy, Real & temperature, Real & dT_dp, Real & dT_dh) const
{
  FPDualReal p = pressure;
  Moose::derivInsert(p.derivatives(), 0, 1.0);
  FPDualReal h = enthalpy;
  Moose::derivInsert(h.derivatives(), 1, 1.0);

  const FPDualReal T = T_from_p_h_ad(p, h);

  temperature = T.value();
  dT_dp = T.derivatives()[0];
  dT_dh = T.derivatives()[1];
}

FPDualReal
Water97FluidProperties::T_from_p_h_ad(const FPDualReal & pressure,
                                      const FPDualReal & enthalpy) const
{
  FPDualReal temperature = 0.0;

  // Determine which region the point is in
  const unsigned int region = inRegionPH(pressure.value(), enthalpy.value());

  switch (region)
  {
    case 1:
      temperature = temperature_from_ph1(pressure, enthalpy);
      break;

    case 2:
    {
      // First, determine which subregion the point is in:
      const unsigned int subregion = subregion2ph(pressure.value(), enthalpy.value());

      if (subregion == 1)
        temperature = temperature_from_ph2a(pressure, enthalpy);
      else if (subregion == 2)
        temperature = temperature_from_ph2b(pressure, enthalpy);
      else
        temperature = temperature_from_ph2c(pressure, enthalpy);
      break;
    }

    case 3:
    {
      // First, determine which subregion the point is in:
      const unsigned int subregion = subregion3ph(pressure.value(), enthalpy.value());

      if (subregion == 1)
        temperature = temperature_from_ph3a(pressure, enthalpy);
      else
        temperature = temperature_from_ph3b(pressure, enthalpy);
      break;
    }

    case 5:
      mooseError("temperature_from_ph() not implemented for region 5");
      break;

    default:
      mooseError("inRegionPH() has given an incorrect region");
  }

  return temperature;
}

FPDualReal
Water97FluidProperties::temperature_from_ph1(const FPDualReal & pressure,
                                             const FPDualReal & enthalpy) const
{
  const FPDualReal pi = pressure / 1.0e6;
  const FPDualReal eta = enthalpy / 2500.0e3;
  FPDualReal sum = 0.0;

  for (std::size_t i = 0; i < _nph1.size(); ++i)
    sum += _nph1[i] * std::pow(pi, _Iph1[i]) * std::pow(eta + 1.0, _Jph1[i]);

  return sum;
}

FPDualReal
Water97FluidProperties::temperature_from_ph2a(const FPDualReal & pressure,
                                              const FPDualReal & enthalpy) const
{
  const FPDualReal pi = pressure / 1.0e6;
  const FPDualReal eta = enthalpy / 2000.0e3;
  FPDualReal sum = 0.0;

  // Factor out the negative in std::pow(eta - 2.1, _Jph2a[i]) to avoid fpe in dbg (see #13163)
  const Real sgn = MathUtils::sign(eta.value() - 2.1);

  for (std::size_t i = 0; i < _nph2a.size(); ++i)
    sum += _nph2a[i] * std::pow(pi, _Iph2a[i]) * std::pow(std::abs(eta - 2.1), _Jph2a[i]) *
           std::pow(sgn, _Jph2a[i]);

  return sum;
}

FPDualReal
Water97FluidProperties::temperature_from_ph2b(const FPDualReal & pressure,
                                              const FPDualReal & enthalpy) const
{
  const FPDualReal pi = pressure / 1.0e6;
  const FPDualReal eta = enthalpy / 2000.0e3;
  FPDualReal sum = 0.0;

  // Factor out the negatives in std::pow(pi - 2.0, _Iph2b[i])* std::pow(eta - 2.6, _Jph2b[i])
  // to avoid fpe in dbg (see #13163)
  const Real sgn0 = MathUtils::sign(pi.value() - 2.0);
  const Real sgn1 = MathUtils::sign(eta.value() - 2.6);

  for (std::size_t i = 0; i < _nph2b.size(); ++i)
    sum += _nph2b[i] * std::pow(std::abs(pi - 2.0), _Iph2b[i]) * std::pow(sgn0, _Iph2b[i]) *
           std::pow(std::abs(eta - 2.6), _Jph2b[i]) * std::pow(sgn1, _Jph2b[i]);

  return sum;
}

FPDualReal
Water97FluidProperties::temperature_from_ph2c(const FPDualReal & pressure,
                                              const FPDualReal & enthalpy) const
{
  const FPDualReal pi = pressure / 1.0e6;
  const FPDualReal eta = enthalpy / 2000.0e3;
  FPDualReal sum = 0.0;

  // Factor out the negative in std::pow(eta - 1.8, _Jph2c[i]) to avoid fpe in dbg (see #13163)
  const Real sgn = MathUtils::sign(eta.value() - 1.8);

  for (std::size_t i = 0; i < _nph2c.size(); ++i)
    sum += _nph2c[i] * std::pow(pi + 25.0, _Iph2c[i]) * std::pow(std::abs(eta - 1.8), _Jph2c[i]) *
           std::pow(sgn, _Jph2c[i]);

  return sum;
}

Real
Water97FluidProperties::b3ab(Real pressure) const
{
  // Check whether the input pressure is within the region of validity of this equation.
  if (pressure < b23p(623.15) || pressure > 100.0e6)
    mooseException(
        name(), ": b3ab(): Pressure ", pressure, "is outside range of 16.529 MPa <= p <= 100 MPa");

  Real pi = pressure / 1.0e6;
  Real eta = 0.201464004206875e4 + 0.374696550136983e1 * pi - 0.219921901054187e-1 * pi * pi +
             0.87513168600995e-4 * pi * pi * pi;

  return eta * 1.0e3;
}

FPDualReal
Water97FluidProperties::temperature_from_ph3a(const FPDualReal & pressure,
                                              const FPDualReal & enthalpy) const
{
  const FPDualReal pi = pressure / 100.0e6;
  const FPDualReal eta = enthalpy / 2300.0e3;
  FPDualReal sum = 0.0;

  for (std::size_t i = 0; i < _nph3a.size(); ++i)
    sum += _nph3a[i] * std::pow(pi + 0.24, _Iph3a[i]) * std::pow(eta - 0.615, _Jph3a[i]);

  return sum * 760.0;
}

FPDualReal
Water97FluidProperties::temperature_from_ph3b(const FPDualReal & pressure,
                                              const FPDualReal & enthalpy) const
{
  const FPDualReal pi = pressure / 100.0e6;
  const FPDualReal eta = enthalpy / 2800.0e3;
  FPDualReal sum = 0.0;

  for (std::size_t i = 0; i < _nph3b.size(); ++i)
    sum += _nph3b[i] * std::pow(pi + 0.298, _Iph3b[i]) * std::pow(eta - 0.72, _Jph3b[i]);

  return sum * 860.0;
}

Real
Water97FluidProperties::henryConstant(Real temperature, const std::vector<Real> & coeffs) const
{
  const Real A = coeffs[0];
  const Real B = coeffs[1];
  const Real C = coeffs[2];

  const Real Tr = temperature / 647.096;
  const Real tau = 1.0 - Tr;

  const Real lnkh =
      A / Tr + B * std::pow(tau, 0.355) / Tr + C * std::pow(Tr, -0.41) * std::exp(tau);

  // The vapor pressure used in this formulation
  const std::vector<Real> a{
      -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
  const std::vector<Real> b{1.0, 1.5, 3.0, 3.5, 4.0, 7.5};
  Real sum = 0.0;

  for (std::size_t i = 0; i < a.size(); ++i)
    sum += a[i] * std::pow(tau, b[i]);

  return 22.064e6 * std::exp(sum / Tr) * std::exp(lnkh);
}

void
Water97FluidProperties::henryConstant(Real temperature,
                                      const std::vector<Real> & coeffs,
                                      Real & Kh,
                                      Real & dKh_dT) const
{
  const Real A = coeffs[0];
  const Real B = coeffs[1];
  const Real C = coeffs[2];

  const Real pc = 22.064e6;
  const Real Tc = 647.096;

  const Real Tr = temperature / Tc;
  const Real tau = 1.0 - Tr;

  const Real lnkh =
      A / Tr + B * std::pow(tau, 0.355) / Tr + C * std::pow(Tr, -0.41) * std::exp(tau);
  const Real dlnkh_dT =
      (-A / Tr / Tr - B * std::pow(tau, 0.355) / Tr / Tr - 0.355 * B * std::pow(tau, -0.645) / Tr -
       0.41 * C * std::pow(Tr, -1.41) * std::exp(tau) - C * std::pow(Tr, -0.41) * std::exp(tau)) /
      Tc;

  // The vapor pressure used in this formulation
  const std::vector<Real> a{
      -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
  const std::vector<Real> b{1.0, 1.5, 3.0, 3.5, 4.0, 7.5};
  Real sum = 0.0;
  Real dsum = 0.0;

  for (std::size_t i = 0; i < a.size(); ++i)
  {
    sum += a[i] * std::pow(tau, b[i]);
    dsum += a[i] * b[i] * std::pow(tau, b[i] - 1.0);
  }

  const Real p = pc * std::exp(sum / Tr);
  const Real dp_dT = -p / Tc / Tr * (sum / Tr + dsum);

  // Henry's constant and its derivative wrt temperature
  Kh = p * std::exp(lnkh);
  dKh_dT = (p * dlnkh_dT + dp_dT) * std::exp(lnkh);
}

DualReal
Water97FluidProperties::henryConstant(const DualReal & temperature,
                                      const std::vector<Real> & coeffs) const
{
  const Real T = temperature.value();
  Real Kh_real = 0.0;
  Real dKh_dT_real = 0.0;
  henryConstant(T, coeffs, Kh_real, dKh_dT_real);

  DualReal Kh = Kh_real;
  Kh.derivatives() = temperature.derivatives() * dKh_dT_real;

  return Kh;
}

unsigned int
Water97FluidProperties::subregion3(const Real pressure, const Real temperature) const
{
  Real pMPa = pressure / 1.0e6;
  const Real P3cd = 19.00881189173929;
  unsigned int subregion = 0;

  if (pMPa > 40.0 && pMPa <= 100.0)
  {
    if (temperature <= tempXY(pressure, AB))
      subregion = 0;
    else // (temperature > tempXY(pressure, AB))
      subregion = 1;
  }
  else if (pMPa > 25.0 && pMPa <= 40.0)
  {
    if (temperature <= tempXY(pressure, CD))
      subregion = 2;
    else if (temperature > tempXY(pressure, CD) && temperature <= tempXY(pressure, AB))
      subregion = 3;
    else if (temperature > tempXY(pressure, AB) && temperature <= tempXY(pressure, EF))
      subregion = 4;
    else // (temperature > tempXY(pressure, EF))
      subregion = 5;
  }
  else if (pMPa > 23.5 && pMPa <= 25.0)
  {
    if (temperature <= tempXY(pressure, CD))
      subregion = 2;
    else if (temperature > tempXY(pressure, CD) && temperature <= tempXY(pressure, GH))
      subregion = 6;
    else if (temperature > tempXY(pressure, GH) && temperature <= tempXY(pressure, EF))
      subregion = 7;
    else if (temperature > tempXY(pressure, EF) && temperature <= tempXY(pressure, IJ))
      subregion = 8;
    else if (temperature > tempXY(pressure, IJ) && temperature <= tempXY(pressure, JK))
      subregion = 9;
    else // (temperature > tempXY(pressure, JK))
      subregion = 10;
  }
  else if (pMPa > 23.0 && pMPa <= 23.5)
  {
    if (temperature <= tempXY(pressure, CD))
      subregion = 2;
    else if (temperature > tempXY(pressure, CD) && temperature <= tempXY(pressure, GH))
      subregion = 11;
    else if (temperature > tempXY(pressure, GH) && temperature <= tempXY(pressure, EF))
      subregion = 7;
    else if (temperature > tempXY(pressure, EF) && temperature <= tempXY(pressure, IJ))
      subregion = 8;
    else if (temperature > tempXY(pressure, IJ) && temperature <= tempXY(pressure, JK))
      subregion = 9;
    else // (temperature > tempXY(pressure, JK))
      subregion = 10;
  }
  else if (pMPa > 22.5 && pMPa <= 23.0)
  {
    if (temperature <= tempXY(pressure, CD))
      subregion = 2;
    else if (temperature > tempXY(pressure, CD) && temperature <= tempXY(pressure, GH))
      subregion = 11;
    else if (temperature > tempXY(pressure, GH) && temperature <= tempXY(pressure, MN))
      subregion = 12;
    else if (temperature > tempXY(pressure, MN) && temperature <= tempXY(pressure, EF))
      subregion = 13;
    else if (temperature > tempXY(pressure, EF) && temperature <= tempXY(pressure, OP))
      subregion = 14;
    else if (temperature > tempXY(pressure, OP) && temperature <= tempXY(pressure, IJ))
      subregion = 15;
    else if (temperature > tempXY(pressure, IJ) && temperature <= tempXY(pressure, JK))
      subregion = 9;
    else // (temperature > tempXY(pressure, JK))
      subregion = 10;
  }
  else if (pMPa > vaporPressure(643.15) * 1.0e-6 &&
           pMPa <= 22.5) // vaporPressure(643.15) = 21.04 MPa
  {
    if (temperature <= tempXY(pressure, CD))
      subregion = 2;
    else if (temperature > tempXY(pressure, CD) && temperature <= tempXY(pressure, QU))
      subregion = 16;
    else if (temperature > tempXY(pressure, QU) && temperature <= tempXY(pressure, RX))
    {
      if (pMPa > 22.11 && pMPa <= 22.5)
      {
        if (temperature <= tempXY(pressure, UV))
          subregion = 20;
        else if (temperature > tempXY(pressure, UV) && temperature <= tempXY(pressure, EF))
          subregion = 21;
        else if (temperature > tempXY(pressure, EF) && temperature <= tempXY(pressure, WX))
          subregion = 22;
        else // (temperature > tempXY(pressure, WX) && temperature <= tempXY(pressure, RX))
          subregion = 23;
      }
      else if (pMPa > 22.064 && pMPa <= 22.11)
      {
        if (temperature <= tempXY(pressure, UV))
          subregion = 20;
        else if (temperature > tempXY(pressure, UV) && temperature <= tempXY(pressure, EF))
          subregion = 24;
        else if (temperature > tempXY(pressure, EF) && temperature <= tempXY(pressure, WX))
          subregion = 25;
        else // (temperature > tempXY(pressure, WX) && temperature <= tempXY(pressure, RX))
          subregion = 23;
      }
      else if (temperature <= vaporTemperature(pressure))
      {
        if (pMPa > 21.93161551 && pMPa <= 22.064)
          if (temperature > tempXY(pressure, QU) && temperature <= tempXY(pressure, UV))
            subregion = 20;
          else
            subregion = 24;
        else // (pMPa > vaporPressure(643.15) * 1.0e-6 && pMPa <= 21.93161551)
          subregion = 20;
      }
      else if (temperature > vaporTemperature(pressure))
      {
        if (pMPa > 21.90096265 && pMPa <= 22.064)
        {
          if (temperature <= tempXY(pressure, WX))
            subregion = 25;
          else
            subregion = 23;
        }
        else
          subregion = 23;
      }
    }
    else if (temperature > tempXY(pressure, RX) && temperature <= tempXY(pressure, JK))
      subregion = 17;
    else
      subregion = 10;
  }
  else if (pMPa > 20.5 &&
           pMPa <= vaporPressure(643.15) * 1.0e-6) // vaporPressure(643.15) = 21.04 MPa
  {
    if (temperature <= tempXY(pressure, CD))
      subregion = 2;
    else if (temperature > tempXY(pressure, CD) && temperature <= vaporTemperature(pressure))
      subregion = 16;
    else if (temperature > vaporTemperature(pressure) && temperature <= tempXY(pressure, JK))
      subregion = 17;
    else // (temperature > tempXY(pressure, JK))
      subregion = 10;
  }
  else if (pMPa > P3cd && pMPa <= 20.5) // P3cd  = 19.00881189173929
  {
    if (temperature <= tempXY(pressure, CD))
      subregion = 2;
    else if (temperature > tempXY(pressure, CD) && temperature <= vaporTemperature(pressure))
      subregion = 18;
    else
      subregion = 19;
  }
  else if (pMPa > vaporPressure(623.15) * 1.0e-6 && pMPa <= P3cd)
  {
    if (temperature < vaporTemperature(pressure))
      subregion = 2;
    else
      subregion = 19;
  }
  else if (pMPa > 22.11 && pMPa <= 22.5)
  {
    if (temperature > tempXY(pressure, QU) && temperature <= tempXY(pressure, UV))
      subregion = 20;
    else if (temperature > tempXY(pressure, UV) && temperature <= tempXY(pressure, EF))
      subregion = 21;
    else if (temperature > tempXY(pressure, EF) && temperature <= tempXY(pressure, WX))
      subregion = 22;
    else // (temperature > tempXY(pressure, WX) && temperature <= tempXY(pressure, RX))
      subregion = 23;
  }
  else if (pMPa > 22.064 && pMPa <= 22.11)
  {
    if (temperature > tempXY(pressure, QU) && temperature <= tempXY(pressure, UV))
      subregion = 20;
    else if (temperature > tempXY(pressure, UV) && temperature <= tempXY(pressure, EF))
      subregion = 24;
    else if (temperature > tempXY(pressure, EF) && temperature <= tempXY(pressure, WX))
      subregion = 25;
    else // (temperature > tempXY(pressure, WX) && temperature <= tempXY(pressure, RX))
      subregion = 23;
  }
  else
    mooseError("subregion3(): Shouldn't have got here!");

  return subregion;
}

unsigned int
Water97FluidProperties::inRegion(const Real pressure, const Real temperature) const
{
  // Valid for 273.15 K <= T <= 1073.15 K, p <= 100 MPa
  //          1073.15 K <= T <= 2273.15 K, p <= 50 Mpa
  if (temperature >= 273.15 && temperature <= 1073.15)
  {
    if (pressure < vaporPressure(273.15) || pressure > 100.0e6)
      mooseException("Pressure ", pressure, " is out of range in ", name(), ": inRegion()");
  }
  else if (temperature > 1073.15 && temperature <= 2273.15)
  {
    if (pressure < 0.0 || pressure > 50.0e6)
      mooseException("Pressure ", pressure, " is out of range in ", name(), ": inRegion()");
  }
  else
    mooseException("Temperature ", temperature, " is out of range in ", name(), ": inRegion()");

  // Determine the phase region that the (P, T) point lies in
  unsigned int region;

  if (temperature >= 273.15 && temperature <= 623.15)
  {
    if (pressure > vaporPressure(temperature) && pressure <= 100.0e6)
      region = 1;
    else
      region = 2;
  }
  else if (temperature > 623.15 && temperature <= 863.15)
  {
    if (pressure <= b23p(temperature))
      region = 2;
    else
      region = 3;
  }
  else if (temperature > 863.15 && temperature <= 1073.15)
    region = 2;
  else
    region = 5;

  return region;
}
