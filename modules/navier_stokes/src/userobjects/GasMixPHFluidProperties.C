/********************************************************************/
/*                   DO NOT MODIFY THIS HEADER                      */
/*   Pronghorn: Coarse-Mesh, Multi-Dimensional, Thermal-Hydraulics  */
/*                                                                  */
/*              (c) 2020 Battelle Energy Alliance, LLC              */
/*                      ALL RIGHTS RESERVED                         */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*               With the U. S. Department of Energy                */
/*                                                                  */
/*               See COPYRIGHT for full restrictions                */
/********************************************************************/

#include "GasMixPHFluidProperties.h"
#include <numeric>

registerMooseObject("NavierStokesApp", GasMixPHFluidProperties);

const Real GasMixPHFluidProperties::R_molar = 8.3144598;

InputParameters
GasMixPHFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();

  params.addClassDescription("Pronghorn class for fluid properties of an binary gas mixture");

  params.addRequiredParam<UserObjectName>(
      "fp_primary", "Name of fluid properties user object for primary gas component");
  params.addRequiredParam<std::vector<UserObjectName>>(
      "fp_secondary", "Name of fluid properties user object(s) for secondary gas component(s)");

  // This is necessary because initialize() must be called before any interface
  // can be used (which can occur as early as initialization of variables).
  params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;  // Not sure why

  return params;
}

GasMixPHFluidProperties::GasMixPHFluidProperties(
    const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    _fp_primary(&getUserObject<SinglePhaseFluidProperties>("fp_primary")),
    _fp_secondary_names(getParam<std::vector<UserObjectName>>("fp_secondary")),
    _n_secondary_gas(_fp_secondary_names.size())
{
  _fp_secondary.resize(_n_secondary_gas);
  for (unsigned int i = 0; i < _n_secondary_gas; i++)
    _fp_secondary[i] = &getUserObjectByName<SinglePhaseFluidProperties>(_fp_secondary_names[i]);
}


GasMixPHFluidProperties::~GasMixPHFluidProperties() {}

Real
GasMixPHFluidProperties::p_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const
{
  Real xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
  Real vi = v / xp, ei =  _fp_primary->e_from_T_v(T, vi);
  Real p = _fp_primary->p_from_v_e(vi, ei);
  for (unsigned int i = 0; i < _n_secondary_gas; i++)
  {
      vi = v / x[i];
      ei = _fp_secondary[i]->e_from_T_v(T, vi);
      p += _fp_secondary[i]->p_from_v_e(vi, ei);
  }
  return p;
}

ADReal
GasMixPHFluidProperties::p_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const
{
  ADReal xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
  ADReal vi = v / xp, ei =  _fp_primary->e_from_T_v(T, vi);
  ADReal p = _fp_primary->p_from_v_e(vi, ei);
  for (unsigned int i = 0; i < _n_secondary_gas; i++)
  {
      vi = v / x[i];
      ei = _fp_secondary[i]->e_from_T_v(T, vi);
      p += _fp_secondary[i]->p_from_v_e(vi, ei);
  }
  return p;
}

void
GasMixPHFluidProperties::p_from_v_e_X(Real v, Real e, const std::vector<Real> & x, Real & p, Real & dp_dv, Real & dp_de, Real & dp_dx) const
{
  p = p_from_v_e_X(v, e, x);
  Real xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
  Real vi = v / xp, ei = _fp_primary->e_from_T_v(T, vi);
  Real dp_dv_pr, dp_dv_sec, dp_de_pr, dp_de_sec, dummy;
  _fp_primary->p_from_v_e(vi, ei, dummy, dp_dv_pr, dp_de_pr);
  dp_dv = dp_dv_pr / xp;
  dp_de = dp_de_pr;
  for (unsigned int i = 0; i < _n_secondary_gas; i++)
  {
      vi = v / x[i];
      ei = _fp_secondary[i]->e_from_T_v(T, vi);
      _fp_secondary[i]->p_from_v_e(vi, ei, dummy, dp_dv_sec, dp_de_sec);
      dp_dv += dp_dv_sec / x[i];
      dp_de += dp_de_sec;
  }
  dp_dx = dp_dv_pr * v / xp / xp - dp_dv_sec * v / x[0] / x[0];
}

void
GasMixPHFluidProperties::p_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x, ADReal & p, ADReal & dp_dv, ADReal & dp_de, ADReal & dp_dx) const
{
  p = p_from_v_e_X(v, e, x);
  ADReal xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
  ADReal vi = v / xp, ei = _fp_primary->e_from_T_v(T, vi);
  ADReal dp_dv_pr, dp_dv_sec, dp_de_pr, dp_de_sec, dummy;
  _fp_primary->p_from_v_e(vi, ei, dummy, dp_dv_pr, dp_de_pr);
  dp_dv = dp_dv_pr / xp;
  dp_de = dp_de_pr;
  for (unsigned int i = 0; i < _n_secondary_gas; i++)
  {
      vi = v / x[i];
      ei = _fp_secondary[i]->e_from_T_v(T, vi);
      _fp_secondary[i]->p_from_v_e(vi, ei, dummy, dp_dv_sec, dp_de_sec);
      dp_dv += dp_dv_sec / x[i];
      dp_de += dp_de_sec;
  }
  dp_dx = dp_dv_pr * v / xp / xp - dp_dv_sec * v / x[0] / x[0];
}

Real
GasMixPHFluidProperties::T_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const
{
  return e / cv_from_v_e_X(v, e, x);
}

ADReal
GasMixPHFluidProperties::T_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const
{
  return e / cv_from_v_e_X(v, e, x);
}

void
GasMixPHFluidProperties::T_from_v_e_X(Real v, Real e, const std::vector<Real> & x, Real & T, Real & dT_dv, Real & dT_de, Real & dT_dx) const
{
  T = T_from_v_e_X(v, e, x);
  Real cv = cv_from_v_e_X(v, e, x);
  Real vi = v / x[0];
  Real dcv_dx = _fp_secondary[0]->cv_from_v_e(vi, e);
  Real de_dx = _fp_secondary[0]->e_from_T_v(T, vi);
  dT_dv = 0.0;
  dT_de = 1.0 / cv;
  dT_dx = (de_dx * cv - dcv_dx * e) / cv / cv;
}

void
GasMixPHFluidProperties::T_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x, ADReal & T, ADReal & dT_dv, ADReal & dT_de, ADReal & dT_dx) const
{
  T = T_from_v_e_X(v, e, x);
  ADReal cv = cv_from_v_e_X(v, e, x);
  ADReal vi = v / x[0];
  ADReal dcv_dx = _fp_secondary[0]->cv_from_v_e(vi, e);
  ADReal de_dx = _fp_secondary[0]->e_from_T_v(T, vi);
  dT_dv = 0.0;
  dT_de = 1.0 / cv;
  dT_dx = (de_dx * cv - dcv_dx * e) / cv / cv;
}

Real
GasMixPHFluidProperties::c_from_v_e_X(Real v, Real e,  const std::vector<Real> & x) const
{
  Real T = T_from_v_e_X(v, e, x);
  Real gamma = cp_from_v_e_X(v, e, x) / cv_from_v_e_X(v, e, x);
  Real R_specific = R_molar / molarMass_from_X(x);
  const Real c2 = gamma * R_specific * T;
  if (c2 < 0)
    mooseException(name() + ": Sound speed squared (gamma * R * T) is negative: c2 = " +
                   Moose::stringify(c2) + ".");
  return std::sqrt(c2);
}

ADReal
GasMixPHFluidProperties::c_from_v_e_X(ADReal v, ADReal e,  const std::vector<ADReal> & x) const
{
  ADReal T = T_from_v_e_X(v, e, x);
  ADReal gamma = cp_from_v_e_X(v, e, x) / cv_from_v_e_X(v, e, x);
  ADReal R_specific = R_molar / molarMass_from_X(x);
  const ADReal c2 = gamma * R_specific * T;
  if (c2 < 0)
    mooseException(name() + ": Sound speed squared (gamma * R * T) is negative: c2 = " +
                   Moose::stringify(c2) + ".");
  return std::sqrt(c2);
}

void
GasMixPHFluidProperties::c_from_v_e_X(Real v, Real e, const std::vector<Real> & x, Real c, Real dc_dv, Real dc_de, Real dc_dx) const
{
  Real T, dT_dv, dT_de, dT_dx;
  T_from_v_e_X(v, e, x, T, dT_dv, dT_de, dT_dx);
  Real gamma = cp_from_v_e_X(v, e, x) / cv_from_v_e_X(v, e, x);
  Real R_specific = R_molar / molarMass_from_X(x);

  c = std::sqrt(gamma * R_specific * T);
  const Real dc_dT = 0.5 / c * gamma * R_specific;

  dc_dv = dc_dT * dT_dv;
  dc_de = dc_dT * dT_de;
  dc_dx = dc_dT * dT_dx;
}

void
GasMixPHFluidProperties::c_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x, ADReal c, ADReal dc_dv, ADReal dc_de, ADReal dc_dx) const
{
  ADReal T, dT_dv, dT_de, dT_dx;
  T_from_v_e_X(v, e, x, T, dT_dv, dT_de, dT_dx);
  ADReal gamma = cp_from_v_e_X(v, e, x) / cv_from_v_e_X(v, e, x);
  ADReal R_specific = R_molar / molarMass_from_X(x);

  c = std::sqrt(gamma * R_specific * T);
  const ADReal dc_dT = 0.5 / c * gamma * R_specific;

  dc_dv = dc_dT * dT_dv;
  dc_de = dc_dT * dT_de;
  dc_dx = dc_dT * dT_dx;
}

Real
GasMixPHFluidProperties::c_from_p_T_x(Real p, Real T, const std::vector<Real> & x) const
{
  Real gamma = cp_from_p_T_X(p, T, x) / cv_from_p_T_X(p, T, x);
  Real R_specific = R_molar / molarMass_from_X(x);
  return std::sqrt(gamma * R_specific* T) ;
}

ADReal
GasMixPHFluidProperties::c_from_p_T_x(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
  ADReal gamma = cp_from_p_T_X(p, T, x) / cv_from_p_T_X(p, T, x);
  ADReal R_specific = R_molar / molarMass_from_X(x);
  return std::sqrt(gamma * R_specific* T) ;
}

Real GasMixPHFluidProperties::cp_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const
{
    Real xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
    Real vi = v / xp, ei =  _fp_primary->e_from_T_v(T, vi);
    Real cp = xp * _fp_primary->cp_from_v_e(vi, ei);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
    {
        vi = v / x[i];
        ei = _fp_secondary[i]->e_from_T_v(T, vi);
        cp += x[i] * _fp_secondary[i]->cp_from_v_e(vi, ei);
    }
    return cp;
}

ADReal GasMixPHFluidProperties::cp_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const
{
    ADReal xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
    ADReal vi = v / xp, ei =  _fp_primary->e_from_T_v(T, vi);
    ADReal cp = xp * _fp_primary->cp_from_v_e(vi, ei);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
    {
        vi = v / x[i];
        ei = _fp_secondary[i]->e_from_T_v(T, vi);
        cp += x[i] * _fp_secondary[i]->cp_from_v_e(vi, ei);
    }
    return cp;
}

Real GasMixPHFluidProperties::cv_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const
{
    Real xp = xp_from_X(x);
    //Real vi = v / xp, ei =  _fp_primary->e_from_T_v(T, vi);
    Real cv = xp * _fp_primary->cv_from_v_e(v, e);
    for (unsigned int i = 0; i < _n_secondary_gas; i++){
        //vi = v / x[i];
        //ei = _fp_secondary[i]->e_from_T_v(T, vi);
        cv += x[i] * _fp_secondary[i]->cv_from_v_e(v, e);
    }
    return cv;
}

ADReal GasMixPHFluidProperties::cv_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const
{
    ADReal xp = xp_from_X(x);
    //ADReal vi = v / xp, ei =  _fp_primary->e_from_p_rho(pi, 1 / vi);
    ADReal cv = xp * _fp_primary->cv_from_v_e(v, e);
    for (unsigned int i = 0; i < _n_secondary_gas; i++){
        //vi = v / x[i];
        //ei = _fp_secondary[i]->e_from_p_rho(pi, 1 / vi);
        cv += x[i] * _fp_secondary[i]->cv_from_v_e(v, e);
    }
    return cv;
}

Real
GasMixPHFluidProperties::mu_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const
{
  Real mu_g = 0.0;
  Real g_ij = 0.0, g_i = 0.0;
  Real mu_ij, M_ij;
  Real multiplier = std::pow(2, 0.5) / 4;
  Real M_g = molarMass_from_X(x);
  Real sum_y = 0.0, T = T_from_v_e_X(v, e, x);
  Real vi = v / xp_from_X(x), ei = _fp_primary->e_from_T_v(T, vi);

  std::vector<Real> mass_mol(_n_secondary_gas+1);
  std::vector<Real> yi(_n_secondary_gas+1);
  std::vector<Real> mui(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  mui[0] = _fp_primary->mu_from_v_e(vi, ei);
  for (unsigned int i = 1; i < _n_secondary_gas+1; i++)
  {
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    vi = v / x[i-1];
    ei = _fp_primary->e_from_T_v(T, vi);
    mui[i] = _fp_secondary[i-1]->mu_from_v_e(vi, ei);
    sum_y += yi[i];
  }
  yi[0] = 1.0 - sum_y;
  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            mu_g += mui[i] / (1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return mu_g;
}

ADReal
GasMixPHFluidProperties::mu_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const
{
  ADReal mu_g = 0.0;
  ADReal g_ij = 0.0, g_i = 0.0;
  ADReal mu_ij, M_ij;
  ADReal multiplier = std::pow(2, 0.5) / 4;
  ADReal M_g = molarMass_from_X(x);
  ADReal sum_y = 0.0, T = T_from_v_e_X(v, e, x);
  ADReal vi = v / xp_from_X(x), ei = _fp_primary->e_from_T_v(T, vi);

  std::vector<ADReal> mass_mol(_n_secondary_gas+1);
  std::vector<ADReal> yi(_n_secondary_gas+1);
  std::vector<ADReal> mui(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  mui[0] = _fp_primary->mu_from_v_e(vi, ei);
  for (unsigned int i = 1; i < _n_secondary_gas+1; i++)
  {
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    vi = v / x[i-1];
    ei = _fp_primary->e_from_T_v(T, vi);
    mui[i] = _fp_secondary[i-1]->mu_from_v_e(vi, ei);
    sum_y += yi[i];
  }
  yi[0] = 1.0 - sum_y;
  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            mu_g += mui[i] / (1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return mu_g;
}

Real
GasMixPHFluidProperties::k_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const
{
  Real k_g = 0.0;
  Real g_ij = 0.0, g_i = 0.0;
  Real mu_ij, M_ij;
  Real multiplier = std::pow(2, 0.5) / 4;
  Real M_g = molarMass_from_X(x);
  Real xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
  Real vi = v / xp_from_X(x), ei = _fp_primary->e_from_T_v(T, vi);

  std::vector<Real> mass_mol(_n_secondary_gas+1);
  std::vector<Real> yi(_n_secondary_gas+1);
  std::vector<Real> mui(_n_secondary_gas+1);
  std::vector<Real> ki(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  yi[0] = xp_from_X(x) * M_g / mass_mol[0];
  mui[0] = _fp_primary->mu_from_v_e(vi, ei);
  ki[0] = _fp_primary->k_from_v_e(vi, ei);

  for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    vi = v / x[i-1];
    ei = _fp_primary->e_from_T_v(T, vi);
    mui[i] = _fp_secondary[i-1]->mu_from_v_e(vi, ei);
    ki[i] = _fp_secondary[i-1]->k_from_v_e(vi, ei);
  }

  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            k_g += ki[i] /(1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return k_g;
}

ADReal
GasMixPHFluidProperties::k_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const
{
  ADReal k_g = 0.0;
  ADReal g_ij = 0.0, g_i = 0.0;
  ADReal mu_ij, M_ij;
  ADReal multiplier = std::pow(2, 0.5) / 4;
  ADReal M_g = molarMass_from_X(x);
  ADReal xp = xp_from_X(x), T = T_from_v_e_X(v, e, x);
  ADReal vi = v / xp_from_X(x), ei = _fp_primary->e_from_T_v(T, vi);

  std::vector<ADReal> mass_mol(_n_secondary_gas+1);
  std::vector<ADReal> yi(_n_secondary_gas+1);
  std::vector<ADReal> mui(_n_secondary_gas+1);
  std::vector<ADReal> ki(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  yi[0] = xp_from_X(x) * M_g / mass_mol[0];
  mui[0] = _fp_primary->mu_from_v_e(vi, ei);
  ki[0] = _fp_primary->k_from_v_e(vi, ei);

  for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    vi = v / x[i-1];
    ei = _fp_primary->e_from_T_v(T, vi);
    mui[i] = _fp_secondary[i-1]->mu_from_v_e(vi, ei);
    ki[i] = _fp_secondary[i-1]->k_from_v_e(vi, ei);
  }

  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            k_g += ki[i] /(1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return k_g;
}

Real
GasMixPHFluidProperties::beta_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
    Real rho, drho_dp, drho_dT;
    rho_from_p_T_X(p, T, x, rho, drho_dp, drho_dT);
    return -drho_dT / rho;
}

Real
GasMixPHFluidProperties:: rho_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
    /*Real xp = xp_from_X(x);
    Real M_g = molarMass_from_X(x);
    Real pi = p * (xp * M_g / _fp_primary->molarMass() );
    Real rho =  xp_from_X(x) / _fp_primary->rho_from_p_T(pi, T);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
    {
        pi = p * (x[i] * M_g / _fp_secondary[i]->molarMass());
        rho += x[i] / _fp_secondary[i]->rho_from_p_T(pi, T);
    }
    return 1 / rho;
    */
    Real R_specific = R_molar / molarMass_from_X(x);
    return p / R_specific / T;
}

ADReal
GasMixPHFluidProperties:: rho_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
    /*ADReal xp = xp_from_X(x);
    ADReal M_g = molarMass_from_X(x);
    ADReal pi = p * (xp * M_g / _fp_primary->molarMass() );
    ADReal rho =  xp_from_X(x) / _fp_primary->rho_from_p_T(pi, T);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
    {
        pi = p * (x[i] * M_g / _fp_secondary[i]->molarMass());
        rho += x[i] / _fp_secondary[i]->rho_from_p_T(pi, T);
    }
    return 1 / rho;
*/
    ADReal R_specific = R_molar / molarMass_from_X(x);
    return p / R_specific / T;
}

void
GasMixPHFluidProperties:: rho_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & rho, Real & drho_dp, Real & drho_dT) const
{
    /*rho =  rho_from_p_T_X(p, T, x);
    // Derivatives calculation
    Real xp = xp_from_X(x);
    Real M_g = molarMass_from_X(x);
    Real p_i = p * (xp * M_g / _fp_primary->molarMass());
    Real rho_i = 0., drhoi_dT = 0., drhoi_dp = 0.;
    _fp_primary->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
    Real dv_dT = xp_from_X(x) * drhoi_dT / std::pow(rho_i, 2);
    Real dv_dp = xp_from_X(x) * drhoi_dp / std::pow(rho_i, 2);
    for (unsigned int i = 0; i < _n_secondary_gas; i++){
        p_i = p * (x[i] * M_g / _fp_secondary[i]->molarMass());
        _fp_secondary[i]->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
        dv_dT += x[i] * drhoi_dT / std::pow(rho_i, 2);
        dv_dp += x[i] * drhoi_dp / std::pow(rho_i, 2);
    }
    drho_dp = std::pow(rho, 2) * dv_dp;
    drho_dT = std::pow(rho, 2) * dv_dT;
    */
    Real R_specific = R_molar / molarMass_from_X(x);
    rho = p / R_specific / T;
    drho_dp = 1 / R_specific / T;
    drho_dT = - p / R_specific / T / T;
}


void
GasMixPHFluidProperties:: rho_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x, ADReal & rho, ADReal & drho_dp, ADReal & drho_dT) const
{
    /*rho =  rho_from_p_T_X(p, T, x);
    // Derivatives calculation
    ADReal xp = xp_from_X(x);
    ADReal M_g = molarMass_from_X(x);
    ADReal p_i = p * (xp * M_g / _fp_primary->molarMass()  );
    ADReal rho_i = 0., drhoi_dT = 0., drhoi_dp = 0.;
    _fp_primary->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
    ADReal dv_dT = xp_from_X(x) * drhoi_dT / std::pow(rho_i, 2);
    ADReal dv_dp = xp_from_X(x) * drhoi_dp / std::pow(rho_i, 2);
    for (unsigned int i = 0; i < _n_secondary_gas; i++){
        p_i = p * (x[i] * M_g / _fp_secondary[i]->molarMass());
        _fp_secondary[i]->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
        dv_dT += x[i] * drhoi_dT / std::pow(rho_i, 2);
        dv_dp += x[i] * drhoi_dp / std::pow(rho_i, 2);
    }
    drho_dp = std::pow(rho, 2) * dv_dp;
    drho_dT = std::pow(rho, 2) * dv_dT;
    */
    ADReal R_specific = R_molar / molarMass_from_X(x);
    rho = p / R_specific / T;
    drho_dp = 1 / R_specific / T;
    drho_dT = - p / R_specific / T / T;
}

void
GasMixPHFluidProperties:: rho_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & rho, Real & drho_dp, Real & drho_dT, Real & drho_dx) const
{
    Real R_specific = R_molar / molarMass_from_X(x);
    rho = p / R_specific / T;
    drho_dp = 1 / R_specific / T;
    drho_dT = -p / R_specific / T / T;
    drho_dx = rho * molarMass_from_X(x) * ( 1/ _fp_primary->molarMass() - 1 / _fp_secondary[0]->molarMass() );

    /*rho =  rho_from_p_T_X(p, T, x);
    // Derivatives calculation
    Real xp = xp_from_X(x);
    Real M_g = molarMass_from_X(x);
    Real p_i = p * (xp * M_g / _fp_primary->molarMass() );
    Real rho_i = 0., drhoi_dT = 0., drhoi_dp = 0.;
    _fp_primary->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
    Real dv_dT = xp_from_X(x) * drhoi_dT / std::pow(rho_i, 2);
    Real dv_dp = xp_from_X(x) * drhoi_dp / std::pow(rho_i, 2);
    for (unsigned int i = 0; i < _n_secondary_gas; i++){
        p_i = p * (x[i] * M_g / _fp_secondary[i]->molarMass());
        _fp_secondary[i]->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
        dv_dT += x[i] * drhoi_dT / std::pow(rho_i, 2);
        dv_dp += x[i] * drhoi_dp / std::pow(rho_i, 2);
    }
    drho_dp = std::pow(rho, 2) * dv_dp;
    drho_dT = std::pow(rho, 2) * dv_dT;
    drho_dx = -std::pow(rho, 2) / rho_i;
    */
}

void
GasMixPHFluidProperties:: rho_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x, ADReal & rho, ADReal & drho_dp, ADReal & drho_dT, ADReal & drho_dx) const
{
    ADReal R_specific = R_molar / molarMass_from_X(x);
    rho = p / R_specific / T;
    drho_dp = 1 / R_specific / T;
    drho_dT = - p / R_specific / T / T;
    drho_dx = rho * molarMass_from_X(x) * ( 1/ _fp_primary->molarMass() - 1 / _fp_secondary[0]->molarMass() );

    /*rho =  rho_from_p_T_X(p, T, x);
    // Derivatives calculation
    ADReal xp = xp_from_X(x);
    ADReal M_g = molarMass_from_X(x);
    ADReal p_i = p * (xp * M_g / _fp_primary->molarMass() );
    ADReal rho_i = 0., drhoi_dT = 0., drhoi_dp = 0.;
    _fp_primary->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
    ADReal dv_dT = xp_from_X(x) * drhoi_dT / std::pow(rho_i, 2);
    ADReal dv_dp = xp_from_X(x) * drhoi_dp / std::pow(rho_i, 2);
    for (unsigned int i = 0; i < _n_secondary_gas; i++){
        p_i = p * (x[i] * M_g / _fp_secondary[i]->molarMass());
        _fp_secondary[i]->rho_from_p_T(p_i, T, rho_i, drhoi_dp, drhoi_dT);
        dv_dT += x[i] * drhoi_dT / std::pow(rho_i, 2);
        dv_dp += x[i] * drhoi_dp / std::pow(rho_i, 2);
    }
    drho_dp = std::pow(rho, 2) * dv_dp;
    drho_dT = std::pow(rho, 2) * dv_dT;
    drho_dx = -std::pow(rho, 2) / rho_i;
    */
}

Real
GasMixPHFluidProperties:: v_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
    return 1. / rho_from_p_T_X(p, T, x);
}

ADReal
GasMixPHFluidProperties:: v_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
    return 1. / rho_from_p_T_X(p, T, x);
}

void
GasMixPHFluidProperties:: v_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & v, Real & dv_dp, Real & dv_dT) const
{
    Real rho, drho_dp, drho_dT;
    rho_from_p_T_X(p, T, x, rho, drho_dp, drho_dT);
    v = 1. / rho;

    Real dv_drho = - 1 / (v * v);
    dv_dp = dv_drho * drho_dp;
    dv_dT = dv_drho * drho_dT;
}

Real
GasMixPHFluidProperties::e_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
    Real xp = xp_from_X(x);
    Real v = v_from_p_T_X(p, T, x);
    Real e = xp * _fp_primary->e_from_T_v(T, v / xp);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
        e += x[i] * _fp_secondary[i]->e_from_T_v(T, v / x[i]);
    return e;
}

ADReal
GasMixPHFluidProperties::e_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
    ADReal xp = xp_from_X(x);
    ADReal v = v_from_p_T_X(p, T, x);
    ADReal e = xp * _fp_primary->e_from_T_v(T, v / xp);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
        e += x[i] * _fp_secondary[i]->e_from_T_v(T, v / x[i]);
    return e;
}

void
GasMixPHFluidProperties::e_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & e, Real & de_dp, Real & de_dT) const
{
    e = e_from_p_T_X(p, T, x);
    de_dp = 0.0;
    de_dT = cv_from_p_T_X(p, T, x); // constant cv
}

void
GasMixPHFluidProperties::e_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & e, Real & de_dp, Real & de_dT, Real & de_dx) const
{
    Real v = v_from_p_T_X(p, T, x);
    Real xp = xp_from_X(x);
    e = e_from_p_T_X(p, T, x);
    de_dp = 0.0;
    de_dT = cv_from_p_T_X(p, T, x); // constant cv
    de_dx = _fp_secondary[0]->e_from_T_v(T, v / x[0]);
}

void
GasMixPHFluidProperties::e_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x, ADReal & e, ADReal & de_dp, ADReal & de_dT, ADReal & de_dx) const
{
    ADReal v = v_from_p_T_X(p, T, x);
    ADReal xp = xp_from_X(x);
    e = e_from_p_T_X(p, T, x);
    de_dp = 0.0;
    de_dT = cv_from_p_T_X(p, T, x); // constant cv
    de_dx = _fp_secondary[0]->e_from_T_v(T, v / x[0]);
}

Real
GasMixPHFluidProperties::e_from_p_rho_X(Real p, Real rho, const std::vector<Real> & x) const
{
    Real R_specific = R_molar / molarMass_from_X(x);
    Real T = p / rho / R_specific;
    Real cv = cv_from_p_T_X(p, T, x);
    return cv * T;
}

ADReal
GasMixPHFluidProperties::e_from_p_rho_X(ADReal p, ADReal rho, const std::vector<ADReal> & x) const
{
    ADReal R_specific = R_molar / molarMass_from_X(x);
    ADReal T =  p / rho / R_specific;
    ADReal cv = cv_from_p_T_X(p, T, x);
    return cv * T;
}


Real
GasMixPHFluidProperties:: h_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
    Real xp = xp_from_X(x);
    Real v = v_from_p_T_X(p, T, x);
    Real h = xp * _fp_primary->h_from_T_v(T, v / xp);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
        h += x[i] * _fp_secondary[i]->h_from_T_v(T, v / x[i]);
    return h;
}

ADReal
GasMixPHFluidProperties:: h_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
    ADReal xp = xp_from_X(x);
    ADReal v = v_from_p_T_X(p, T, x);
    ADReal h = xp * _fp_primary->h_from_T_v(T, v / xp);
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
        h += x[i] * _fp_secondary[i]->h_from_T_v(T, v / x[i]);
    return h;
}

void
GasMixPHFluidProperties:: h_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & h, Real & dh_dp, Real & dh_dT) const
{
    h = h_from_p_T_X(p, T, x);
    dh_dp = 0.;
    dh_dT = cp_from_p_T_X(p, T, x); // holds because cp is constant
}

Real GasMixPHFluidProperties::cp_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
    return h_from_p_T_X(p, T, x) / T;
}

ADReal GasMixPHFluidProperties::cp_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
    return h_from_p_T_X(p, T, x) / T;
}

void
GasMixPHFluidProperties::cp_from_p_T_X(Real p, Real T, const std::vector<Real> & x,Real & cp, Real & dcp_dp, Real & dcp_dT) const
{
    cp = cp_from_p_T_X(p, T, x);
    dcp_dp = 0.0;
    dcp_dT = 0.0;
}

Real GasMixPHFluidProperties::cv_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
    return e_from_p_T_X(p, T, x) / T;
}

ADReal GasMixPHFluidProperties::cv_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
    return e_from_p_T_X(p, T, x) / T;
}

void
GasMixPHFluidProperties::cv_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & cv, Real & dcv_dp, Real & dcv_dT) const
{
    cv = cv_from_p_T_X(p, T, x);
    dcv_dp = 0.0;
    dcv_dT = 0.0;
}

Real
GasMixPHFluidProperties::mu_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
  Real mu_g = 0.0;
  Real g_ij = 0.0, g_i = 0.0;
  Real mu_ij, M_ij;
  Real multiplier =  std::pow(2, 0.5) / 4;
  Real M_g = molarMass_from_X(x);

  std::vector<Real> mass_mol(_n_secondary_gas+1);
  std::vector<Real> yi(_n_secondary_gas+1);
  std::vector<Real> mui(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  yi[0] = xp_from_X(x) * M_g / mass_mol[0];
  mui[0] = _fp_primary->mu_from_p_T(p * yi[0], T);
  for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    mui[i] = _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T);
  }

  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            mu_g += mui[i] / (1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return mu_g;
}

ADReal
GasMixPHFluidProperties::mu_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
  ADReal mu_g = 0.0;
  ADReal g_ij = 0.0, g_i = 0.0;
  ADReal mu_ij, M_ij;
  ADReal multiplier =  std::pow(2, 0.5) / 4;
  ADReal M_g = molarMass_from_X(x);

  std::vector<ADReal> mass_mol(_n_secondary_gas+1);
  std::vector<ADReal> yi(_n_secondary_gas+1);
  std::vector<ADReal> mui(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  yi[0] = xp_from_X(x) * M_g / mass_mol[0];
  mui[0] = _fp_primary->mu_from_p_T(p * yi[0], T);
  for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    mui[i] = _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T);
  }

  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            mu_g += mui[i] / (1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return mu_g;
}

void
GasMixPHFluidProperties::mu_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
    // Value of mu
    mu = mu_from_p_T_X(p, T, x);
    // Calculate the derivatives
    Real g_ij = 0.0, g_i = 0.0;
    Real mu_ij, M_ij;
    Real dmui_p, dmui_T;
    Real dgij_dp, dgij_dT, dgi_dp, dgi_dT;
    Real dgp_i, dgT_i;
    Real dummy;
    Real multiplier= std::pow(2, 0.5) / 4;
    Real M_g = molarMass_from_X(x);

    std::vector<Real> mass_mol(_n_secondary_gas+1);
    std::vector<Real> yi(_n_secondary_gas+1);
    std::vector<Real> mui(_n_secondary_gas+1);
    std::vector<Real> dmui_dp(_n_secondary_gas+1);
    std::vector<Real> dmui_dT(_n_secondary_gas+1);

    mass_mol[0] = _fp_primary->molarMass();
    yi[0] = xp_from_X(x) * M_g / mass_mol[0];
    mui[0] = _fp_primary->mu_from_p_T(p * yi[0], T);
    _fp_primary->mu_from_p_T(p * yi[0], T, dummy, dmui_p, dmui_T);
    dmui_dp[0] = dmui_p;
    dmui_dT[0] = dmui_T;
    for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
        mass_mol[i] = _fp_secondary[i-1]->molarMass();
        yi[i] = x[i-1] * M_g / mass_mol[i];
        mui[i] = _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T);
        _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T, dummy, dmui_p, dmui_T);
        dmui_dp[i] = dmui_p;
        dmui_dT[i] = dmui_T;
    }

    for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            dgp_i = 0.0;
            dgT_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];

                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);

                    dgij_dp = multiplier * (1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)) / (std::pow((1 + M_ij), 0.5) *std::pow(mu_ij, 0.5) * std::pow(M_ij, 0.25));
                    dgij_dp = dgij_dp * (dmui_dp[i]* mui[j] - mui[i] * dmui_dp[j]) / mui[j] / mui[j];

                    dgij_dT = multiplier * (1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)) / (std::pow((1 + M_ij), 0.5) *std::pow(mu_ij, 0.5) * std::pow(M_ij, 0.25));
                    dgij_dT = dgij_dT * (dmui_dT[i]* mui[j] - mui[i] * dmui_dT[j]) / mui[j] / mui[j];

                    dgp_i += dgij_dp * yi[j];
                    dgT_i += dgij_dT * yi[j];
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            dmu_dp += dmuii_calc(mui[i], dmui_dp[i], g_i / yi[i], dgp_i / yi[i]);
            dmu_dT += dmuii_calc(mui[i], dmui_dT[i], g_i / yi[i], dgT_i / yi[i]);
        }
        else{
            continue;
        }
    }
}

Real
GasMixPHFluidProperties::k_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const
{
  Real k_g = 0.0;
  Real g_ij = 0.0, g_i = 0.0;
  Real mu_ij, M_ij;
  Real multiplier = std::pow(2, 0.5) / 4;
  Real M_g = molarMass_from_X(x);

  std::vector<Real> mass_mol(_n_secondary_gas+1);
  std::vector<Real> yi(_n_secondary_gas+1);
  std::vector<Real> mui(_n_secondary_gas+1);
  std::vector<Real> ki(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  yi[0] = xp_from_X(x) * M_g / mass_mol[0];
  mui[0] = _fp_primary->mu_from_p_T(p * yi[0], T);
  ki[0] = _fp_primary->k_from_p_T(p * yi[0], T);

  for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    mui[i] = _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T);
    ki[i] = _fp_secondary[i-1]->k_from_p_T(p * yi[i], T);
  }

  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }

            k_g += ki[i] /(1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return k_g;
}

ADReal
GasMixPHFluidProperties::k_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const
{
  ADReal k_g = 0.0;
  ADReal g_ij = 0.0, g_i = 0.0;
  ADReal mu_ij, M_ij;
  ADReal multiplier = std::pow(2, 0.5) / 4;
  ADReal M_g = molarMass_from_X(x);

  std::vector<ADReal> mass_mol(_n_secondary_gas+1);
  std::vector<ADReal> yi(_n_secondary_gas+1);
  std::vector<ADReal> mui(_n_secondary_gas+1);
  std::vector<ADReal> ki(_n_secondary_gas+1);

  mass_mol[0] = _fp_primary->molarMass();
  yi[0] = xp_from_X(x) * M_g / mass_mol[0];
  mui[0] = _fp_primary->mu_from_p_T(p * yi[0], T);
  ki[0] = _fp_primary->k_from_p_T(p * yi[0], T);

  for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
    mass_mol[i] = _fp_secondary[i-1]->molarMass();
    yi[i] = x[i-1] * M_g / mass_mol[i];
    mui[i] = _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T);
    ki[i] = _fp_secondary[i-1]->k_from_p_T(p * yi[i], T);
  }

  for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];
                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }

            k_g += ki[i] /(1 + g_i / yi[i]);
        }
        else{
            continue;
        }
  }
  return k_g;
}

void
GasMixPHFluidProperties::k_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & k, Real & dk_dp, Real & dk_dT) const
{
    k = k_from_p_T_X(p, T, x);
    // Derivatives
    Real g_ij = 0.0, g_i = 0.0;
    Real mu_ij, M_ij;
    Real dmui_p, dmui_T;
    Real dki_p, dki_T;
    Real dgij_dp, dgij_dT, dgi_dp, dgi_dT;
    Real dgp_i, dgT_i;
    Real dummy;
    Real multiplier=  std::pow(2, 0.5) / 4;
    Real M_g = molarMass_from_X(x);

    std::vector<Real> mass_mol(_n_secondary_gas+1);
    std::vector<Real> yi(_n_secondary_gas+1);
    std::vector<Real> mui(_n_secondary_gas+1);
    std::vector<Real> dmui_dp(_n_secondary_gas+1);
    std::vector<Real> dmui_dT(_n_secondary_gas+1);
    std::vector<Real> ki(_n_secondary_gas+1);
    std::vector<Real> dki_dp(_n_secondary_gas+1);
    std::vector<Real> dki_dT(_n_secondary_gas+1);

    mass_mol[0] = _fp_primary->molarMass();
    yi[0] = xp_from_X(x) * M_g / mass_mol[0];
    mui[0] = _fp_primary->mu_from_p_T(p * yi[0], T);
    _fp_primary->mu_from_p_T(p * yi[0], T, dummy, dmui_p, dmui_T);
    dmui_dp[0] = dmui_p;
    dmui_dT[0] = dmui_T;
    _fp_primary->k_from_p_T(p * yi[0], T, dummy, dki_p, dki_T);
    dki_dp[0] = dki_p;
    dki_dT[0] = dki_T;

    for (unsigned int i = 1; i < _n_secondary_gas+1; i++){
        mass_mol[i] = _fp_secondary[i-1]->molarMass();
        yi[i] = x[i-1] * M_g / mass_mol[i];
        mui[i] = _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T);
        _fp_secondary[i-1]->mu_from_p_T(p * yi[i], T, dummy, dmui_p, dmui_T);
        dmui_dp[i] = dmui_p;
        dmui_dT[i] = dmui_T;
        _fp_secondary[i-1]->k_from_p_T(p * yi[i], T, dummy, dki_p, dki_T);
        dki_dp[i] = dki_p;
        dki_dT[i] = dki_T;
    }

    for (unsigned int i = 0; i < _n_secondary_gas+1; i++){
        if (yi[i] != 0){
            g_i = 0.0;
            dgp_i = 0.0;
            dgT_i = 0.0;
            for (unsigned int j = 0; j < _n_secondary_gas+1; j++){
                if (j != i ){
                    mu_ij = mui[i] / mui[j];
                    M_ij = mass_mol[i] / mass_mol[j];

                    g_ij = multiplier / std::pow((1 + M_ij), 0.5) * std::pow((1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)),2);

                    dgij_dp = multiplier * (1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)) / (std::pow((1 + M_ij), 0.5) *std::pow(mu_ij, 0.5) * std::pow(M_ij, 0.25));
                    dgij_dp = dgij_dp * (dmui_dp[i]* mui[j] - mui[i] * dmui_dp[j]) / mui[j] / mui[j];

                    dgij_dT = multiplier * (1 + std::pow(mu_ij, 0.5) / std::pow(M_ij, 0.25)) / (std::pow((1 + M_ij), 0.5) *std::pow(mu_ij, 0.5) * std::pow(M_ij, 0.25));
                    dgij_dT = dgij_dT * (dmui_dT[i]* mui[j] - mui[i] * dmui_dT[j]) / mui[j] / mui[j];

                    dgp_i += dgij_dp * yi[j];
                    dgT_i += dgij_dT * yi[j];
                    g_i += g_ij * yi[j];
                }
                else{
                    continue;
                }
            }
            dk_dp += dmuii_calc(ki[i], dki_dp[i], g_i / yi[i], dgp_i / yi[i]);
            dk_dT += dmuii_calc(ki[i], dki_dT[i], g_i / yi[i], dgT_i / yi[i]);
        }
        else{
            continue;
        }
    }
}

Real
GasMixPHFluidProperties::dmuii_calc(Real mui, Real dmui, Real gi, Real dgi) const
{
  Real dmuii = (dmui * (1 + gi) - mui * dgi) / (1 + gi) / (1 + gi);
  return dmuii;
}

Real
GasMixPHFluidProperties::xp_from_X(const std::vector<Real> & x) const
{
    Real xp = 1.0;
    for (unsigned int i = 0; i < _n_secondary_gas; i++ ){
        xp -= x[i];
    }
    return xp;
}

ADReal
GasMixPHFluidProperties::xp_from_X(const std::vector<ADReal> & x) const
{
    ADReal xp = 1.0;
    for (unsigned int i = 0; i < _n_secondary_gas; i++ )
        xp -= x[i];
    return xp;
}

Real
GasMixPHFluidProperties::molarMass_from_X(const std::vector<Real> & x) const
{
    Real xp = xp_from_X(x);
    Real Mg = xp / _fp_primary->molarMass();
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
        Mg += x[i] / _fp_secondary[i]->molarMass();
    return 1 / Mg;
}

ADReal
GasMixPHFluidProperties::molarMass_from_X(const std::vector<ADReal> & x) const
{
    ADReal xp = xp_from_X(x);
    ADReal Mg = xp/_fp_primary->molarMass();
    for (unsigned int i = 0; i < _n_secondary_gas; i++)
        Mg += x[i] / _fp_secondary[i]->molarMass();
    return 1 / Mg;
}
