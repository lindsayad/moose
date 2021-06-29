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

#pragma once

#include "FluidProperties.h"
#include "SinglePhaseFluidProperties.h"
 /**
 * Class for computing the fluid properties of an arbitrary gas mixture
 * given the properties of the N mixture gaseous components and the mass
 * fractions.
 */
class GasMixPHFluidProperties : public SinglePhaseFluidProperties
{
public:
  static InputParameters validParams();
   ~GasMixPHFluidProperties();

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

   GasMixPHFluidProperties(const InputParameters & parameters);

  /**
   * Temperature from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @param[in] x   mass fractions
   * @return temperature (K)
   */
    virtual Real p_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const override;
    virtual ADReal p_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const override;
    virtual void p_from_v_e_X(Real v,  Real e, const std::vector<Real> & x, Real & p, Real & dp_dv, Real & dp_de, Real & dp_dx) const override;
    virtual void p_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x, ADReal & p, ADReal & dp_dv, ADReal & dp_de, ADReal & dp_dx) const override;

  /**
   * Temperature from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @param[in] x   mass fractions
   * @return temperature (K)
   */
    virtual Real T_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const override;
    virtual ADReal T_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const override;
    virtual void T_from_v_e_X(Real v, Real e,  const std::vector<Real> & x, Real & T, Real & dT_dv, Real & dT_de, Real & dT_dx) const override;
    virtual void T_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x, ADReal & T, ADReal & dT_dv, ADReal & dT_de, ADReal & dT_dx) const override;

    virtual Real c_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const override;
    virtual ADReal c_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const override;
    virtual void c_from_v_e_X(Real v, Real e, const std::vector<Real> & x, Real c, Real dc_dv, Real dc_de, Real dc_dx) const override;
    virtual void c_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x, ADReal c, ADReal dc_dv, ADReal dc_de, ADReal dc_dx) const override;

    virtual Real c_from_p_T_x(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal c_from_p_T_x(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;
  /**
   * Isobaric specific heat from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @param[in] x   mass fractions
   * @return isobaric specific heat (J/kg.K)
   */
    virtual Real cp_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const override;
    virtual ADReal cp_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const override;

  /**
   * Isochoric specific heat from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @param[in] x   mass fractions
   * @return isochoric specific heat (J/kg.K)
   */
    virtual Real cv_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const override;
    virtual ADReal cv_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const override;

   /**
   * Dynamic viscosity from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @param[in] x   mass fractions
   * @return dynamic viscosity (Pa.s)
   */
    virtual Real mu_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const override;
    virtual ADReal mu_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const override;
   /**
   * Thermal conductivity from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @param[in] x   mass fractions
   * @return thermal conductivity (W/m.K)
   */
    virtual Real k_from_v_e_X(Real v, Real e, const std::vector<Real> & x) const override;
    virtual ADReal k_from_v_e_X(ADReal v, ADReal e, const std::vector<ADReal> & x) const override;

  /**
   * Density from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @param[in] x   mass fractions
   * @return density (kg/m$^3$)
   */
    virtual Real beta_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual Real rho_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal rho_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

  /**
   * Density and its derivatives from pressure and temperature
   *
   * @param[in] p          pressure (Pa)
   * @param[in] T          temperature (K)
   * @param[in] x          mass fractions
   * @param[out] rho       density (kg/m$^3$)
   * @param[out] drho_dp   derivative of density w.r.t. pressure
   * @param[out] drho_dT   derivative of density w.r.t. temperature
   */
    virtual void rho_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & rho, Real & drho_dp, Real & drho_dT) const override;
    virtual void rho_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x, ADReal & rho, ADReal & drho_dp, ADReal & drho_dT) const override;
    virtual void rho_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & rho, Real & drho_dp, Real & drho_dT, Real & drho_dx) const override;
    virtual void rho_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x, ADReal & rho, ADReal & drho_dp, ADReal & drho_dT, ADReal & drho_dx) const override;
   /**
   * Specific volume from pressure and temperature
   * @param p    pressure (Pa)
   * @param T    temperature (K)
   * @param x    mass fractions
   * @return specific volume (m$^3$/kg)
   */
    virtual Real v_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal v_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

  /**
   * Specific volume and its derivatives from pressure and temperature
   * @param[in] p       pressure (Pa)
   * @param[in] T       temperature (K)
   * @param[in] x       mass fractions
   * @param[out] v      specific volume (m$^3$/kg)
   * @param[out] dv_dp  derivative of specific volume with respect to pressure
   * @param[out] dv_dT  derivative of specific volume with respect to temperature
   */
    virtual void v_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & v, Real & dv_dp, Real & dv_dT) const override;
  /**
   * Specific internal energy from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @param[in] x   mass fractions
   * @return specific internal energy (J/kg)
   */
    virtual Real e_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal e_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

  /**
   * Specific internal energy and its derivatives from pressure and temperature
   *
   * @param[in] p        pressure (Pa)
   * @param[in] T        temperature (K)
   * @param[in] x        mass fractions
   * @param[out] e       specific internal energy (J/kg)
   * @param[out] de_dp   derivative of specific internal energy w.r.t. pressure
   * @param[out] de_dT   derivative of specific internal energy w.r.t. temperature
   */
    virtual void e_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & e, Real & de_dp, Real & de_dT) const override;
    virtual void e_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & e, Real & de_dp, Real & de_dT, Real & de_dx) const override;
    virtual void e_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x, ADReal & e, ADReal & de_dp, ADReal & de_dT, ADReal & de_dx) const override;

    virtual Real e_from_p_rho_X(Real p, Real rho, const std::vector<Real> & x) const override;
    virtual ADReal e_from_p_rho_X(ADReal p, ADReal rho, const std::vector<ADReal> & x) const override;

  /**
   * Specific enthalpy from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @param[in] x   mass fractions
   * @return specific enthalpy (J/kg)
   */
    virtual Real h_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal h_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

  /**
   * Specific enthalpy and its derivatives from pressure and temperature
   *
   * @param[in] p        pressure (Pa)
   * @param[in] T        temperature (K)
   * @param[in] x        mass fractions
   * @param[out] h       specific enthalpy (J/kg)
   * @param[out] dh_dp   derivative of specific enthalpy w.r.t. pressure
   * @param[out] dh_dT   derivative of specific enthalpy w.r.t. temperature
   */
    virtual void h_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & h, Real & dh_dp, Real & dh_dT) const override;


  //ADReal beta_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

  /**
   * Isobaric specific heat capacity from pressure and temperature
   *
   * @param p   pressure (Pa)
   * @param T   temperature (K)
   * @param x   mass fractions
   * @return isobaric specific heat (J/kg/.K)
   */
    virtual Real cp_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal cp_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

  /**
   * Isobaric specific heat capacity and its derivatives from pressure and temperature
   *
   * @param[in] p       pressure (Pa)
   * @param[in] T       temperature (K)
   * @param[in] x           mass fractions
   * @param[out] cp     isobaric specific heat (J/kg/K)
   * @param[out] dcp_dp derivative of isobaric specific heat w.r.t. pressure (J/kg/K/Pa)
   */
    virtual void cp_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & cp, Real & dcp_dp, Real & dcp_dT) const override;

  /**
   * Isochoric specific heat capacity from pressure and temperature
   *
   * @param p   pressure (Pa)
   * @param T   temperature (K)
   * @param x   mass fractions
   * @return isochoric specific heat (J/kg.K)
   */
    virtual Real cv_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal cv_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

   /**
   * Isochoric specific heat capacity and its derivatives from pressure and temperature
   *
   * @param[in] p       pressure (Pa)
   * @param[in] T       temperature (K)
   * @param[in] x           mass fractions
   * @param[out] cv     isochoric specific heat (J/kg/K)
   * @param[out] dcv_dp derivative of isochoric specific heat w.r.t. pressure (J/kg/K/Pa)
   */
    virtual void cv_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & cv, Real & dcv_dp, Real & dcv_dT) const override;

    /**
     * Dynamic viscosity from pressure and temperature
     *
     * @param p   pressure (Pa)
     * @param T   temperature (K)
     * @param x   mass fractions
     * @return dynamic viscosity (Pa.s)
     */
    virtual Real mu_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
    virtual ADReal mu_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;

    /**
     * Dynamic viscosity and its derivatives wrt pressure and temperature
     *
     * @param p             pressure (Pa)
     * @param T             temperature (K)
     * @param x             mass fractions
     * @param[out] mu       viscosity (Pa.s)
     * @param[out] dmu_dp   derivative of viscosity wrt pressure
     * @param[out] dmu_dT   derivative of viscosity wrt temperature
     */
   virtual void mu_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & mu, Real & dmu_drho, Real & dmu_dT) const override;
   /**
   * Thermal conductivity from pressure and temperature
   *
   * @param p   pressure (Pa)
   * @param T   temperature (K)
   * @param x   mass fractions
   * @return thermal conductivity  (W/m.K)
   */
   virtual Real k_from_p_T_X(Real p, Real T, const std::vector<Real> & x) const override;
   virtual ADReal k_from_p_T_X(ADReal p, ADReal T, const std::vector<ADReal> & x) const override;
   /**
   * Thermal conductivity and its derivatives wrt pressure and temperature
   *
   * @param p           pressure (Pa)
   * @param T           temperature (K)
   * @param x           mass fractions
   * @param[out]  k     thermal conductivity  (W/m.K)
   * @param[out]  dk_dp derivative of thermal conductivity wrt pressure
   * @param[out]  dk_dT derivative of thermal conductivity wrt temperature
   */
    virtual void k_from_p_T_X(Real p, Real T, const std::vector<Real> & x, Real & k, Real & dk_dp, Real & dk_dT) const override;
   /**
   * Primary mass fractions from secondary gases mass fractions vectoe
   *
   * @param[in] x mass fractions
   * @return    mass fractions of primary species
   */
    virtual Real xp_from_X(const std::vector<Real> & x) const override;
    virtual ADReal xp_from_X(const std::vector<ADReal> & x) const override;
   /**
   * molar mass from the molar mass of the constituing gas species
   *
   * @param[in] x  mass fractions
   * @return    molar mass (kg/mol)
   */
   virtual Real molarMass_from_X(const std::vector<Real> & x) const override;
   virtual ADReal molarMass_from_X(const std::vector<ADReal> & x) const override;


  protected:
  /// Primary vapor fluid properties
  const SinglePhaseFluidProperties * const _fp_primary;
  /// Secondary vapor fluid properties
  std::vector<const SinglePhaseFluidProperties *> _fp_secondary;
  /// Names of secondary vapor fluid properties
  const std::vector<UserObjectName> _fp_secondary_names;
  /// Number of secondary vapors
  const unsigned int _n_secondary_gas;
  /// molar (or universal) gas constant
  static const Real R_molar;

  Real dmuii_calc(Real mui, Real dmui, Real gi, Real dgi) const;
};
