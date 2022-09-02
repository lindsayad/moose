//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "InterfaceKernel.h"
#include "JvarMapInterface.h"
using namespace std;
//c
//c This kernel implements DG (IIPG) method for single phase and fully saturated Darcy flow
//c
class PorousFlowFullySaturatedDarcyDG: public JvarMapKernelInterface<InterfaceKernel>
{
public:
  static InputParameters validParams();
  PorousFlowFullySaturatedDarcyDG(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;
  Real computeAverageFlux();
  Real computePorePressureJump();
  Real computeAverageConduct();

  //c
  //c  stabilized parameter (/lambda)
  //c
  Real _stabilized_para;

  //c
  //c "_L" and "_R" indicates variables defined for the left and
  //c the right side of the interface, respectively. _L is the primary
  //c side and _R is the associate side.
  //c
  /**
   * The mobility of the fluid = density / viscosity
   * from the left-side and right-side
   */

  virtual Real mobility_L() const;
  virtual Real mobility_R() const;

  /// If true then the mobility contains the fluid density, otherwise it doesn't
  const bool _multiply_by_density;

  /// Permeability of porous material
  const MaterialProperty<RealTensorValue> & _permeability_L;
  const MaterialProperty<RealTensorValue> & _permeability_R;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _density_L;
  const MaterialProperty<std::vector<Real>> & _density_R;

  /// Viscosity of the fluid at the qp
  const MaterialProperty<std::vector<Real>> & _viscosity_L;
  const MaterialProperty<std::vector<Real>> & _viscosity_R;

  /// Quadpoint pore pressure in each phase
  const MaterialProperty<std::vector<Real>> & _pp_L;
  const MaterialProperty<std::vector<Real>> & _pp_R;

  /// Gradient of the pore pressure in each phase
  const MaterialProperty<std::vector<RealGradient>> & _grad_p_L;
  const MaterialProperty<std::vector<RealGradient>> & _grad_p_R;

};
