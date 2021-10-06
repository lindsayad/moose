//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVMatAdvectionOutflowBC.h"
#include "INSFVFullyDevelopedFlowBC.h"
#include "INSFVMomentumResidualObject.h"

class INSFVVelocityVariable;

/**
 * A class for finite volume fully developed outflow boundary conditions for the momentum equation
 * It advects momentum at the outflow, and may replace outlet pressure boundary conditions
 * when selecting a mean-pressure approach.
 */
class INSFVMomentumAdvectionOutflowBC : public FVMatAdvectionOutflowBC,
                                        public INSFVFullyDevelopedFlowBC,
                                        public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  INSFVMomentumAdvectionOutflowBC(const InputParameters & params);

  void gatherRCData(const Elem &) override {}
  void gatherRCData(const FaceInfo &) override;

protected:
  virtual ADReal computeQpResidual() override;

  /// x-velocity
  const INSFVVelocityVariable * const _u_var;
  /// y-velocity
  const INSFVVelocityVariable * const _v_var;
  /// z-velocity
  const INSFVVelocityVariable * const _w_var;

  /// the dimension of the simulation
  const unsigned int _dim;

  /// The density
  const Moose::Functor<ADReal> & _rho;

  /// whether we are computing Rhie-Chow data
  bool _computing_rc_data;

  /// A member to hold computation of the Rhie-Chow coefficient
  ADReal _a = 0;
};
