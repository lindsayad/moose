//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CGMassBC.h"
#include "Function.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", CGMassBC);

InputParameters
CGMassBC::validParams()
{
  InputParameters params = ADDGConvectionBC::validParams();
  // Note: we (arbitrarily) multiply the mass conservation term by -1 so that it matches the -p(div
  // v) term in the momentum equation.  This makes the matrix symmetric
  params.set<MaterialPropertyName>("advected_quantity") = "-1";
  params.suppressParameter<FunctionName>("primal_dirichlet_value");
  params.suppressParameter<MaterialPropertyName>(NS::density);
  params.suppressParameter<MaterialPropertyName>(NS::cp);
  return params;
}

CGMassBC::CGMassBC(const InputParameters & parameters) : ADDGConvectionBC(parameters) {}
