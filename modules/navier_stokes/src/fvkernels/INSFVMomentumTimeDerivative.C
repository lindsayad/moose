//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumTimeDerivative.h"

#include "NS.h"

registerMooseObject("NavierStokesApp", INSFVMomentumTimeDerivative);

InputParameters
INSFVMomentumTimeDerivative::validParams()
{
  InputParameters params = FVTimeKernel::validParams();
  params.addClassDescription(
      "Adds the time derivative term to the incompressible Navier-Stokes momentum equation.");
  params.addRequiredParam<Real>(NS::density, "The value for the density");
  params.declareControllable(NS::density);
  return params;
}

INSFVMomentumTimeDerivative::INSFVMomentumTimeDerivative(const InputParameters & params)
  : FVTimeKernel(params), INSFVMomentumResidualObject(*this), _rho(getParam<Real>(NS::density))
{
}

ADReal
INSFVMomentumTimeDerivative::computeQpResidual(const Elem & elem)
{
  return _rho * _var.dot(&elem);
}

ADReal
INSFVMomentumTimeDerivative::computeQpResidual()
{
  return computeQpResidual(*_current_elem);
}

void
INSFVMomentumTimeDerivative::gatherRCData(const Elem & elem)
{
  const auto saved_do_derivatives = ADReal::do_derivatives;
  // We rely on derivative indexing
  ADReal::do_derivatives = true;
  const auto residual = computeQpResidual(elem) * _assembly.elementVolume(&elem);
  const auto dof_number = elem.dof_number(_sys.number(), _var.number(), 0);
  const Real a = residual.derivatives()[dof_number];
  ADReal::do_derivatives = saved_do_derivatives;

  _rc_uo.addToA(&elem, _index, a);
}
