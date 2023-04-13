//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGConvectionBC.h"
#include "NS.h"
#include "Function.h"

registerMooseObject("NavierStokesApp", ADDGConvectionBC);

InputParameters
ADDGConvectionBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addParam<MaterialPropertyName>(
      "velocity_mat_prop",
      "Velocity vector as a material property. Should be provided when we want the velocity value "
      "to be determined implicitly (e.g. we don't have a Dirichlet condition)");
  params.addParam<FunctionName>("velocity_function",
                                "Function describing the values of velocity on the boundary.");
  params.addClassDescription("DG for convection");
  params.addParam<MaterialPropertyName>("advected_quantity",
                                        "An optional material property to be advected. If not "
                                        "supplied, then the variable will be used.");
  params.addParam<FunctionName>("primal_dirichlet_value",
                                "The value of the primal variable on the boundary.");
  params.addParam<MaterialPropertyName>(
      NS::density,
      "The density for potential use with the 'primal_dirichlet_value' parameter when advecting "
      "momentum or energy");
  params.addParam<MaterialPropertyName>(NS::cp,
                                        "The specific heat capcity for potential use with the "
                                        "'primal_dirichlet_value' parameter when advecting energy");
  return params;
}

ADDGConvectionBC::ADDGConvectionBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _velocity_mat_prop(isParamValid("velocity_mat_prop")
                           ? &getADMaterialProperty<RealVectorValue>("velocity_mat_prop")
                           : nullptr),
    _velocity_function(isParamValid("velocity_function") ? &getFunction("velocity_function")
                                                         : nullptr),
    _adv_quant(isParamValid("advected_quantity")
                   ? getADMaterialProperty<Real>("advected_quantity").get()
                   : _u),
    _primal_dirichlet(
        isParamValid("primal_dirichlet_value") ? &getFunction("primal_dirichlet_value") : nullptr),
    _rho(isParamValid(NS::density) ? &getADMaterialProperty<Real>(NS::density) : nullptr),
    _cp(isParamValid(NS::cp) ? &getADMaterialProperty<Real>(NS::cp) : nullptr)
{
  if (_rho && !_primal_dirichlet)
    paramError(NS::density,
               "Density should only be provided when 'primal_dirichlet_value' is provided");
  if (_cp && (!_primal_dirichlet || !_rho))
    paramError(NS::cp,
               "Specific heat capacity should only be provided when 'primal_dirichlet_value' and "
               "'density' are provided");
  if (static_cast<bool>(_primal_dirichlet) + isParamValid("advected_quantity") != 1)
    mooseError("Exactly one of 'primal_dirichlet_value' or 'advected_quantity' should be provided");
  if (static_cast<bool>(_velocity_mat_prop) + static_cast<bool>(_velocity_function) != 1)
    mooseError("Exactly one of 'velocity_mat_prop' or 'velocity_function' should be provided");
}

ADReal
ADDGConvectionBC::computeQpResidual()
{
  const auto vdotn =
      (_velocity_mat_prop ? (*_velocity_mat_prop)[_qp]
                          : ADRealVectorValue(_velocity_function->vectorValue(_t, _q_point[_qp]))) *
      _normals[_qp];

  if (_primal_dirichlet)
  {
    ADReal adv_quant = _primal_dirichlet->value(_t, _q_point[_qp]);
    if (_rho)
      adv_quant *= (*_rho)[_qp];
    if (_cp)
      adv_quant *= (*_cp)[_qp];
    return _test[_i][_qp] * vdotn * adv_quant;
  }
  else
    return _test[_i][_qp] * vdotn * _adv_quant[_qp];
}
