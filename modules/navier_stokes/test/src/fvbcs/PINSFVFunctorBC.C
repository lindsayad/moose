//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVFunctorBC.h"
#include "NS.h"
#include "SinglePhaseFluidProperties.h"
#include "Function.h"

registerMooseObject("NavierStokesTestApp", PINSFVFunctorBC);

InputParameters
PINSFVFunctorBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params += INSFVFlowBC::validParams();
  params += INSFVMomentumResidualObject::validParams();
  params.makeParamNotRequired<MooseEnum>("momentum_component");
  params.addClassDescription("Computes the residual of the advective and pressure term (the latter "
                             "when this object is added for the momentum equation) on a boundary.");
  params.addParam<MaterialPropertyName>(NS::density, NS::density, "The density material property");
  MooseEnum eqn("mass momentum");
  params.addRequiredParam<MooseEnum>("eqn", eqn, "The equation you're solving.");
  params.addRequiredCoupledVar(NS::superficial_velocity_x,
                               "The x-component of the superficial velocity");
  params.addCoupledVar(NS::superficial_velocity_y, "The y-component of the superficial velocity");
  params.addCoupledVar(NS::superficial_velocity_z, "The z-component of the superficial velocity");
  params.addRequiredCoupledVar(NS::pressure, "The pressure");
  params.addRequiredCoupledVar(NS::porosity, "The porosity");
  MooseEnum advected_interp_method("average upwind", "upwind");
  params.addParam<MooseEnum>("advected_interp_method",
                             advected_interp_method,
                             "The interpolation to use for the advected quantity. Options are "
                             "'upwind' and 'average', with the default being 'upwind'.");
  return params;
}

PINSFVFunctorBC::PINSFVFunctorBC(const InputParameters & params)
  : FVFluxBC(params),
    INSFVFlowBC(params),
    INSFVMomentumResidualObject(*this),
    _sup_vel_x(getFunctor<ADReal>(NS::superficial_velocity_x)),
    _sup_vel_y(isParamValid(NS::superficial_velocity_y)
                   ? &getFunctor<ADReal>(NS::superficial_velocity_y)
                   : nullptr),
    _sup_vel_z(isParamValid(NS::superficial_velocity_z)
                   ? &getFunctor<ADReal>(NS::superficial_velocity_z)
                   : nullptr),
    _pressure(getFunctor<ADReal>(NS::pressure)),
    _rho(getFunctor<ADReal>(NS::density)),
    _eps(getFunctor<ADReal>(NS::porosity)),
    _eqn(getParam<MooseEnum>("eqn")),
    _index(getParam<MooseEnum>("momentum_component"))
{
  if ((_eqn == "momentum") && !isParamValid("momentum_component"))
    paramError("eqn",
               "If 'momentum' is specified for 'eqn', then you must provide a parameter "
               "value for 'momentum_component'");
  if ((_eqn != "momentum") && isParamValid("momentum_component"))
    paramError("momentum_component",
               "'momentum_component' should not be specified when the 'eqn' is not 'momentum'");

  using namespace Moose::FV;

  const auto & advected_interp_method = getParam<MooseEnum>("advected_interp_method");
  if (advected_interp_method == "average")
    _advected_interp_method = InterpMethod::Average;
  else if (advected_interp_method == "upwind")
    _advected_interp_method = InterpMethod::Upwind;
  else
    mooseError("Unrecognized interpolation type ",
               static_cast<std::string>(advected_interp_method));
}

ADReal
PINSFVFunctorBC::computeQpResidual()
{
  const auto ft = _face_info->faceType(_var.name());
  const bool out_of_elem = (ft == FaceInfo::VarFaceNeighbors::ELEM);
  const auto normal = out_of_elem ? _face_info->normal() : Point(-_face_info->normal());

  // No interpolation on a boundary so argument values to fi_elem_is_upwind do not
  // matter
  const auto face = std::make_tuple(
      _face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains());
  const VectorValue<ADReal> sup_vel(_sup_vel_x(face),
                                    _sup_vel_y ? (*_sup_vel_y)(face) : ADReal(0),
                                    _sup_vel_z ? (*_sup_vel_z)(face) : ADReal(0));
  const auto rho = _rho(face);

  if (_eqn == "mass")
    return rho * sup_vel * normal;
  else if (_eqn == "momentum")
  {
    const auto eps = _eps(face);
    const auto rhou = sup_vel(_index) / eps * rho;
    return rhou * sup_vel * normal + eps * _pressure(face) * normal(_index);
  }
  else
    mooseError("Unrecognized equation type ", _eqn);
}

void
PINSFVFunctorBC::gatherRCData(const FaceInfo & /*fi*/)
{
  // This is copied from PINSFVMomentumAdvection for flow boundaries. In the future we should
  // probably make this look more like the above computeQpResidual

  // const Elem & elem = fi.elem();
  // const Elem * const neighbor = fi.neighborPtr();
  // const Point & normal = fi.normal();
  // Real coord;
  // coordTransformFactor(_subproblem, elem.subdomain_id(), fi.faceCentroid(), coord);
  // const auto surface_vector = normal * fi.faceArea() * coord;

  // auto ft = fi.faceType(_var.name());
  // mooseAssert((ft == FaceInfo::VarFaceNeighbors::ELEM) ||
  //                 (ft == FaceInfo::VarFaceNeighbors::NEIGHBOR),
  //             "Do we want to allow internal flow boundaries?");
  // const bool var_on_elem_side = ft == FaceInfo::VarFaceNeighbors::ELEM;
  // Real residual_sign = var_on_elem_side ? 1. : -1.;
  // const Elem * const boundary_elem = var_on_elem_side ? &elem : neighbor;

  // mooseAssert(boundary_elem, "the boundary elem should be non-null");
  // const auto boundary_face = std::make_tuple(
  //     &fi, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id());
  // const auto face_rho = _rho(boundary_face);
  // const auto face_eps = _eps(boundary_face);

  // ADRealVectorValue face_velocity(_sup_vel_x(boundary_face));
  // if (_sup_vel_y)
  //   face_velocity(1) = (*_sup_vel_y)(boundary_face);
  // if (_sup_vel_z)
  //   face_velocity(2) = (*_sup_vel_z)(boundary_face);

  // const auto advection_coeffs =
  //     Moose::FV::interpCoeffs(_advected_interp_method, fi, true, face_velocity);
  // const auto advection_coeff = var_on_elem_side ? advection_coeffs.first :
  // advection_coeffs.second; const ADReal coeff =
  //     face_rho * face_velocity / face_eps * surface_vector * advection_coeff * residual_sign;

  // _rc_uo.addToA(boundary_elem, _index, coeff);
}
