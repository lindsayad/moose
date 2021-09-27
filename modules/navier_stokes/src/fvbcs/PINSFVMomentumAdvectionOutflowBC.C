//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVMomentumAdvectionOutflowBC.h"
#include "PINSFVSuperficialVelocityVariable.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", PINSFVMomentumAdvectionOutflowBC);

InputParameters
PINSFVMomentumAdvectionOutflowBC::validParams()
{
  InputParameters params = INSFVMomentumAdvectionOutflowBC::validParams();
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity auxiliary variable");
  params.addClassDescription("Outflow boundary condition for advecting momentum. This will impose "
                             "a zero normal gradient on the boundary velocity.");
  return params;
}

PINSFVMomentumAdvectionOutflowBC::PINSFVMomentumAdvectionOutflowBC(const InputParameters & params)
  : INSFVMomentumAdvectionOutflowBC(params), _eps(getFunctor<ADReal>("porosity"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("PINSFV is not supported by local AD indexing. In order to use PINSFV, please run the "
             "configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif

  if (!dynamic_cast<const PINSFVSuperficialVelocityVariable *>(_u_var))
    paramError("u", "the u velocity must be a PINSFVSuperficialVelocityVariable.");

  if (_dim >= 2 && !dynamic_cast<const PINSFVSuperficialVelocityVariable *>(_v_var))
    paramError("v",
               "In two or more dimensions, the v velocity must be supplied and it must be an "
               "PINSFVSuperficialVelocityVariable.");

  if (_dim >= 3 && !dynamic_cast<const PINSFVSuperficialVelocityVariable *>(_w_var))
    paramError("w",
               "In three-dimensions, the w velocity must be supplied and it must be an "
               "PINSFVSuperficialVelocityVariable.");
}

ADReal
PINSFVMomentumAdvectionOutflowBC::computeQpResidual()
{
#ifdef MOOSE_GLOBAL_AD_INDEXING
  using namespace Moose::FV;

  const auto elem_face = elemFromFace();
  const auto neighbor_face = neighborFromFace();

  ADRealVectorValue v(_u_var->getBoundaryFaceValue(*_face_info));
  if (_v_var)
    v(1) = _v_var->getBoundaryFaceValue(*_face_info);
  if (_w_var)
    v(2) = _w_var->getBoundaryFaceValue(*_face_info);

  ADReal adv_quant_boundary;
  interpolate(_advected_interp_method,
              adv_quant_boundary,
              _adv_quant(elem_face) / _eps(elem_face),
              _adv_quant(neighbor_face) / _eps(neighbor_face),
              v,
              *_face_info,
              true);
  mooseAssert(_normal * v >= 0,
              "This boundary condition is for outflow but the flow is in the opposite direction of "
              "the boundary normal");
  return _normal * v * adv_quant_boundary;
#else
  mooseError("PINSFV is not supported by local AD indexing. In order to use PINSFV, please run the "
             "configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif
}

void
PINSFVMomentumAdvectionOutflowBC::gatherRCData(const FaceInfo & /*fi*/)
{
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
  // const auto face_rho = _rho(std::make_tuple(
  //     &fi, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id()));
  // const auto face_eps = _eps(std::make_tuple(
  //     &fi, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id()));

  // ADRealVectorValue face_velocity(_u_var->getBoundaryFaceValue(fi));
  // if (_v_var)
  //   face_velocity(1) = _v_var->getBoundaryFaceValue(fi);
  // if (_w_var)
  //   face_velocity(2) = _w_var->getBoundaryFaceValue(fi);

  // const auto advection_coeffs =
  //     Moose::FV::interpCoeffs(_advected_interp_method, fi, true, face_velocity);
  // const auto advection_coeff = var_on_elem_side ? advection_coeffs.first :
  // advection_coeffs.second; const ADReal coeff =
  //     face_rho * face_velocity / face_eps * surface_vector * advection_coeff * residual_sign;

  // _rc_uo.addToA(boundary_elem, _index, coeff);
}
