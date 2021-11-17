//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVMomentumAdvection.h"
#include "INSFVPressureVariable.h"
#include "PINSFVSuperficialVelocityVariable.h"
#include "FVUtils.h"
#include "MathFVUtils.h"
#include "NS.h"
#include "INSFVRhieChowInterpolator.h"

registerMooseObject("NavierStokesApp", PINSFVMomentumAdvection);

InputParameters
PINSFVMomentumAdvection::validParams()
{
  auto params = INSFVMomentumAdvection::validParams();
  params.addClassDescription("Object for advecting superficial momentum, e.g. rho*u_d, "
                             "in the porous media momentum equation");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity");
  params.addParam<bool>(
      "smooth_porosity", false, "Whether the porosity field is smooth or has discontinuities");

  return params;
}

PINSFVMomentumAdvection::PINSFVMomentumAdvection(const InputParameters & params)
  : INSFVMomentumAdvection(params),
    _eps(getFunctor<ADReal>(NS::porosity)),
    _smooth_porosity(getParam<bool>("smooth_porosity"))
{
  if (!dynamic_cast<const PINSFVSuperficialVelocityVariable *>(_u_var))
    mooseError("PINSFVMomentumAdvection may only be used with a superficial advective velocity, "
               "of variable type PINSFVSuperficialVelocityVariable.");
}

void
PINSFVMomentumAdvection::interpolate(Moose::FV::InterpMethod m, ADRealVectorValue & v)
{
  const Elem * const elem = &_face_info->elem();
  const Elem * const neighbor = _face_info->neighborPtr();

  if (onBoundary(*_face_info))
  {
    v(0) = _u_var->getBoundaryFaceValue(*_face_info);
    if (_v_var)
      v(1) = _v_var->getBoundaryFaceValue(*_face_info);
    if (_w_var)
      v(2) = _w_var->getBoundaryFaceValue(*_face_info);

    return;
  }

  const auto elem_face = elemFromFace();
  const auto neighbor_face = neighborFromFace();

  Moose::FV::interpolate(
      Moose::FV::InterpMethod::Average, v, _vel(elem_face), _vel(neighbor_face), *_face_info, true);

  if (m == Moose::FV::InterpMethod::Average)
    return;

  // Avoid computing a pressure gradient near porosity jumps
  if (!_smooth_porosity)
    if (MetaPhysicL::raw_value(_eps.gradient(elem)).norm() > 1e-12 ||
        MetaPhysicL::raw_value(_eps.gradient(neighbor)).norm() > 1e-12)
      return;

  mooseAssert((neighbor && this->hasBlocks(neighbor->subdomain_id())),
              "We should be on an internal face...");

  // Get pressure gradient. This is the uncorrected gradient plus a correction from cell centroid
  // values on either side of the face
  const VectorValue<ADReal> & grad_p = _p_var->adGradSln(*_face_info);

  // Get uncorrected pressure gradient. This will use the element centroid gradient if we are
  // along a boundary face
  const VectorValue<ADReal> & unc_grad_p = _p_var->uncorrectedAdGradSln(*_face_info);

  const Point & elem_centroid = _face_info->elemCentroid();
  const Point & neighbor_centroid = _face_info->neighborCentroid();
  Real elem_volume = _face_info->elemVolume();
  Real neighbor_volume = _face_info->neighborVolume();

  // Now we need to perform the computations of D
  const VectorValue<ADReal> & elem_a = _rc_uo.rcCoeff(elem);

  mooseAssert(_subproblem.getCoordSystem(elem->subdomain_id()) ==
                  _subproblem.getCoordSystem(neighbor->subdomain_id()),
              "Coordinate systems must be the same between the two elements");

  Real coord;
  coordTransformFactor(_subproblem, elem->subdomain_id(), elem_centroid, coord);

  elem_volume *= coord;

  VectorValue<ADReal> elem_D = 0;
  for (const auto i : make_range(_dim))
  {
    mooseAssert(elem_a(i).value() != 0, "We should not be dividing by zero");
    elem_D(i) = elem_volume / elem_a(i);
  }

  VectorValue<ADReal> face_D;

  const VectorValue<ADReal> & neighbor_a = _rc_uo.rcCoeff(neighbor);

  coordTransformFactor(_subproblem, neighbor->subdomain_id(), neighbor_centroid, coord);
  neighbor_volume *= coord;

  VectorValue<ADReal> neighbor_D = 0;
  for (const auto i : make_range(_dim))
  {
    mooseAssert(neighbor_a(i).value() != 0, "We should not be dividing by zero");
    neighbor_D(i) = neighbor_volume / neighbor_a(i);
  }
  Moose::FV::interpolate(
      Moose::FV::InterpMethod::Average, face_D, elem_D, neighbor_D, *_face_info, true);

  // evaluate face porosity, see (18) in Hanimann 2021 or (11) in Nordlund 2016
  const auto face_eps = _eps(Moose::FV::makeCDFace(*_face_info, faceArgSubdomains()));

  const auto & b1 = _rc_uo.getB1(*_face_info);
  const auto & b3 = _rc_uo.getB3(*_face_info);

  // perform the pressure correction
  for (const auto i : make_range(_dim))
    v(i) += -face_D(i) * face_eps * (grad_p(i) - unc_grad_p(i)) + face_D(i) * (b1(i) - b3(i));
}

ADReal
PINSFVMomentumAdvection::computeQpResidual()
{
  ADRealVectorValue v;
  ADReal adv_quant_interface;

  const auto elem_face = elemFromFace();
  const auto neighbor_face = neighborFromFace();

  // Superficial velocity interpolation
  this->interpolate(_velocity_interp_method, v);

  const auto interp_coeffs = Moose::FV::interpCoeffs(_advected_interp_method, *_face_info, true, v);

  _ae = _normal * v * _rho(elem_face) * interp_coeffs.first / _eps(elem_face);
  // Minus sign because we apply a minus sign to the residual in computeResidual
  _an = -_normal * v * _rho(neighbor_face) * interp_coeffs.second / _eps(neighbor_face);

  return _ae * _var(elem_face) - _an * _var(neighbor_face);
}
