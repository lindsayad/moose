//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVMomentumDiffusion.h"
#include "PINSFVSuperficialVelocityVariable.h"
#include "NS.h"
#include "INSFVRhieChowInterpolator.h"

registerMooseObject("NavierStokesApp", PINSFVMomentumDiffusion);

InputParameters
PINSFVMomentumDiffusion::validParams()
{
  auto params = INSFVMomentumDiffusion::validParams();
  params.addClassDescription("Viscous diffusion term, div(mu grad(u_d / eps)), in the porous media "
                             "incompressible Navier-Stokes momentum equation.");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity auxiliary variable");
  params.addParam<bool>(
      "smooth_porosity", false, "Whether to include the diffusion porosity gradient term");
  params.addParam<MooseFunctorName>(NS::superficial_velocity, "The superficial velocity");
  return params;
}

PINSFVMomentumDiffusion::PINSFVMomentumDiffusion(const InputParameters & params)
  : INSFVMomentumDiffusion(params),
    _eps(getFunctor<ADReal>(NS::porosity)),
    _vel(isParamValid(NS::superficial_velocity)
             ? &getFunctor<ADRealVectorValue>(NS::superficial_velocity)
             : nullptr),
    _smooth_porosity(getParam<bool>("smooth_porosity"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("PINSFV is not supported by local AD indexing. In order to use PINSFV, please run "
             "the configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif
  if (!dynamic_cast<PINSFVSuperficialVelocityVariable *>(&_var))
    mooseError("PINSFVMomentumDiffusion may only be used with a superficial velocity "
               "variable, of variable type PINSFVSuperficialVelocityVariable.");

  // Check that the parameters required for the porosity gradient term are set by the user
  if (_smooth_porosity && (!parameters().isParamSetByUser("momentum_component") ||
                           !isParamValid(NS::superficial_velocity)))
    paramError("smooth_porosity",
               "The porosity gradient diffusion term requires specifying "
               "both the momentum component and a superficial velocity material property.");
}

ADReal
PINSFVMomentumDiffusion::computeQpResidual()
{
#ifdef MOOSE_GLOBAL_AD_INDEXING
  using namespace Moose::FV;

  const auto elem_face = elemFromFace();
  const auto neighbor_face = neighborFromFace();
  const auto mu_elem = _mu(elem_face);
  const auto mu_neighbor = _mu(neighbor_face);
  const auto eps_elem = _eps(elem_face);
  const auto eps_neighbor = _eps(neighbor_face);

  // Compute the diffusion driven by the velocity gradient
  // Interpolate viscosity divided by porosity on the face
  ADReal mu_eps_face;
  interpolate(Moose::FV::InterpMethod::Average,
              mu_eps_face,
              mu_elem / eps_elem,
              mu_neighbor / eps_neighbor,
              *_face_info,
              true);

  // Compute face superficial velocity gradient
  auto dudn =
      _var.gradient(Moose::FV::makeCDFace(*_face_info, faceArgSubdomains())) * _face_info->normal();

  if (_computing_rc_data)
  {
    if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
        _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    {
      const auto dof_number = _face_info->elem().dof_number(_sys.number(), _var.number(), 0);
      // A gradient is a linear combination of degrees of freedom so it's safe to straight-up index
      // into the derivatives vector at the dof we care about
      _ae = dudn.derivatives()[dof_number];
      _ae *= -mu_eps_face;
    }
    if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
        _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    {
      const auto dof_number = _face_info->neighbor().dof_number(_sys.number(), _var.number(), 0);
      _an = dudn.derivatives()[dof_number];
      _an *= mu_eps_face;
    }
  }

  // First term of residual
  ADReal residual = mu_eps_face * dudn;

  if (_smooth_porosity)
  {
    // Get the face porosity gradient separately
    const auto & grad_eps_face =
        MetaPhysicL::raw_value(_eps.gradient(Moose::FV::makeCDFace(*_face_info)));

    ADRealVectorValue term_elem = mu_elem / eps_elem / eps_elem * grad_eps_face;
    ADRealVectorValue term_neighbor = mu_neighbor / eps_neighbor / eps_neighbor * grad_eps_face;

    const auto vel_elem = (*_vel)(elem_face);
    const auto vel_neighbor = (*_vel)(neighbor_face);

    for (int i = 0; i < LIBMESH_DIM; i++)
    {
      term_elem(i) *= vel_elem(i);
      term_neighbor(i) *= vel_neighbor(i);
    }

    // Interpolate to get the face value
    ADRealVectorValue term_face;
    interpolate(
        Moose::FV::InterpMethod::Average, term_face, term_elem, term_neighbor, *_face_info, true);
    residual -= term_face * _normal;
  }
  return -residual;
#else
  return 0;
#endif
}
