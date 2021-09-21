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
  params.addRequiredCoupledVar(NS::porosity, "Porosity auxiliary variable");
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
    _eps_var(dynamic_cast<const MooseVariableFVReal *>(getFieldVar(NS::porosity, 0))),
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

void
PINSFVMomentumDiffusion::gatherRCData(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;

  const Elem & elem = fi.elem();
  const Elem * const neighbor = fi.neighborPtr();
  const Point & normal = fi.normal();
  Real coord;
  coordTransformFactor(_subproblem, elem.subdomain_id(), fi.faceCentroid(), coord);
  const auto surface_vector = normal * fi.faceArea() * coord;

  if (onBoundary(fi))
  {
    auto ft = fi.faceType(_var.name());
    mooseAssert((ft == FaceInfo::VarFaceNeighbors::ELEM) ||
                    (ft == FaceInfo::VarFaceNeighbors::NEIGHBOR),
                "Do we want to allow internal boundaries?");
    const bool var_on_elem_side = ft == FaceInfo::VarFaceNeighbors::ELEM;
    const Elem * const boundary_elem = var_on_elem_side ? &elem : neighbor;
    const Point & boundary_elem_centroid =
        var_on_elem_side ? fi.elemCentroid() : fi.neighborCentroid();

    mooseAssert(boundary_elem, "the boundary elem should be non-null");
    mooseAssert(this->hasBlocks(boundary_elem->subdomain_id()),
                "This object should exist on our boundary element subdomain ID");
    const auto face_mu = _mu(std::make_tuple(
        &fi, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id()));
    const auto face_eps = _eps(std::make_tuple(
        &fi, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id()));

    // Find the boundary id that has an associated INSFV boundary condition
    // if a face has more than one bc_id
    for (const auto bc_id : fi.boundaryIDs())
    {
      if (_no_slip_wall_boundaries.find(bc_id) != _no_slip_wall_boundaries.end())
      {
        // Need to account for viscous shear stress from wall
        const ADReal coeff = face_mu / face_eps * surface_vector.norm() /
                             std::abs((fi.faceCentroid() - boundary_elem_centroid) * normal) *
                             (1 - normal(_index) * normal(_index));
        _rc_uo.addToA(boundary_elem, _index, coeff);
        return;
      }

      if (_flow_boundaries.find(bc_id) != _flow_boundaries.end())
      {
        if (_fully_developed_flow_boundaries.find(bc_id) == _fully_developed_flow_boundaries.end())
        {
          // We are not on a fully developed flow boundary, so we have a viscous term
          // contribution. This term is slightly modified relative to the internal face term.
          // Instead of the distance between elem and neighbor centroid, we just have the distance
          // between the elem and face centroid. Specifically, the term below is the result of
          // Moukalled 8.80, 8.82, and the orthogonal correction approach equation for E_f,
          // equation 8.89. So relative to the internal face viscous term, we have substituted
          // eqn. 8.82 for 8.78
          // Note: If mu is an effective diffusivity, this should not be divided by face_eps
          const ADReal coeff = face_mu / face_eps * surface_vector.norm() /
                               (fi.faceCentroid() - boundary_elem_centroid).norm();
          _rc_uo.addToA(boundary_elem, _index, coeff);
        }
        return;
      }

      if (_slip_wall_boundaries.find(bc_id) != _slip_wall_boundaries.end())
        mooseError("Slip wall boundaries should have a flux bc such that we should never get here");

      if (_symmetry_boundaries.find(bc_id) != _symmetry_boundaries.end())
        mooseError("Symmetry boundaries should have a flux bc such that we should never get here");
    }

    const auto bc_id = *fi.boundaryIDs().begin();
    mooseError("The INSFVMomentumAdvection object ",
               this->name(),
               " is not completely bounded by INSFVBCs. Please examine surface ",
               bc_id,
               " and your FVBCs blocks.");
  }

  // Else we are on an internal face

  const auto face_mu = _mu(std::make_tuple(
      &fi, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains(&fi)));
  const auto face_eps = _eps(std::make_tuple(
      &fi, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains(&fi)));

  // Now add the viscous flux. Note that this includes only the orthogonal component! See
  // Moukalled equations 8.80, 8.78, and the orthogonal correction approach equation for
  // E_f, equation 8.69
  const ADReal coeff = face_mu / face_eps * surface_vector.norm() /
                       (fi.neighborCentroid() - fi.elemCentroid()).norm();
  _rc_uo.addToA(&elem, _index, coeff);
  _rc_uo.addToA(neighbor, _index, coeff);
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
  auto dudn = gradUDotNormal();

  // First term of residual
  ADReal residual = mu_eps_face * dudn;

  if (_smooth_porosity)
  {
    // Get the face porosity gradient separately
    const auto & grad_eps_face = MetaPhysicL::raw_value(_eps_var->adGradSln(*_face_info));

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
