//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumDiffusion.h"
#include "INSFVRhieChowInterpolator.h"

registerMooseObject("NavierStokesApp", INSFVMomentumDiffusion);

InputParameters
INSFVMomentumDiffusion::validParams()
{
  auto params = FVFluxKernel::validParams();
  params += INSFVResidualObject::validParams();
  params.addRequiredParam<MooseFunctorName>(NS::mu, "The viscosity");
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

INSFVMomentumDiffusion::INSFVMomentumDiffusion(const InputParameters & params)
  : FVFluxKernel(params), INSFVResidualObject(*this), _mu(getFunctor<ADReal>(NS::mu))
{
}

void
INSFVMomentumDiffusion::gatherRCData(const FaceInfo & fi)
{
  const Elem & elem = fi.elem();
  const Elem * const neighbor = fi.neighborPtr();
  const Point & normal = fi.normal();
  Real coord;
  coordTransformFactor(_subproblem, elem.subdomain_id(), fi.faceCentroid(), coord);
  const auto surface_vector = normal * fi.faceArea() * coord;

  const auto face_mu = _mu(std::make_tuple(
      &fi, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains(&fi)));

  if (onBoundary(fi))
  {
    auto ft = fi.faceType(_var.name());
    const Elem * const boundary_elem = (ft == FaceInfo::VarFaceNeighbors::ELEM) ? &elem : neighbor;
    const Point & boundary_elem_centroid =
        (ft == FaceInfo::VarFaceNeighbors::ELEM) ? fi.elemCentroid() : fi.neighborCentroid();

    // Find the boundary id that has an associated INSFV boundary condition
    // if a face has more than one bc_id
    for (const auto bc_id : fi.boundaryIDs())
    {
      if (_no_slip_wall_boundaries.find(bc_id) != _no_slip_wall_boundaries.end())
      {
        // Need to account for viscous shear stress from wall
        const ADReal coeff = face_mu * surface_vector.norm() /
                             std::abs((fi.faceCentroid() - boundary_elem_centroid) * normal) *
                             (1 - normal(_index) * normal(_index));
        _rc_uo.addToA(boundary_elem, _index, coeff);
        return;
      }

      if (_slip_wall_boundaries.find(bc_id) != _slip_wall_boundaries.end())
        // In the case of a slip wall we do not have viscous shear stress from the wall, so our
        // contribution to the coefficient is zero
        return;

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
          const ADReal coeff =
              face_mu * surface_vector.norm() / (fi.faceCentroid() - boundary_elem_centroid).norm();
          _rc_uo.addToA(boundary_elem, _index, coeff);
        }
        return;
      }

      if (_symmetry_boundaries.find(bc_id) != _symmetry_boundaries.end())
      {
        // Moukalled eqns. 15.154 - 15.156
        const ADReal coeff = 2. * face_mu * surface_vector.norm() /
                             std::abs((fi.faceCentroid() - boundary_elem_centroid) * normal) *
                             normal(_index) * normal(_index);
        _rc_uo.addToA(boundary_elem, _index, coeff);
        return;
      }
    }

    mooseError("The INSFVMomentumAdvection object ",
               this->name(),
               " is not completely bounded by INSFVBCs. Please examine sideset ",
               *fi.boundaryIDs().begin(),
               " and your FVBCs blocks.");
  }

  // Else we are on an internal face

  // Now add the viscous flux. Note that this includes only the orthogonal component! See
  // Moukalled equations 8.80, 8.78, and the orthogonal correction approach equation for
  // E_f, equation 8.69
  const ADReal coeff =
      face_mu * surface_vector.norm() / (fi.neighborCentroid() - fi.elemCentroid()).norm();
  _rc_uo.addToA(&elem, _index, coeff);
  _rc_uo.addToA(neighbor, _index, coeff);
}

ADReal
INSFVMomentumDiffusion::computeQpResidual()
{
  auto dudn = gradUDotNormal();

  const auto face_mu = _mu(std::make_tuple(
      _face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains()));

  return -face_mu * dudn;
}
