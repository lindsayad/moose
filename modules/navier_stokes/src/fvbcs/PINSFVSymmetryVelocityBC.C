//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVSymmetryVelocityBC.h"

registerMooseObject("NavierStokesApp", PINSFVSymmetryVelocityBC);

InputParameters
PINSFVSymmetryVelocityBC::validParams()
{
  InputParameters params = INSFVSymmetryVelocityBC::validParams();
  params.addClassDescription(
      "Implements a free slip boundary condition using a penalty formulation.");
  params.addRequiredCoupledVar("porosity", "The porosity.");
  return params;
}

PINSFVSymmetryVelocityBC::PINSFVSymmetryVelocityBC(const InputParameters & params)
  : INSFVSymmetryVelocityBC(params), _eps(getFunctor<ADReal>("porosity"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("PINSFV is not supported by local AD indexing. In order to use PINSFV, please run "
             "the configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif
}

ADReal
PINSFVSymmetryVelocityBC::computeQpResidual()
{
  const Elem & elem = _face_info->elem();
  const Elem * const neighbor = _face_info->neighborPtr();

  auto ft = _face_info->faceType(_var.name());
  const bool var_on_elem_side = ft == FaceInfo::VarFaceNeighbors::ELEM;
  const Elem * const boundary_elem = var_on_elem_side ? &elem : neighbor;
  mooseAssert(boundary_elem, "the boundary elem should be non-null");

  const auto face_eps = _eps(std::make_tuple(
      _face_info, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id()));
  return INSFVSymmetryVelocityBC::computeQpResidual() / face_eps;
}

void
PINSFVSymmetryVelocityBC::gatherRCData(const FaceInfo & /*fi*/)
{
  // const Elem & elem = fi.elem();
  // const Elem * const neighbor = fi.neighborPtr();
  // const Point & normal = fi.normal();
  // Real coord;
  // coordTransformFactor(_subproblem, elem.subdomain_id(), fi.faceCentroid(), coord);
  // const auto surface_vector = normal * fi.faceArea() * coord;

  // auto ft = fi.faceType(_var.name());
  // const bool var_on_elem_side = ft == FaceInfo::VarFaceNeighbors::ELEM;
  // const Elem * const boundary_elem = var_on_elem_side ? &elem : neighbor;
  // const Point & boundary_elem_centroid =
  //     var_on_elem_side ? fi.elemCentroid() : fi.neighborCentroid();

  // mooseAssert(boundary_elem, "the boundary elem should be non-null");
  // const auto face_mu = _mu(std::make_tuple(
  //     &fi, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id()));
  // const auto face_eps = _eps(std::make_tuple(
  //     &fi, Moose::FV::LimiterType::CentralDifference, true, boundary_elem->subdomain_id()));

  // // Moukalled eqns. 15.154 - 15.156, adapted for porosity
  // const ADReal coeff = 2. * face_mu / face_eps * surface_vector.norm() /
  //                      std::abs((fi.faceCentroid() - boundary_elem_centroid) * normal) *
  //                      normal(_index) * normal(_index);
  // _rc_uo.addToA(boundary_elem, _index, coeff);
}
