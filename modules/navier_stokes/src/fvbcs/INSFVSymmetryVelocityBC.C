//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVSymmetryVelocityBC.h"

registerMooseObject("NavierStokesApp", INSFVSymmetryVelocityBC);

InputParameters
INSFVSymmetryVelocityBC::validParams()
{
  InputParameters params = INSFVSymmetryBC::validParams();
  params += INSFVMomentumResidualObject::validParams();
  params.addClassDescription(
      "Implements a free slip boundary condition using a penalty formulation.");
  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", 0, "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", 0, "The velocity in the z direction.");
  params.addRequiredParam<MaterialPropertyName>("mu", "The viscosity");
  return params;
}

INSFVSymmetryVelocityBC::INSFVSymmetryVelocityBC(const InputParameters & params)
  : INSFVSymmetryBC(params),
    INSFVMomentumResidualObject(*this),
    _u_functor(getFunctor<ADReal>("u")),
    _v_functor(getFunctor<ADReal>("v")),
    _w_functor(getFunctor<ADReal>("w")),
    _mu(getFunctor<ADReal>("mu")),
    _dim(_subproblem.mesh().dimension())
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("INSFV is not supported by local AD indexing. In order to use INSFV, please run the "
             "configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif
}

ADReal
INSFVSymmetryVelocityBC::computeQpResidual()
{
  const bool use_elem = _face_info->faceType(_var.name()) == FaceInfo::VarFaceNeighbors::ELEM;
  const auto normal = use_elem ? _face_info->normal() : Point(-_face_info->normal());
  const Point & cell_centroid =
      use_elem ? _face_info->elemCentroid() : _face_info->neighborCentroid();
  _u_eval = use_elem ? _u_functor(&_face_info->elem()) : _u_functor(_face_info->neighborPtr());
  _v_eval = use_elem ? _v_functor(&_face_info->elem()) : _v_functor(_face_info->neighborPtr());
  _w_eval = use_elem ? _w_functor(&_face_info->elem()) : _w_functor(_face_info->neighborPtr());

  // Evaluate viscosity on the face
  const auto mu_b = use_elem ? _mu(std::make_tuple(_face_info,
                                                   Moose::FV::LimiterType::CentralDifference,
                                                   true,
                                                   _face_info->elem().subdomain_id()))
                             : _mu(std::make_tuple(_face_info,
                                                   Moose::FV::LimiterType::CentralDifference,
                                                   true,
                                                   _face_info->neighborPtr()->subdomain_id()));

  const auto d_perpendicular = std::abs((_face_info->faceCentroid() - cell_centroid) * normal);

  // See Moukalled 15.150. Recall that we multiply by the area in the base class, so S_b ->
  // normal.norm() here

  ADReal v_dot_n = _u_eval * normal(0);
  if (_dim > 1)
    v_dot_n += _v_eval * normal(1);
  if (_dim > 2)
    v_dot_n += _w_eval + normal(2);

  return 2. * mu_b * normal.norm() / d_perpendicular * v_dot_n * normal(_index);
}

void
INSFVSymmetryVelocityBC::gatherRCData(const FaceInfo & fi)
{
  _face_info = &fi;
  const auto residual = fi.faceArea() * fi.faceCoord() * computeQpResidual();
  const auto & var = [this]() {
    switch (_index)
    {
      case 0:
        return _u_eval;
      case 1:
        return _v_eval;
      case 2:
        return _w_eval;
      default:
        mooseError("Unrecognized index value");
    }
  }();

  auto * const boundary_elem = (fi.faceType(_var.name()) == FaceInfo::VarFaceNeighbors::ELEM)
                                   ? &fi.elem()
                                   : fi.neighborPtr();

  _rc_uo.addToA(boundary_elem, _index, residual / var);
}
