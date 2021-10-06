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
  params += INSFVMomentumResidualObject::validParams();
  params.addRequiredParam<MaterialPropertyName>("mu", "The viscosity");
  // Yea this thing causes an insane stencil when doing Rhie-Chow
  params.set<unsigned short>("ghost_layers") = 3;
  return params;
}

INSFVMomentumDiffusion::INSFVMomentumDiffusion(const InputParameters & params)
  : FVFluxKernel(params),
    INSFVMomentumResidualObject(*this),
    _mu(getFunctor<ADReal>("mu")),
    _computing_rc_data(false)
{
}

void
INSFVMomentumDiffusion::gatherRCData(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;

  _face_info = &fi;
  _normal = fi.normal();
  _face_type = fi.faceType(_var.name());

  _computing_rc_data = true;
  const auto saved_do_derivatives = ADReal::do_derivatives;
  // We rely on derivative indexing
  ADReal::do_derivatives = true;
  const auto residual = fi.faceArea() * fi.faceCoord() * computeQpResidual();
  _computing_rc_data = false;
  ADReal::do_derivatives = saved_do_derivatives;

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    _rc_uo.addToA(&fi.elem(), _index, _ae);
  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    _rc_uo.addToA(fi.neighborPtr(), _index, _an);
}

ADReal
INSFVMomentumDiffusion::computeQpResidual()
{
  const auto face = std::make_tuple(
      _face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains());
  const auto dudn = _var.gradient(face) * _face_info->normal();
  const auto face_mu = _mu(face);

  if (_computing_rc_data)
  {
    if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
        _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    {
      const auto dof_number = _face_info->elem().dof_number(_sys.number(), _var.number(), 0);
      // A gradient is a linear combination of degrees of freedom so it's safe to straight-up index
      // into the derivatives vector at the dof we care about
      _ae = dudn.derivatives()[dof_number];
      _ae *= -face_mu;
    }
    if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
        _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    {
      const auto dof_number = _face_info->neighbor().dof_number(_sys.number(), _var.number(), 0);
      _an = dudn.derivatives()[dof_number];
      _an *= face_mu;
    }
  }

  return -face_mu * dudn;
}
