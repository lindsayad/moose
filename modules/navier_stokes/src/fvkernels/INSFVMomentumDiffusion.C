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
  : FVFluxKernel(params), INSFVMomentumResidualObject(*this), _mu(getFunctor<ADReal>("mu"))
{
}

void
INSFVMomentumDiffusion::gatherRCData(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;

  _face_info = &fi;
  _normal = fi.normal();

  const auto residual = fi.faceArea() * fi.faceCoord() * computeQpResidual();

  const auto ft = fi.faceType(_var.name());
  if (ft == FaceInfo::VarFaceNeighbors::ELEM || ft == FaceInfo::VarFaceNeighbors::BOTH)
  {
    const auto dof_number = fi.elem().dof_number(_sys.number(), _var.number(), 0);
    const auto var_value = _var(&fi.elem());
    const auto a = residual / var_value;
    // residual contribution of this kernel to the elem element
    _rc_uo.addToA(&fi.elem(), _index, residual / _var(&fi.elem()));
  }
  if (ft == FaceInfo::VarFaceNeighbors::NEIGHBOR || ft == FaceInfo::VarFaceNeighbors::BOTH)
  {
    const auto dof_number = fi.neighbor().dof_number(_sys.number(), _var.number(), 0);
    const auto var_value = _var(&fi.neighbor());
    const auto a = -residual / var_value;
    // residual contribution of this kernel to the neighbor element
    _rc_uo.addToA(fi.neighborPtr(), _index, -residual / _var(fi.neighborPtr()));
  }
}

ADReal
INSFVMomentumDiffusion::computeQpResidual()
{
  const auto face = std::make_tuple(
      _face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains());
  const auto dudn = _var.gradient(face) * _face_info->normal();
  const auto face_mu = _mu(face);

  return -face_mu * dudn;
}
