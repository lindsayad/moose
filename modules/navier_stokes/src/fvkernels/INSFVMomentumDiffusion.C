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
#include "NS.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", INSFVMomentumDiffusion);

InputParameters
INSFVMomentumDiffusion::validParams()
{
  auto params = FVFluxKernel::validParams();
  params += INSFVMomentumResidualObject::validParams();
  params.addRequiredParam<MaterialPropertyName>(NS::mu, "The viscosity");
  params.addClassDescription(
      "Implements the Laplace form of the viscous strees in the Navier-Stokes equation.");
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

INSFVMomentumDiffusion::INSFVMomentumDiffusion(const InputParameters & params)
  : FVFluxKernel(params),
    INSFVMomentumResidualObject(*this),
    _mu(getFunctor<ADReal>(NS::mu)),
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
  // Fill-in the coefficients _ae and _an (but without multiplication by A)
  computeQpResidual();
  _computing_rc_data = false;

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    _rc_uo.addToA(&fi.elem(), _index, _ae * (fi.faceArea() * fi.faceCoord()));
  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    _rc_uo.addToA(fi.neighborPtr(), _index, _an * (fi.faceArea() * fi.faceCoord()));
}

ADReal
INSFVMomentumDiffusion::computeQpResidual()
{
  const auto face = Moose::FV::makeCDFace(*_face_info, faceArgSubdomains());
  const auto dudn = gradUDotNormal();
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
