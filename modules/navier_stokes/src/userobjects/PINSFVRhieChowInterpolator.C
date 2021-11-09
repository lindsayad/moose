//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVRhieChowInterpolator.h"
#include "Reconstructions.h"

registerMooseObject("NavierStokesApp", PINSFVRhieChowInterpolator);

InputParameters
PINSFVRhieChowInterpolator::validParams()
{
  auto params = INSFVRhieChowInterpolator::validParams();
  params.addParam<MooseFunctorName>(NS::porosity, "The porosity");
  params.addParam<unsigned int>(
      "reconstructions", 0, "The number of reconstructions to perform on the porosity");
  return params;
}

PINSFVRhieChowInterpolator::PINSFVRhieChowInterpolator(const InputParameters & params)
  : INSFVRhieChowInterpolator(params),
    _eps(isParamValid(NS::porosity)
             ? &const_cast<Moose::Functor<ADReal> &>(getFunctor<ADReal>(NS::porosity))
             : nullptr),
    _rec(getParam<unsigned int>("reconstructions")),
    _done(false)
{
  if (_rec && !_eps)
    paramError("reconstructions",
               "If a non-zero number of 'reconstructions' is applied, then the parameter '",
               NS::porosity,
               "' must be supplied.");
}

void
PINSFVRhieChowInterpolator::residualSetup()
{
  INSFVRhieChowInterpolator::residualSetup();

  if (!_rec || _done)
    return;

  CellCenteredMapFunctor<ADReal, std::unordered_map<dof_id_type, ADReal>> reconstructed_eps(
      _moose_mesh, true);
  ADReal::do_derivatives = true;
  Moose::FV::reconstruct(reconstructed_eps, *_eps, _rec, false, false, _moose_mesh);
  ADReal::do_derivatives = false;

  (*_eps) = std::move(reconstructed_eps);

  _done = true;
}
