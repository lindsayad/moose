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
#include "NS.h"

registerMooseObject("NavierStokesApp", PINSFVRhieChowInterpolator);

InputParameters
PINSFVRhieChowInterpolator::validParams()
{
  auto params = INSFVRhieChowInterpolator::validParams();
  params.addParam<MooseFunctorName>(NS::porosity, "The porosity");
  params.addParam<unsigned short>(
      "reconstructions", 0, "The number of reconstructions to perform on the porosity");
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC,
      [](const InputParameters & obj_params, InputParameters & rm_params) {
        rm_params.set<unsigned short>("layers") = obj_params.get<unsigned short>("reconstructions");
        rm_params.set<bool>("use_displaced_mesh") = obj_params.get<bool>("use_displaced_mesh");
      });
  return params;
}

PINSFVRhieChowInterpolator::PINSFVRhieChowInterpolator(const InputParameters & params)
  : INSFVRhieChowInterpolator(params),
    _eps(isParamValid(NS::porosity)
             ? &const_cast<Moose::Functor<ADReal> &>(getFunctor<ADReal>(NS::porosity))
             : nullptr),
    _rec(getParam<unsigned short>("reconstructions"))
{
  if (_rec)
  {
    if (_eps)
    {
      if (_eps->wrapsType<MooseVariableBase>())
        paramError(
            NS::porosity,
            "If we are reconstructing porosity, then the input porosity to this user object cannot "
            "be a Moose variable. There are issues with reconstructing Moose variables: 1) initial "
            "conditions are run after use object initial setup 2) reconstructing from a variable "
            "requires ghosting the solution vectors 3) it's difficult to restrict the face "
            "informations we evaluate reconstructions on such that we never query an algebraically "
            "remote element due to things like two-term extrapolated boundary faces which trigger "
            "gradient evaluations which trigger neighbor element evaluation");
    }
    else
      paramError("reconstructions",
                 "If a non-zero number of 'reconstructions' is applied, then the parameter '",
                 NS::porosity,
                 "' must be supplied.");
  }
}

void
PINSFVRhieChowInterpolator::interpolatorSetup()
{
  INSFVRhieChowInterpolator::interpolatorSetup();

  if (!_rec)
    return;

  const auto & all_fi = _moose_mesh.allFaceInfo();
  _geometric_fi.reserve(all_fi.size());

  for (const auto & fi : all_fi)
    if (isFaceGeometricallyRelevant(fi))
      _geometric_fi.push_back(&fi);

  _geometric_fi.shrink_to_fit();

  CellCenteredMapFunctor<ADReal, std::unordered_map<dof_id_type, ADReal>> reconstructed_eps(
      _moose_mesh, true);
  ADReal::do_derivatives = true;
  Moose::FV::reconstruct(reconstructed_eps, *_eps, _rec, false, false, _geometric_fi, *this);
  ADReal::do_derivatives = false;

  (*_eps) = std::move(reconstructed_eps);
}
