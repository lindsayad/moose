//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVBodyForce.h"
#include "Function.h"

registerMooseObject("NavierStokesApp", INSFVBodyForce);

InputParameters
INSFVBodyForce::validParams()
{
  auto params = INSFVElementalKernel::validParams();
  params.addClassDescription("Body force that contributes to the Rhie-Chow interpolation");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function", "1", "A function that describes the body force");
  params.addParam<PostprocessorName>(
      "postprocessor", 1, "A postprocessor whose value is multiplied by the body force");
  params.declareControllable("value");
  return params;
}

INSFVBodyForce::INSFVBodyForce(const InputParameters & parameters)
  : INSFVElementalKernel(parameters),
    _scale(getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor"))
{
}

void
INSFVBodyForce::gatherRCData(const Elem & elem)
{
  _rc_uo.addToB(
      &elem, _index, -_scale * _postprocessor * _function.value(_t, elem.vertex_average()));
}
