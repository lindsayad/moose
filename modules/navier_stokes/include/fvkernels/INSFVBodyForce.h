//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVBodyForce.h"
#include "INSFVResidualObject.h"

/**
 * Body force that contributes to the Rhie-Chow interpolation
 */
class INSFVBodyForce : public FVBodyForce, public INSFVResidualObject
{
public:
  INSFVBodyForce(const InputParameters & params);
  static InputParameters validParams();

  // requires RC implementation
  void gatherRCData(const Elem &) override {}

  void gatherRCData(const FaceInfo &) override {}
};
