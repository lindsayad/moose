//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InputParameters.h"
#include "INSFVRhieChowInterpolator.h"

class MooseObject;
class FaceInfo;

class INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  template <typename T>
  INSFVMomentumResidualObject(T & obj);
  virtual void gatherRCData(const Elem & elem) = 0;
  virtual void gatherRCData(const FaceInfo & fi) = 0;

  virtual ~INSFVMomentumResidualObject() = default;

protected:
  INSFVRhieChowInterpolator & _rc_uo;

  /// index x|y|z
  const unsigned int _index;
};

template <typename T>
INSFVMomentumResidualObject::INSFVMomentumResidualObject(T & obj)
  : _rc_uo(const_cast<INSFVRhieChowInterpolator &>(
        obj.template getUserObject<INSFVRhieChowInterpolator>("rhie_chow_user_object"))),
    _index(obj.template getParam<MooseEnum>("momentum_component"))
{
}
