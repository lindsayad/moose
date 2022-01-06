//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVRhieChowInterpolator.h"
#include "CellCenteredMapFunctor.h"
#include <unordered_map>

class PINSFVRhieChowInterpolator : public INSFVRhieChowInterpolator
{
public:
  static InputParameters validParams();
  PINSFVRhieChowInterpolator(const InputParameters & params);

  void meshChanged() override;
  void residualSetup() override;

protected:
  const Moose::FunctorImpl<ADReal> & epsilon(THREAD_ID tid) const override;

  Moose::Functor<ADReal> & _eps;
  std::vector<const Moose::Functor<ADReal> *> _epss;
  const unsigned short _smoothing_layers;
  std::vector<const FaceInfo *> _geometric_fi;

  CellCenteredMapFunctor<ADReal, std::unordered_map<dof_id_type, ADReal>> _smoothed_eps;

private:
  void pinsfvSetup();
  bool _initial_setup_done = false;
};

inline const Moose::FunctorImpl<ADReal> &
PINSFVRhieChowInterpolator::epsilon(const THREAD_ID tid) const
{
  return *_epss[tid];
}
