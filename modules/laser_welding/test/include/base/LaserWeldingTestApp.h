//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef LASER_WELDINGTESTAPP_H
#define LASER_WELDINGTESTAPP_H

#include "MooseApp.h"

class LaserWeldingTestApp;

template <>
InputParameters validParams<LaserWeldingTestApp>();

class LaserWeldingTestApp : public MooseApp
{
public:
  LaserWeldingTestApp(InputParameters parameters);
  virtual ~LaserWeldingTestApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs = false);
};

#endif /* LASER_WELDINGTESTAPP_H */
