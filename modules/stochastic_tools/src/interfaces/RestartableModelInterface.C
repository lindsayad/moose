//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RestartableModelInterface.h"

InputParameters
RestartableModelInterface::validParams()
{
  return emptyInputParameters();
}

RestartableModelInterface::RestartableModelInterface(const MooseObject * object,
                                                     const std::string & meta_data_name)
  : _model_meta_data_name(meta_data_name), _object(object)
{
  _object->getMooseApp().registerRestartableDataMapName(_model_meta_data_name, _object->name());
}
