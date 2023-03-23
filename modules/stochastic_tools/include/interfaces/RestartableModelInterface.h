//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RestartableDataIO.h"
#include "StochasticToolsApp.h"

class RestartableModelInterface
{
public:
  static InputParameters validParams();

  RestartableModelInterface(const MooseObject * object, const std::string & meta_data_name);

  ///@{
  /**
   * Declare model data for loading from file as well as restart
   */
  // MOOSEDOCS_BEGIN
  template <typename T>
  T & declareModelData(const std::string & data_name);

  template <typename T>
  T & declareModelData(const std::string & data_name, const T & value);
  // MOOSEDOCS_END
  ///@}

  ///@{
  /**
   * Declare model data for loading from file as well as restart
   */
  template <typename T>
  const T & getModelData(const std::string & data_name) const;
  ///@}

  template <typename T>
  T & setModelData(const std::string & data_name);

  const std::string & modelMetaDataName() const { return _model_meta_data_name; }

  const std::string _model_meta_data_name;

private:
  /**
   * Internal function used by public declareModelData methods.
   */
  template <typename T>
  RestartableData<T> & declareModelDataHelper(const std::string & data_name);

  /**
   * Internal function used by public declareModelData methods.
   */
  template <typename T>
  RestartableData<T> & getModelDataHelper(const std::string & data_name) const;

  /// Reference to FEProblemBase instance
  const MooseObject * _object;
};

template <typename T>
T &
RestartableModelInterface::declareModelData(const std::string & data_name)
{
  RestartableData<T> & data_ref = declareModelDataHelper<T>(data_name);
  return data_ref.set();
}

template <typename T>
T &
RestartableModelInterface::declareModelData(const std::string & data_name, const T & value)
{
  RestartableData<T> & data_ref = declareModelDataHelper<T>(data_name);
  data_ref.set() = value;
  return data_ref.set();
}

template <typename T>
RestartableData<T> &
RestartableModelInterface::declareModelDataHelper(const std::string & data_name)
{
  auto data_ptr = std::make_unique<RestartableData<T>>(data_name, nullptr);
  RestartableDataValue & value = _object->getMooseApp().registerRestartableData(
      data_name, std::move(data_ptr), 0, false, _model_meta_data_name);
  RestartableData<T> & data_ref = static_cast<RestartableData<T> &>(value);
  return data_ref;
}

template <typename T>
const T &
RestartableModelInterface::getModelData(const std::string & data_name) const
{
  RestartableData<T> & data_ref = getModelDataHelper<T>(data_name);
  return data_ref.get();
}

template <typename T>
T &
RestartableModelInterface::setModelData(const std::string & data_name)
{
  RestartableData<T> & data_ref = getModelDataHelper<T>(data_name);
  return data_ref.set();
}

template <typename T>
RestartableData<T> &
RestartableModelInterface::getModelDataHelper(const std::string & data_name) const
{
  auto data_ptr = std::make_unique<RestartableData<T>>(data_name, nullptr);
  RestartableDataValue & value = _object->getMooseApp().registerRestartableData(
      data_name, std::move(data_ptr), 0, true, _model_meta_data_name);
  RestartableData<T> & data_ref = static_cast<RestartableData<T> &>(value);
  return data_ref;
}
