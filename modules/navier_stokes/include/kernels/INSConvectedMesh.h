//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INSCONVECTEDMESH_H
#define INSCONVECTEDMESH_H

#include "Kernel.h"

// Forward Declaration
class INSConvectedMesh;

template <>
InputParameters validParams<INSConvectedMesh>();

/**
 * This calculates the time derivative for a coupled variable
 **/
class INSConvectedMesh : public Kernel
{
public:
  INSConvectedMesh(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _disp_x_dot;
  const VariableValue & _d_disp_x_dot;
  const unsigned int _disp_x_id;
  const VariableValue & _disp_y_dot;
  const VariableValue & _d_disp_y_dot;
  const unsigned int _disp_y_id;
  const VariableValue & _disp_z_dot;
  const VariableValue & _d_disp_z_dot;
  const unsigned int _disp_z_id;

  const MaterialProperty<Real> & _rho;
};

#endif // INSCONVECTEDMESH}_H
