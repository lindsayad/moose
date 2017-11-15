/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifndef VECTORFEWAVE_H
#define VECTORFEWAVE_H

#include "VectorKernel.h"
#include "MaterialProperty.h"

// Forward Declaration
class VectorFEWave;

template <>
InputParameters validParams<VectorFEWave>();

class VectorFEWave : public VectorKernel
{
public:
  VectorFEWave(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Function & _x_ffn;
  Function & _y_ffn;
  Function & _z_ffn;
};

#endif // VECTORFEWAVE_H
