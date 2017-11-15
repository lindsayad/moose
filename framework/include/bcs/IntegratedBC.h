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

#ifndef INTEGRATEDBC_H
#define INTEGRATEDBC_H

#include "IntegratedBCBase.h"

// Forward declarations
class IntegratedBC;

template <>
InputParameters validParams<IntegratedBC>();

/**
 * Base class for deriving any boundary condition of a integrated type
 */
class IntegratedBC : public IntegratedBCBase
{
public:
  IntegratedBC(const InputParameters & parameters);

  MooseVariable & variable() { return _var; }

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  /**
   * Computes d-ivar-residual / d-jvar...
   */
  virtual void computeJacobianBlock(unsigned int jvar) override;
  /**
   * Computes jacobian block with respect to a scalar variable
   * @param jvar The number of the scalar variable
   */
  void computeJacobianBlockScalar(unsigned int jvar) override;

protected:
  MooseVariable & _var;

  /// normals at quadrature points
  const MooseArray<Point> & _normals;

  // shape functions

  /// shape function values (in QPs)
  const VariablePhiValue & _phi;
  /// gradients of shape functions (in QPs)
  const VariablePhiGradient & _grad_phi;

  // test functions

  /// test function values (in QPs)
  const VariableTestValue & _test;
  /// gradients of test functions  (in QPs)
  const VariableTestGradient & _grad_test;

  // unknown

  /// the values of the unknown variable this BC is acting on
  const VariableValue & _u;
  /// the gradient of the unknown variable this BC is acting on
  const VariableGradient & _grad_u;
};

#endif /* INTEGRATEDBC_H */
