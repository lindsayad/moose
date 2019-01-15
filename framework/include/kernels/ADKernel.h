//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADKERNEL_H
#define ADKERNEL_H

#include "KernelBase.h"

#define usingTemplKernelMembers(type)                                                              \
  using ADKernelTempl<type, compute_stage>::_test;                                                 \
  using ADKernelTempl<type, compute_stage>::_qp;                                                   \
  using ADKernelTempl<type, compute_stage>::_i;                                                    \
  using ADKernelTempl<type, compute_stage>::_j;                                                    \
  using ADKernelTempl<type, compute_stage>::_u;                                                    \
  using ADKernelTempl<type, compute_stage>::_var;                                                  \
  using ADKernelTempl<type, compute_stage>::_grad_test;                                            \
  using ADKernelTempl<type, compute_stage>::_grad_u;                                               \
  using ADKernelTempl<type, compute_stage>::_JxW;                                                  \
  using ADKernelTempl<type, compute_stage>::_coord;                                                \
  using ADKernelTempl<type, compute_stage>::_q_point;                                              \
  using ADKernelTempl<type, compute_stage>::_local_re;                                             \
  using ADKernelTempl<type, compute_stage>::_local_ke;                                             \
  using ADKernelTempl<type, compute_stage>::_qrule;                                                \
  using ADKernelTempl<type, compute_stage>::_has_save_in;                                          \
  using ADKernelTempl<type, compute_stage>::_save_in;                                              \
  using ADKernelTempl<type, compute_stage>::_has_diag_save_in;                                     \
  using ADKernelTempl<type, compute_stage>::_diag_save_in;                                         \
  using ADKernelTempl<type, compute_stage>::_current_elem_volume;                                  \
  using ADKernelTempl<type, compute_stage>::_sys;                                                  \
  using ADKernelTempl<type, compute_stage>::_assembly;                                             \
  using ADKernelTempl<type, compute_stage>::_current_elem;                                         \
  using ADKernelTempl<type, compute_stage>::_t;                                                    \
  using ADKernelTempl<type, compute_stage>::_dt;                                                   \
  using ADKernelTempl<type, compute_stage>::coupled;                                               \
  using ADKernelTempl<type, compute_stage>::coupledComponents;                                     \
  using ADKernelTempl<type, compute_stage>::getBlockCoordSystem;                                   \
  using ADKernelTempl<type, compute_stage>::precalculateResidual;                                  \
  using ADKernelTempl<type, compute_stage>::prepareVectorTag;                                      \
  using ADKernelTempl<type, compute_stage>::prepareMatrixTag;                                      \
  using ADKernelTempl<type, compute_stage>::accumulateTaggedLocalResidual;                         \
  using ADKernelTempl<type, compute_stage>::accumulateTaggedLocalMatrix;                           \
  using ADKernelTempl<type, compute_stage>::variable;                                              \
  using ADKernelTempl<type, compute_stage>::paramError;                                            \
  using ADKernelTempl<type, compute_stage>::isParamValid;                                          \
  using ADKernelTempl<type, compute_stage>::isCoupled;                                             \
  using ADKernelTempl<type, compute_stage>::getFunction

#define usingKernelMembers usingTemplKernelMembers(Real)
#define usingVectorKernelMembers usingTemplKernelMembers(RealVectorValue)

// forward declarations
template <typename, ComputeStage>
class ADKernelTempl;
template <ComputeStage compute_stage>
using ADKernel = ADKernelTempl<Real, compute_stage>;
template <ComputeStage compute_stage>
using ADVectorKernel = ADKernelTempl<RealVectorValue, compute_stage>;

declareADValidParams(ADKernel);
declareADValidParams(ADVectorKernel);

template <typename T, ComputeStage compute_stage>
class ADKernelTempl : public KernelBase, public MooseVariableInterface<T>
{
public:
  ADKernelTempl(const InputParameters & parameters);

  virtual ~ADKernelTempl();

  // See KernelBase base for documentation of these overridden methods
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;
  virtual void computeADOffDiagJacobian() override;
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

  virtual MooseVariableFE<T> & variable() override { return _var; }

protected:
  /// Compute this Kernel's contribution to the residual at the current quadrature point
  virtual ADResidual computeQpResidual() = 0;

  /// This is a regular kernel so we cast to a regular MooseVariable
  MooseVariableFE<T> & _var;

  /// the current test function
  const ADTemplateVariableTestValue & _test;

  /// gradient of the test function
  const typename VariableTestGradientType<compute_stage, T>::type & _grad_test;

  /// Holds the solution at current quadrature points
  const ADTemplateVariableValue & _u;

  /// Holds the solution gradient at the current quadrature points
  const ADTemplateVariableGradient & _grad_u;
};

#endif /* ADKERNEL_H */
