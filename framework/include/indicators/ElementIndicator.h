//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ELEMENTINDICATOR_H
#define ELEMENTINDICATOR_H

#include "Indicator.h"
#include "TransientInterface.h"
#include "PostprocessorInterface.h"
#include "Coupleable.h"
#include "ScalarCoupleable.h"
#include "MooseVariableInterface.h"
#include "MaterialPropertyInterface.h"

// Forward declarations
class ElementIndicator;
template <typename>
class MooseVariableFEImpl;
typedef MooseVariableFEImpl<Real> MooseVariable;
typedef MooseVariableFEImpl<VectorValue<Real>> VectorMooseVariable;

template <>
InputParameters validParams<ElementIndicator>();

class ElementIndicator : public Indicator,
                         public TransientInterface,
                         public PostprocessorInterface,
                         public Coupleable,
                         public ScalarCoupleable,
                         public MooseVariableInterface<Real>
{
public:
  ElementIndicator(const InputParameters & parameters);

protected:
  MooseVariableFE & _field_var;

  const Elem *& _current_elem;
  /// Volume of the current element
  const Real & _current_elem_volume;

  unsigned int _qp;
  const MooseArray<Point> & _q_point;
  QBase *& _qrule;
  const MooseArray<Real> & _JxW;
  const MooseArray<Real> & _coord;

  MooseVariable & _var;

  /// Holds the solution at current quadrature points
  const VariableValue & _u;

  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  /// Time derivative of u
  const VariableValue & _u_dot;

  /// Derivative of u_dot wrt u
  const VariableValue & _du_dot_du;

  /// Holds local indicator entries as their accumulated by this ElementIndicator
  DenseVector<Number> _local_indtr;
};

#endif /* ELEMENTINDICATOR_H */
