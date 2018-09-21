//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTERESIDUALTHREAD_H
#define COMPUTERESIDUALTHREAD_H

#include "ThreadedElementLoop.h"
#include "MooseObjectTagWarehouse.h"

#include "libmesh/elem_range.h"

// Forward declarations
class FEProblemBase;
class NonlinearSystemBase;
class IntegratedBCBase;
class DGKernel;
class InterfaceKernel;
class TimeKernel;
class KernelBase;
class Kernel;

class ComputeResidualThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  ComputeResidualThread(FEProblemBase & fe_problem, const std::set<TagID> & tags);

  // Splitting Constructor
  ComputeResidualThread(ComputeResidualThread & x, Threads::split split);

  virtual ~ComputeResidualThread();

  virtual void subdomainChanged() override;
  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInterface(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void postElement(const Elem * /*elem*/) override;
  virtual void post() override;

  void join(const ComputeResidualThread & /*y*/);

protected:
  NonlinearSystemBase & _nl;
  const std::set<TagID> & _tags;
  unsigned int _num_cached;

  /// Reference to BC storage structures
  const MooseObjectWarehouse<IntegratedBCBase> & _integrated_bcs;

  /// Reference to DGKernel storage structure
  const MooseObjectWarehouse<DGKernel> & _dg_kernels;

  /// Reference to interface kernel storage structure
  const MooseObjectWarehouse<InterfaceKernel> & _interface_kernels;

  ///@{
  /// Reference to Kernel storage structures
  MooseObjectTagWarehouse<KernelBase> & _kernels;

  MooseObjectWarehouse<KernelBase> * _tag_kernels;
  ///@}
};

#endif // COMPUTERESIDUALTHREAD_H
