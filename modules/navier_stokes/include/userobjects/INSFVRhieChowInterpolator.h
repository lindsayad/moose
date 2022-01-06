//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "TaggingInterface.h"
#include "BlockRestrictable.h"
#include "ADReal.h"
#include "MooseTypes.h"
#include "CellCenteredMapFunctor.h"
#include "VectorComponentFunctor.h"
#include "libmesh/vector_value.h"
#include "libmesh/id_types.h"
#include "libmesh/stored_range.h"
#include <unordered_map>
#include <set>
#include <unordered_set>

class MooseMesh;
class INSFVVelocityVariable;
class INSFVPressureVariable;
namespace libMesh
{
class Elem;
class MeshBase;
}

class INSFVRhieChowInterpolator : public GeneralUserObject,
                                  public TaggingInterface,
                                  public BlockRestrictable
{
public:
  static InputParameters validParams();
  INSFVRhieChowInterpolator(const InputParameters & params);

  void addToA(const libMesh::Elem * elem, unsigned int component, const ADReal & value);
  void addToB(const libMesh::Elem * elem, unsigned int component, const ADReal & value);
  VectorValue<ADReal>
  getVelocity(Moose::FV::InterpMethod m, const FaceInfo & fi, THREAD_ID tid) const;

  void initialSetup() override;
  void meshChanged() override;

  void initialize() override final;
  void execute() override final;
  void finalize() override final;

protected:
  /**
   * @return whether this face is geometrically relevant to us
   */
  bool isFaceGeometricallyRelevant(const FaceInfo & fi) const;

  /**
   * A virtual method that allows us to only implement getVelocity once for free and porous flows
   */
  virtual const Moose::FunctorImpl<ADReal> & epsilon(THREAD_ID tid) const;

  MooseMesh & _moose_mesh;
  const libMesh::MeshBase & _mesh;
  const unsigned int _dim;
  /// A functor for computing the (non-RC corrected) velocity
  std::vector<std::unique_ptr<PiecewiseByBlockLambdaFunctor<ADRealVectorValue>>> _vel;

  INSFVPressureVariable * const _p;
  INSFVVelocityVariable * const _u;
  INSFVVelocityVariable * const _v;
  INSFVVelocityVariable * const _w;

  std::vector<MooseVariableFVReal *> _ps;
  std::vector<MooseVariableFVReal *> _us;
  std::vector<MooseVariableFVReal *> _vs;
  std::vector<MooseVariableFVReal *> _ws;

  std::unique_ptr<ConstElemRange> _elem_range;

  std::vector<const FaceInfo *> _evaluable_fi;

  std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>> _a;
  CellCenteredMapFunctor<libMesh::VectorValue<ADReal>,
                         std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>>>
      _b;
  // Here the suffix on _b refers to the number of bar operations we've performed
  CellCenteredMapFunctor<libMesh::VectorValue<ADReal>,
                         std::unordered_map<dof_id_type, libMesh::VectorValue<ADReal>>>
      _b2;

private:
  void insfvSetup();
  void finalizeAData();
  void computeFirstAndSecondOverBars();
  void computeThirdOverBar();
  void applyBData();
  void finalizeBData();

  std::vector<unsigned int> _var_numbers;
  std::unordered_set<const Elem *> _elements_to_push_pull;

  SystemBase & _sys;
  const VectorValue<ADReal> _example;
  const bool _standard_body_forces;

  VectorComponentFunctor<ADReal> _bx;
  VectorComponentFunctor<ADReal> _by;
  VectorComponentFunctor<ADReal> _b2x;
  VectorComponentFunctor<ADReal> _b2y;

  /// The subdomain ids this object operates on
  const std::set<SubdomainID> _sub_ids;

  /// Mutex that prevents multiple threads from saving into the 'a' coefficients at the same time
  Threads::spin_mutex _a_mutex;

  /// Mutex that prevents multiple threads from saving into the 'b' coefficients  at the same time
  Threads::spin_mutex _b_mutex;

  /// A unity functor used in the epsilon virtual method
  const Moose::ConstantFunctor<ADReal> _unity_functor{1};
};

inline const Moose::FunctorImpl<ADReal> & INSFVRhieChowInterpolator::epsilon(THREAD_ID) const
{
  return _unity_functor;
}

inline void
INSFVRhieChowInterpolator::addToA(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  Threads::spin_mutex::scoped_lock lock(_a_mutex);

  if (elem->processor_id() != this->processor_id())
    _elements_to_push_pull.insert(elem);

  _a[elem->id()](component) += value;
}

inline void
INSFVRhieChowInterpolator::addToB(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  mooseAssert(elem->processor_id() == this->processor_id(), "Sources should be local");

  Threads::spin_mutex::scoped_lock lock(_b_mutex);
  // We have our users write their RC data imagining that they've moved all terms to the LHS, but
  // the balance in Moukalled assumes that the body forces are on the RHS with positive sign, e.g.
  // 0 = -\nabla p + \mathbf{B}, so we must apply a minus sign here
  _b[elem->id()](component) -= value;
}
