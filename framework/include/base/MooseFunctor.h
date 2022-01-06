//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <tuple>

#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseError.h"
#include "Limiter.h"

#include "libmesh/elem.h"
#include "libmesh/quadrature.h"
#include "libmesh/remote_elem.h"
#include "libmesh/tensor_tools.h"

#include <unordered_map>
#include <functional>

class FaceInfo;

namespace Moose
{
struct ElemArg
{
  const libMesh::Elem * elem;
  bool correct_skewness;
  bool apply_gradient_to_skewness;
};

struct ElemFromFaceArg
{
  const libMesh::Elem * elem;
  const FaceInfo * fi;
  bool correct_skewness;
  bool apply_gradient_to_skewness;
  SubdomainID sub_id;

  ElemArg makeElem() const { return {elem, correct_skewness, apply_gradient_to_skewness}; }
};

struct FaceArg
{
  const FaceInfo * fi;
  Moose::FV::LimiterType limiter_type;
  bool elem_is_upwind;
  bool correct_skewness;
  bool apply_gradient_to_skewness;
  SubdomainID elem_sub_id;
  SubdomainID neighbor_sub_id;

  ElemArg makeElem() const { return {&fi->elem(), correct_skewness, apply_gradient_to_skewness}; }
  ElemArg makeNeighbor() const
  {
    return {fi->neighborPtr(), correct_skewness, apply_gradient_to_skewness};
  }
  ElemFromFaceArg elemFromFace() const
  {
    return {&fi->elem(), fi, correct_skewness, apply_gradient_to_skewness, elem_sub_id};
  }
  ElemFromFaceArg neighborFromFace() const
  {
    return {fi->neighborPtr(), fi, correct_skewness, apply_gradient_to_skewness, neighbor_sub_id};
  }
};

struct SingleSidedFaceArg
{
  const FaceInfo * fi;
  Moose::FV::LimiterType limiter_type;
  bool elem_is_upwind;
  bool correct_skewness;
  bool apply_gradient_to_skewness;
  SubdomainID sub_id;

  ElemArg makeElem() const { return {&fi->elem(), correct_skewness, apply_gradient_to_skewness}; }
  ElemArg makeNeighbor() const
  {
    return {fi->neighborPtr(), correct_skewness, apply_gradient_to_skewness};
  }
};

/**
 * Argument for requesting functor evaluation at a quadrature point location in an element. Data
 * in the argument:
 * - The element containing the quadrature point
 * - The quadrature point index, e.g. if there are \p n quadrature points, we are requesting the\n
 *   evaluation of the ith point
 * - The quadrature rule that can be used to initialize the functor on the given element
 */
using ElemQpArg = std::tuple<const libMesh::Elem *, unsigned int, const QBase *>;

/**
 * Argument for requesting functor evaluation at quadrature point locations on an element side.
 * Data in the argument:
 * - The element
 * - The element side on which the quadrature points are located
 * - The quadrature point index, e.g. if there are \p n quadrature points, we are requesting the\n
 *   evaluation of the ith point
 * - The quadrature rule that can be used to initialize the functor on the given element and side
 */
using ElemSideQpArg = std::tuple<const libMesh::Elem *, unsigned int, unsigned int, const QBase *>;

/**
 * Base class template for functor objects. This class template defines various \p operator()
 * overloads that allow a user to evaluate the functor at arbitrary geometric locations. This
 * template is meant to enable highly flexible on-the-fly variable and material property
 * evaluations
 */
template <typename T>
class FunctorImpl
{
public:
  using FunctorType = FunctorImpl<T>;
  using FunctorReturnType = T;
  using ValueType = T;
  using GradientType = typename libMesh::TensorTools::IncrementRank<T>::type;
  using DotType = ValueType;

  virtual ~FunctorImpl() = default;
  FunctorImpl() : _clearance_schedule({EXEC_ALWAYS}) {}
  FunctorImpl(const std::set<ExecFlagType> & clearance_schedule)
    : _clearance_schedule(clearance_schedule)
  {
  }

  ///@{
  /**
   * Same as their \p evaluate overloads with the same arguments but allows for caching
   * implementation. These are the methods a user will call in their code
   */
  ValueType operator()(const ElemArg & elem, unsigned int state = 0) const;
  ValueType operator()(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const;
  ValueType operator()(const FaceArg & face, unsigned int state = 0) const;
  ValueType operator()(const SingleSidedFaceArg & face, unsigned int state = 0) const;
  ValueType operator()(const ElemQpArg & qp, unsigned int state = 0) const;
  ValueType operator()(const ElemSideQpArg & qp, unsigned int state = 0) const;
  ///@}

  ///@{
  /**
   * Same as their \p evaluateGradient overloads with the same arguments but allows for caching
   * implementation. These are the methods a user will call in their code
   */
  GradientType gradient(const ElemArg & elem, unsigned int state = 0) const;
  GradientType gradient(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const;
  GradientType gradient(const FaceArg & face, unsigned int state = 0) const;
  GradientType gradient(const SingleSidedFaceArg & face, unsigned int state = 0) const;
  GradientType gradient(const ElemQpArg & qp, unsigned int state = 0) const;
  GradientType gradient(const ElemSideQpArg & qp, unsigned int state = 0) const;
  ///@}

  ///@{
  /**
   * Same as their \p evaluateDot overloads with the same arguments but allows for caching
   * implementation. These are the methods a user will call in their code
   */
  DotType dot(const ElemArg & elem, unsigned int state = 0) const;
  DotType dot(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const;
  DotType dot(const FaceArg & face, unsigned int state = 0) const;
  DotType dot(const SingleSidedFaceArg & face, unsigned int state = 0) const;
  DotType dot(const ElemQpArg & qp, unsigned int state = 0) const;
  DotType dot(const ElemSideQpArg & qp, unsigned int state = 0) const;
  ///@}

  virtual void residualSetup();
  virtual void jacobianSetup();
  virtual void timestepSetup();

  /**
   * Set how often to clear the functor evaluation cache
   */
  void setCacheClearanceSchedule(const std::set<ExecFlagType> & clearance_schedule);

  /**
   * Returns whether this is an extrapolated boundary face
   */
  virtual bool isExtrapolatedBoundaryFace(const FaceInfo &) const { mooseError("not implemented"); }

protected:
  /**
   * Evaluate the functor with a given element. Some example implementations of this method
   * could compute an element-average or evaluate at the element centroid
   */
  virtual ValueType evaluate(const ElemArg & elem, unsigned int state) const = 0;

  /**
   * @param elem_from_face See the \p ElemFromFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const ElemFromFaceArg & elem_from_face, unsigned int state) const = 0;

  /**
   * @param face See the \p FaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const FaceArg & face, unsigned int state) const = 0;

  /**
   * @param face See the \p SingleSidedFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const SingleSidedFaceArg & face, unsigned int state) const = 0;

  /**
   * @param qp See the \p ElemQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const ElemQpArg & qp, unsigned int state) const = 0;

  /**
   * @param side_qp See the \p ElemSideQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const ElemSideQpArg & side_qp, unsigned int state) const = 0;

  /**
   * Evaluate the functor gradient with a given element. Some example implementations of this
   * method could compute an element-average or evaluate at the element centroid
   */
  virtual GradientType evaluateGradient(const ElemArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param elem_from_face See the \p ElemFromFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const ElemFromFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p FaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const FaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p SingleSidedFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const SingleSidedFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param qp See the \p ElemQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const ElemQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param side_qp See the \p ElemSideQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const ElemSideQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * Evaluate the functor time derivative with a given element. Some example implementations of
   * this method could compute an element-average or evaluate at the element centroid
   */
  virtual DotType evaluateDot(const ElemArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param elem_from_face See the \p ElemFromFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const ElemFromFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p FaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const FaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p SingleSidedFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const SingleSidedFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param qp See the \p ElemQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const ElemQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param side_qp See the \p ElemSideQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const ElemSideQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

private:
  /**
   * clear cache data
   */
  void clearCacheData();

  /**
   * check a qp cache and if invalid then evaluate
   */
  template <typename SpaceArg, typename TimeArg>
  ValueType queryQpCache(unsigned int qp,
                         const QBase & qrule,
                         std::vector<std::pair<bool, T>> & qp_cache_data,
                         const SpaceArg & space,
                         const TimeArg & time) const;

  /// How often to clear the material property cache
  std::set<ExecFlagType> _clearance_schedule;

  // Data for traditional element-quadrature point property evaluations which are useful for
  // caching implementation

  /// Current key for qp map cache
  mutable dof_id_type _current_qp_map_key = DofObject::invalid_id;

  /// Current value for qp map cache
  mutable std::vector<std::pair<bool, T>> * _current_qp_map_value = nullptr;

  /// Cached element quadrature point functor property evaluations. The map key is the element
  /// id. The map values should have size corresponding to the number of quadrature points on the
  /// element. The vector elements are pairs. The first member of the pair indicates whether a
  /// cached value has been computed. The second member of the pair is the (cached) value. If the
  /// boolean is false, then the value cannot be trusted
  mutable std::unordered_map<dof_id_type, std::vector<std::pair<bool, T>>> _qp_to_value;

  // Data for traditional element-side-quadrature point property evaluations which are useful for
  // caching implementation

  /// Current key for side-qp map cache
  mutable dof_id_type _current_side_qp_map_key = DofObject::invalid_id;

  /// Current value for side-qp map cache
  mutable std::vector<std::vector<std::pair<bool, T>>> * _current_side_qp_map_value = nullptr;

  /// Cached element quadrature point functor property evaluations. The map key is the element
  /// id. The map values are a multi-dimensional vector (or vector of vectors) with the first index
  /// corresponding to the side and the second index corresponding to the quadrature point
  /// index. The elements returned after double indexing are pairs. The first member of the pair
  /// indicates whether a cached value has been computed. The second member of the pair is the
  /// (cached) value. If the boolean is false, then the value cannot be trusted
  mutable std::unordered_map<dof_id_type, std::vector<std::vector<std::pair<bool, T>>>>
      _side_qp_to_value;

  template <typename>
  friend class Functor;
};

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const ElemArg & elem, const unsigned int state) const
{
  return evaluate(elem, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const ElemFromFaceArg & elem_from_face, const unsigned int state) const
{
  return evaluate(elem_from_face, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const FaceArg & face, const unsigned int state) const
{
  return evaluate(face, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const SingleSidedFaceArg & face, const unsigned int state) const
{
  return evaluate(face, state);
}

template <typename T>
template <typename SpaceArg, typename TimeArg>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::queryQpCache(const unsigned int qp,
                             const QBase & qrule,
                             std::vector<std::pair<bool, T>> & qp_cache_data,
                             const SpaceArg & space,
                             const TimeArg & time) const
{
  // Check and see whether we even have sized for this quadrature point. If we haven't then we
  // must evaluate
  if (qp >= qp_cache_data.size())
  {
    qp_cache_data.resize(qrule.n_points(), std::make_pair(false, T()));
    auto & pr = qp_cache_data[qp];
    pr.second = evaluate(space, time);
    pr.first = true;
    return pr.second;
  }

  // We've already sized for this qp, so let's see whether we have a valid cache value
  auto & pr = qp_cache_data[qp];
  if (pr.first)
    return pr.second;

  // No valid cache value so evaluate
  pr.second = evaluate(space, time);
  pr.first = true;
  return pr.second;
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const ElemQpArg & elem_qp, const unsigned int state) const
{
  if (_clearance_schedule.count(EXEC_ALWAYS))
    return evaluate(elem_qp, state);

  const auto elem_id = std::get<0>(elem_qp)->id();
  if (elem_id != _current_qp_map_key)
  {
    _current_qp_map_key = elem_id;
    _current_qp_map_value = &_qp_to_value[elem_id];
  }
  auto & qp_data = *_current_qp_map_value;
  const auto qp = std::get<1>(elem_qp);
  const auto * const qrule = std::get<2>(elem_qp);
  mooseAssert(qrule, "qrule must be non-null");

  return queryQpCache(qp, *qrule, qp_data, elem_qp, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const ElemSideQpArg & elem_side_qp, const unsigned int state) const
{
  if (_clearance_schedule.count(EXEC_ALWAYS))
    return evaluate(elem_side_qp, state);

  const Elem * const elem = std::get<0>(elem_side_qp);
  mooseAssert(elem, "elem must be non-null");
  const auto elem_id = elem->id();
  if (elem_id != _current_side_qp_map_key)
  {
    _current_side_qp_map_key = elem_id;
    _current_side_qp_map_value = &_side_qp_to_value[elem_id];
  }
  auto & side_qp_data = *_current_side_qp_map_value;
  const auto side = std::get<1>(elem_side_qp);
  const auto qp = std::get<2>(elem_side_qp);
  const auto * const qrule = std::get<3>(elem_side_qp);
  mooseAssert(qrule, "qrule must be non-null");

  // Check and see whether we even have sized for this side
  if (side >= side_qp_data.size())
    side_qp_data.resize(elem->n_sides());

  // Ok we were sized enough for our side
  auto & qp_data = side_qp_data[side];
  return queryQpCache(qp, *qrule, qp_data, elem_side_qp, state);
}

template <typename T>
void
FunctorImpl<T>::setCacheClearanceSchedule(const std::set<ExecFlagType> & clearance_schedule)
{
  _clearance_schedule = clearance_schedule;
}

template <typename T>
void
FunctorImpl<T>::clearCacheData()
{
  for (auto & map_pr : _qp_to_value)
    for (auto & pr : map_pr.second)
      pr.first = false;

  for (auto & map_pr : _side_qp_to_value)
  {
    auto & side_vector = map_pr.second;
    for (auto & qp_vector : side_vector)
      for (auto & pr : qp_vector)
        pr.first = false;
  }

  _current_qp_map_key = DofObject::invalid_id;
  _current_qp_map_value = nullptr;
  _current_side_qp_map_key = DofObject::invalid_id;
  _current_side_qp_map_value = nullptr;
}

template <typename T>
void
FunctorImpl<T>::timestepSetup()
{
  if (_clearance_schedule.count(EXEC_TIMESTEP_BEGIN))
    clearCacheData();
}

template <typename T>
void
FunctorImpl<T>::residualSetup()
{
  if (_clearance_schedule.count(EXEC_LINEAR))
    clearCacheData();
}

template <typename T>
void
FunctorImpl<T>::jacobianSetup()
{
  if (_clearance_schedule.count(EXEC_NONLINEAR))
    clearCacheData();
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const ElemArg & elem, const unsigned int state) const
{
  return evaluateGradient(elem, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const ElemFromFaceArg & elem_from_face, const unsigned int state) const
{
  return evaluateGradient(elem_from_face, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const FaceArg & face, const unsigned int state) const
{
  return evaluateGradient(face, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const SingleSidedFaceArg & face, const unsigned int state) const
{
  return evaluateGradient(face, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const ElemQpArg & elem_qp, const unsigned int state) const
{
  return evaluateGradient(elem_qp, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const ElemSideQpArg & elem_side_qp, const unsigned int state) const
{
  return evaluateGradient(elem_side_qp, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const ElemArg & elem, const unsigned int state) const
{
  return evaluateDot(elem, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const ElemFromFaceArg & elem_from_face, const unsigned int state) const
{
  return evaluateDot(elem_from_face, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const FaceArg & face, const unsigned int state) const
{
  return evaluateDot(face, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const SingleSidedFaceArg & face, const unsigned int state) const
{
  return evaluateDot(face, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const ElemQpArg & elem_qp, const unsigned int state) const
{
  return evaluateDot(elem_qp, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const ElemSideQpArg & elem_side_qp, const unsigned int state) const
{
  return evaluateDot(elem_side_qp, state);
}

/**
 * A non-templated base class for functors that allow an owner object to hold
 * different class template instantiations of \p Functor in a single container
 */
class FunctorBase
{
public:
  FunctorBase() = default;
  virtual ~FunctorBase() = default;

  ///@{
  /**
   * Virtual methods meant to be used for handling functor evaluation cache clearance
   */
  virtual void timestepSetup() = 0;
  virtual void residualSetup() = 0;
  virtual void jacobianSetup() = 0;
  virtual bool wrapsNull() const = 0;
  virtual std::string returnType() const = 0;
  ///@}
};

template <typename T>
class NullFunctor;

/**
 * This is a wrapper that forwards calls to the implementation,
 * which can be switched out at any time without disturbing references to
 * FunctorImpl. Implementation motivated by https://stackoverflow.com/a/65455485/4493669
 */
template <typename T>
class Functor final : public FunctorImpl<T>, public FunctorBase
{
public:
  using typename Moose::FunctorImpl<T>::ValueType;
  using typename Moose::FunctorImpl<T>::GradientType;
  using typename Moose::FunctorImpl<T>::DotType;

  /**
   * Construct wrapper from wrapped object
   */
  Functor(const FunctorImpl<T> & wrapped) : FunctorBase(), _wrapped(&wrapped) {}
  Functor(FunctorImpl<T> &&) = delete;

  Functor(std::unique_ptr<FunctorImpl<T>> && wrapped)
    : FunctorBase(), _owned(std::move(wrapped)), _wrapped(_owned.get())
  {
  }

  /**
   * Assign our wrapped object to be something new and release our previously wrapped object
   */
  void assign(const FunctorImpl<T> & wrapped)
  {
    _owned.reset();
    _wrapped = &wrapped;
  }

  void assign(std::unique_ptr<FunctorImpl<T>> && wrapped)
  {
    _owned.reset();
    _owned = std::move(wrapped);
    _wrapped = _owned.get();
  }

  Functor(const Functor &) = delete;
  Functor(Functor &&) = delete;
  Functor & operator=(const Functor &) = delete;
  Functor & operator=(Functor &&) = delete;

  virtual ~Functor() = default;

  bool wrapsNull() const override { return wrapsType<NullFunctor<T>>(); }
  std::string returnType() const override { return libMesh::demangle(typeid(T).name()); }

  /**
   * Tests whether the wrapped object is of the requested type
   */
  template <typename T2>
  bool wrapsType() const
  {
    return dynamic_cast<const T2 *>(_wrapped);
  }

  void timestepSetup() override
  {
    if (_owned)
      _owned->timestepSetup();
  }
  void residualSetup() override
  {
    if (_owned)
      _owned->residualSetup();
  }
  void jacobianSetup() override
  {
    if (_owned)
      _owned->jacobianSetup();
  }

protected:
  ///@{
  /**
   * Forward calls to wrapped object
   */
  ValueType evaluate(const ElemArg & elem, unsigned int state = 0) const override
  {
    return _wrapped->operator()(elem, state);
  }
  ValueType evaluate(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const override
  {
    return _wrapped->operator()(elem_from_face, state);
  }
  ValueType evaluate(const FaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->operator()(face, state);
  }
  ValueType evaluate(const SingleSidedFaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->operator()(face, state);
  }
  ValueType evaluate(const ElemQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->operator()(qp, state);
  }
  ValueType evaluate(const ElemSideQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->operator()(qp, state);
  }

  GradientType evaluateGradient(const ElemArg & elem, unsigned int state = 0) const override
  {
    return _wrapped->gradient(elem, state);
  }
  GradientType evaluateGradient(const ElemFromFaceArg & elem_from_face,
                                unsigned int state = 0) const override
  {
    return _wrapped->gradient(elem_from_face, state);
  }
  GradientType evaluateGradient(const FaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->gradient(face, state);
  }
  GradientType evaluateGradient(const SingleSidedFaceArg & face,
                                unsigned int state = 0) const override
  {
    return _wrapped->gradient(face, state);
  }
  GradientType evaluateGradient(const ElemQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->gradient(qp, state);
  }
  GradientType evaluateGradient(const ElemSideQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->gradient(qp, state);
  }

  DotType evaluateDot(const ElemArg & elem, unsigned int state = 0) const override
  {
    return _wrapped->dot(elem, state);
  }
  DotType evaluateDot(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const override
  {
    return _wrapped->dot(elem_from_face, state);
  }
  DotType evaluateDot(const FaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->dot(face, state);
  }
  DotType evaluateDot(const SingleSidedFaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->dot(face, state);
  }
  DotType evaluateDot(const ElemQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->dot(qp, state);
  }
  DotType evaluateDot(const ElemSideQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->dot(qp, state);
  }
  ///@}

private:
  /// Our wrapped object
  std::unique_ptr<FunctorImpl<T>> _owned;
  const FunctorImpl<T> * _wrapped;

  friend class ::SubProblem;
};

/**
 * Class template for creating constants
 */
template <typename T>
class ConstantFunctor final : public FunctorImpl<T>
{
public:
  using typename FunctorImpl<T>::FunctorType;
  using typename FunctorImpl<T>::FunctorReturnType;
  using typename FunctorImpl<T>::ValueType;
  using typename FunctorImpl<T>::GradientType;
  using typename FunctorImpl<T>::DotType;

  ConstantFunctor(const ValueType & value) : _value(value) {}
  ConstantFunctor(ValueType && value) : _value(value) {}

private:
  ValueType evaluate(const ElemArg &, unsigned int) const override { return _value; }
  ValueType evaluate(const ElemFromFaceArg &, unsigned int) const override { return _value; }
  ValueType evaluate(const FaceArg &, unsigned int) const override { return _value; }
  ValueType evaluate(const SingleSidedFaceArg &, unsigned int) const override { return _value; }
  ValueType evaluate(const ElemQpArg &, unsigned int) const override { return _value; }
  ValueType evaluate(const ElemSideQpArg &, unsigned int) const override { return _value; }

  GradientType evaluateGradient(const ElemArg &, unsigned int) const override { return 0; }
  GradientType evaluateGradient(const ElemFromFaceArg &, unsigned int) const override { return 0; }
  GradientType evaluateGradient(const FaceArg &, unsigned int) const override { return 0; }
  GradientType evaluateGradient(const SingleSidedFaceArg &, unsigned int) const override
  {
    return 0;
  }
  GradientType evaluateGradient(const ElemQpArg &, unsigned int) const override { return 0; }
  GradientType evaluateGradient(const ElemSideQpArg &, unsigned int) const override { return 0; }

  DotType evaluateDot(const ElemArg &, unsigned int) const override { return 0; }
  DotType evaluateDot(const ElemFromFaceArg &, unsigned int) const override { return 0; }
  DotType evaluateDot(const FaceArg &, unsigned int) const override { return 0; }
  DotType evaluateDot(const SingleSidedFaceArg &, unsigned int) const override { return 0; }
  DotType evaluateDot(const ElemQpArg &, unsigned int) const override { return 0; }
  DotType evaluateDot(const ElemSideQpArg &, unsigned int) const override { return 0; }

private:
  ValueType _value;
};

/**
 * A functor that serves as a placeholder during the simulation setup phase if a functor consumer
 * requests a functor that has not yet been constructed.
 */
template <typename T>
class NullFunctor final : public FunctorImpl<T>
{
public:
  using typename FunctorImpl<T>::FunctorType;
  using typename FunctorImpl<T>::FunctorReturnType;
  using typename FunctorImpl<T>::ValueType;
  using typename FunctorImpl<T>::GradientType;
  using typename FunctorImpl<T>::DotType;

private:
  ValueType evaluate(const ElemArg &, unsigned int) const override
  {
    mooseError("we should never get here. If you have, contact a MOOSE developer and tell them "
               "they've written broken code");
  }
  ValueType evaluate(const ElemFromFaceArg &, unsigned int) const override
  {
    mooseError("we should never get here. If you have, contact a MOOSE developer and tell them "
               "they've written broken code");
  }
  ValueType evaluate(const FaceArg &, unsigned int) const override
  {
    mooseError("we should never get here. If you have, contact a MOOSE developer and tell them "
               "they've written broken code");
  }
  ValueType evaluate(const SingleSidedFaceArg &, unsigned int) const override
  {
    mooseError("we should never get here. If you have, contact a MOOSE developer and tell them "
               "they've written broken code");
  }
  ValueType evaluate(const ElemQpArg &, unsigned int) const override
  {
    mooseError("we should never get here. If you have, contact a MOOSE developer and tell them "
               "they've written broken code");
  }
  ValueType evaluate(const ElemSideQpArg &, unsigned int) const override
  {
    mooseError("we should never get here. If you have, contact a MOOSE developer and tell them "
               "they've written broken code");
  }
};
}
