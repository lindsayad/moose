//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseError.h"
#include "FunctorInterface.h"
#include "Moose.h"
#include "SetupInterface.h"

#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"

#include <unordered_map>
#include <functional>

/**
 * A material property that is evaluated on-the-fly via calls to various overloads of \p operator()
 */
template <typename T>
class FunctorMaterialProperty : public FunctorInterface<T>

{
public:
  FunctorMaterialProperty(const std::string & name) : _name(name) {}

  /**
   * Set the functor that will be used in calls to \p evaluate overloads
   * @param mesh The mesh that the functor is defined on
   * @param block_ids The block/subdomain IDs that the user-provided functor is valid for
   * @param my_lammy The functor that defines this how this object is evaluated
   */
  template <typename PolymorphicLambda>
  void setFunctor(const MooseMesh & mesh,
                  const std::set<SubdomainID> & block_ids,
                  PolymorphicLambda my_lammy);

private:
  using typename FunctorInterface<T>::FaceArg;
  using typename FunctorInterface<T>::ElemAndFaceArg;
  using typename FunctorInterface<T>::QpArg;
  using typename FunctorInterface<T>::FunctorType;
  using typename FunctorInterface<T>::FunctorReturnType;

  using ElemFn = std::function<T(const Elem * const &, const unsigned int &)>;
  using ElemAndFaceFn = std::function<T(const ElemAndFaceArg &, const unsigned int &)>;
  using FaceFn = std::function<T(const FaceArg &, const unsigned int &)>;
  using QpFn = std::function<T(const QpArg &, const unsigned int &)>;
  using TQpFn = std::function<T(const std::tuple<Moose::ElementType, unsigned int, SubdomainID> &,
                                const unsigned int &)>;

  T evaluate(const Elem * const & elem, unsigned int state) const override final;
  T evaluate(const ElemAndFaceArg & elem_and_face, unsigned int state) const override final;
  T evaluate(const FaceArg & face, unsigned int state) const override final;
  T evaluate(const QpArg & qp, unsigned int state) const override final;
  T evaluate(const std::tuple<Moose::ElementType, unsigned int, SubdomainID> & tqp,
             unsigned int state) const override final;

  /// Functors that return element average values (or cell centroid values or whatever the
  /// implementer wants to return for a given element argument)
  std::unordered_map<SubdomainID, ElemFn> _elem_functor;

  /// Functors that return the value on the requested element that will perform any necessary
  /// ghosting operations if this object is not technically defined on the requested subdomain
  std::unordered_map<SubdomainID, ElemAndFaceFn> _elem_and_face_functor;

  /// Functors that return potentially limited interpolations at faces
  std::unordered_map<SubdomainID, FaceFn> _face_functor;

  /// Functors that will index elemental data at a provided quadrature point index
  std::unordered_map<SubdomainID, QpFn> _qp_functor;

  /// Functors that will index elemental, neighbor, or lower-dimensional data at a provided
  /// quadrature point index
  std::unordered_map<SubdomainID, TQpFn> _tqp_functor;

  /// The name of this object
  std::string _name;
};

template <typename T>
template <typename PolymorphicLambda>
void
FunctorMaterialProperty<T>::setFunctor(const MooseMesh & mesh,
                                       const std::set<SubdomainID> & block_ids,
                                       PolymorphicLambda my_lammy)
{
  auto add_lammy = [this, my_lammy](const SubdomainID block_id) {
    auto pr = _elem_functor.emplace(block_id, my_lammy);
    if (!pr.second)
      mooseError("No insertion for the functor material property '",
                 _name,
                 "' for block id ",
                 block_id,
                 ". Another material must already declare this property on that block.");
    _elem_and_face_functor.emplace(block_id, my_lammy);
    _face_functor.emplace(block_id, my_lammy);
    _qp_functor.emplace(block_id, my_lammy);
    _tqp_functor.emplace(block_id, my_lammy);
  };

  for (const auto block_id : block_ids)
  {
    if (block_id == Moose::ANY_BLOCK_ID)
    {
      const auto & inner_block_ids = mesh.meshSubdomains();
      for (const auto inner_block_id : inner_block_ids)
        add_lammy(inner_block_id);
    }
    else
      add_lammy(block_id);
  }
}

template <typename T>
T
FunctorMaterialProperty<T>::evaluate(const Elem * const & elem, unsigned int state) const
{
  mooseAssert(elem && elem != libMesh::remote_elem,
              "The element must be non-null and non-remote in functor material properties");
  auto it = _elem_functor.find(elem->subdomain_id());
  mooseAssert(it != _elem_functor.end(), "The provided subdomain ID doesn't exist in the map!");
  return it->second(elem, state);
}

template <typename T>
T
FunctorMaterialProperty<T>::evaluate(const ElemAndFaceArg & elem_and_face, unsigned int state) const
{
  mooseAssert((std::get<0>(elem_and_face) && std::get<0>(elem_and_face) != libMesh::remote_elem) ||
                  std::get<1>(elem_and_face),
              "The element must be non-null and non-remote or the face must be non-null in functor "
              "material properties");
  auto it = _elem_and_face_functor.find(std::get<2>(elem_and_face));
  mooseAssert(it != _elem_and_face_functor.end(),
              "The provided subdomain ID doesn't exist in the map!");
  return it->second(elem_and_face, state);
}

template <typename T>
T
FunctorMaterialProperty<T>::evaluate(const FaceArg & face, unsigned int state) const
{
  mooseAssert(std::get<0>(face), "FaceInfo must be non-null");
  auto it = _face_functor.find(std::get<3>(face));
  mooseAssert(it != _face_functor.end(), "The provided subdomain ID doesn't exist in the map!");
  return it->second(face, state);
}

template <typename T>
T
FunctorMaterialProperty<T>::evaluate(const QpArg & elem_and_qp, unsigned int state) const
{
  auto it = _qp_functor.find(elem_and_qp.first->subdomain_id());
  mooseAssert(it != _qp_functor.end(),
              "The provided element has a subdomain ID that doesn't exist in the map!");
  return it->second(elem_and_qp, state);
}

template <typename T>
T
FunctorMaterialProperty<T>::evaluate(
    const std::tuple<Moose::ElementType, unsigned int, SubdomainID> & tqp, unsigned int state) const
{
  auto it = _tqp_functor.find(std::get<2>(tqp));
  mooseAssert(it != _tqp_functor.end(), "The provided subdomain ID doesn't exist in the map!");
  return it->second(tqp, state);
}