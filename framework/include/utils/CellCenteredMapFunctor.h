//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseFunctor.h"
#include "GreenGaussGradient.h"
#include "MathFVUtils.h"
#include "libmesh/utility.h"
#include "libmesh/type_tensor.h"
#include "libmesh/compare_types.h"
#include "libmesh/threads.h"

template <typename T, typename T2, typename std::enable_if<ScalarTraits<T>::value, int>::type = 0>
inline TypeVector<typename CompareTypes<T, T2>::supertype>
outer_product(const T & a, const TypeVector<T2> & b)
{
  TypeVector<typename CompareTypes<T, T2>::supertype> ret;
  for (unsigned int i = 0; i < LIBMESH_DIM; i++)
    ret(i) = a * b(i);

  return ret;
}

template <typename T, typename Map>
class CellCenteredMapFunctor : public Moose::FunctorImpl<T>, public Map
{
public:
  using typename Moose::FunctorImpl<T>::ValueType;
  using typename Moose::FunctorImpl<T>::GradientType;
  using typename Moose::FunctorImpl<T>::DotType;

  CellCenteredMapFunctor(const MooseMesh & mesh,
                         const bool nonorthgonal_correction,
                         const bool map_filled = false)
    : _mesh(mesh), _nonorthgonal_correction(nonorthgonal_correction), _map_filled(map_filled)
  {
  }

  bool isExtrapolatedBoundaryFace(const FaceInfo & fi) const override { return !fi.neighborPtr(); }

  void mapFilled(const bool map_filled) { _map_filled = map_filled; }

private:
  const MooseMesh & _mesh;
  const bool _nonorthgonal_correction;
  bool _map_filled;
  mutable Threads::spin_mutex _map_mutex;

  ValueType evaluate(const Moose::ElemArg & elem_arg, unsigned int) const override final
  {
    const Elem * const elem = elem_arg.elem;
    if (_map_filled)
      return libmesh_map_find(*this, elem->id());

    Threads::spin_mutex::scoped_lock lock(_map_mutex);
    return (*const_cast<CellCenteredMapFunctor *>(this))[elem->id()];
  }

  ValueType evaluate(const Moose::FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *face.fi;
    mooseAssert(face.limiter_type == Moose::FV::LimiterType::CentralDifference,
                "this implementation currently only supports linear interpolations");
    if (fi.neighborPtr())
      return Moose::FV::linearInterpolation(*this, face);
    else
    {
      const auto elem_arg = face.makeElem();
      const auto elem_value = (*this)(elem_arg);
      // Two term expansion
      return elem_value + this->gradient(elem_arg) * (fi.faceCentroid() - fi.elemCentroid());
    }
  }

  using Moose::FunctorImpl<T>::evaluateGradient;

  GradientType evaluateGradient(const Moose::ElemArg & elem_arg, unsigned int) const override final
  {
    return Moose::FV::greenGaussGradient(elem_arg, *this, true, _mesh);
  }

  GradientType evaluateGradient(const Moose::FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *face.fi;
    const auto elem_arg = face.makeElem();
    const auto elem_gradient = this->gradient(elem_arg);
    if (fi.neighborPtr())
    {
      const auto neighbor_arg = face.makeNeighbor();
      const auto linear_interp_gradient =
          fi.gC() * elem_gradient + (1 - fi.gC()) * this->gradient(neighbor_arg);
      if (_nonorthgonal_correction)
        return linear_interp_gradient +
               outer_product(((*this)(neighbor_arg) - (*this)(elem_arg)) / fi.dCFMag() -
                                 linear_interp_gradient * fi.eCF(),
                             fi.eCF());
      else
        return linear_interp_gradient;
    }
    else
      // One term expansion
      return elem_gradient;
  }

  ValueType evaluate(const Moose::ElemFromFaceArg & elem_from_face, unsigned int) const override
  {
    const auto * elem = elem_from_face.elem;
    if (!elem)
      elem = &elem_from_face.fi->elem();
    return (*this)(Moose::ElemArg(
        {elem, elem_from_face.correct_skewness, elem_from_face.apply_gradient_to_skewness}));
  }

  ValueType evaluate(const Moose::SingleSidedFaceArg & ssf, unsigned int) const override
  {
    return (*this)(Moose::FV::makeCDFace(*ssf.fi));
  }

  ValueType evaluate(const Moose::ElemQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }

  ValueType evaluate(const Moose::ElemSideQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }
};
