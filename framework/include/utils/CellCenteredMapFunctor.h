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
  using typename Moose::FunctorImpl<T>::FaceArg;
  using typename Moose::FunctorImpl<T>::SingleSidedFaceArg;
  using typename Moose::FunctorImpl<T>::ElemFromFaceArg;
  using typename Moose::FunctorImpl<T>::ElemQpArg;
  using typename Moose::FunctorImpl<T>::ElemSideQpArg;
  using typename Moose::FunctorImpl<T>::ValueType;
  using typename Moose::FunctorImpl<T>::GradientType;
  using typename Moose::FunctorImpl<T>::DotType;

  CellCenteredMapFunctor(const MooseMesh & mesh,
                         const bool nonorthgonal_correction,
                         const bool correct_skewness = false)
    : _mesh(mesh),
      _nonorthgonal_correction(nonorthgonal_correction),
      _correct_skewness(correct_skewness)
  {
  }

  bool isExtrapolatedBoundaryFace(const FaceInfo & fi) const override { return !fi.neighborPtr(); }

  using Moose::FunctorImpl<T>::operator();
  ValueType operator()(const FaceInfo & fi) const { return (*this)(makeFace(fi)); }

  using Moose::FunctorImpl<T>::gradient;
  GradientType gradient(const FaceInfo & fi) const { return gradient(makeFace(fi)); }

private:
  FaceArg makeFace(const FaceInfo & fi) const { return Moose::FV::makeCDFace(fi); }

  const MooseMesh & _mesh;
  const bool _nonorthgonal_correction;
  const bool _correct_skewness;

  ValueType evaluate(const libMesh::Elem * const & elem, unsigned int) const override final
  {
    return (*const_cast<CellCenteredMapFunctor *>(this))[elem->id()];
  }

  ValueType evaluate(const FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *std::get<0>(face);
    const auto elem_value = (*this)(&fi.elem());
    if (fi.neighborPtr())
      return fi.gC() * elem_value + (1 - fi.gC()) * (*this)(fi.neighborPtr());
    else
      // Two term expansion
      return elem_value + this->gradient(&fi.elem()) * (fi.faceCentroid() - fi.elemCentroid());
  }

  using Moose::FunctorImpl<T>::evaluateGradient;

  GradientType evaluateGradient(const libMesh::Elem * const & elem,
                                unsigned int) const override final
  {
    return Moose::FV::greenGaussGradient(elem, *this, true, _correct_skewness, _mesh);
  }

  GradientType evaluateGradient(const FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *std::get<0>(face);
    const auto elem_gradient = this->gradient(&fi.elem());
    if (fi.neighborPtr())
    {
      const auto linear_interp_gradient =
          fi.gC() * elem_gradient + (1 - fi.gC()) * this->gradient(fi.neighborPtr());
      if (_nonorthgonal_correction)
        return linear_interp_gradient +
               outer_product(((*this)(fi.neighborPtr()) - (*this)(&fi.elem())) / fi.dCFMag() -
                                 linear_interp_gradient * fi.eCF(),
                             fi.eCF());
      else
        return linear_interp_gradient;
    }
    else
      // One term expansion
      return elem_gradient;
  }

  ValueType evaluate(const ElemFromFaceArg & elem_from_face, unsigned int) const override
  {
    const auto * elem_arg = std::get<0>(elem_from_face);
    if (!elem_arg)
      elem_arg = &std::get<1>(elem_from_face)->elem();
    return (*this)(elem_arg);
  }

  ValueType evaluate(const SingleSidedFaceArg & ssf, unsigned int) const override
  {
    return (*this)(makeFace(*std::get<0>(ssf)));
  }

  ValueType evaluate(const ElemQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }

  ValueType evaluate(const ElemSideQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }
};
