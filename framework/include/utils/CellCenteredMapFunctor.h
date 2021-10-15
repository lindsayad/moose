//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseFunctor.h"
#include "GreenGaussGradient.h"
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
class CellCenteredMapFunctor : public Moose::Functor<T>, public Map
{
public:
  using typename Moose::Functor<T>::FaceArg;
  using typename Moose::Functor<T>::SingleSidedFaceArg;
  using typename Moose::Functor<T>::ElemFromFaceArg;
  using typename Moose::Functor<T>::ElemQpArg;
  using typename Moose::Functor<T>::ElemSideQpArg;
  using typename Moose::Functor<T>::ValueType;
  using typename Moose::Functor<T>::GradientType;
  using typename Moose::Functor<T>::DotType;

  CellCenteredMapFunctor(const MooseMesh & mesh, const bool correct_face)
    : _mesh(mesh), _correct_face(correct_face)
  {
  }

  bool isExtrapolatedBoundaryFace(const FaceInfo & fi) const override { return !fi.neighborPtr(); }

private:
  const MooseMesh & _mesh;
  const bool _correct_face;

  ValueType evaluate(const libMesh::Elem * const & elem, unsigned int) const override final
  {
    return libmesh_map_find(*this, elem->id());
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

  using Moose::Functor<T>::evaluateGradient;

  GradientType evaluateGradient(const libMesh::Elem * const & elem,
                                unsigned int) const override final
  {
    return Moose::FV::greenGaussGradient(elem, *this, true, _mesh, Moose::COORD_XYZ);
  }

  GradientType evaluateGradient(const FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *std::get<0>(face);
    const auto elem_gradient = this->gradient(&fi.elem());
    if (fi.neighborPtr())
    {
      const auto linear_interp_gradient =
          fi.gC() * elem_gradient + (1 - fi.gC()) * this->gradient(fi.neighborPtr());
      if (_correct_face)
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

  ValueType evaluate(const ElemFromFaceArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }

  ValueType evaluate(const SingleSidedFaceArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }

  ValueType evaluate(const ElemQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }

  ValueType evaluate(const ElemSideQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }

  ValueType evaluate(const std::tuple<Moose::ElementType, unsigned int, SubdomainID> &,
                     unsigned int) const override
  {
    mooseError("not implemented");
  }
};
