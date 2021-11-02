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

template <typename T>
class VectorComponentFunctor : public Moose::FunctorImpl<T>
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
  using VectorArg = typename libMesh::TensorTools::IncrementRank<T>::type;
  using VectorFunctor = Moose::FunctorImpl<VectorArg>;

  VectorComponentFunctor(const VectorFunctor & vector, const unsigned int component)
    : _vector(vector), _component(component)
  {
  }

  bool isExtrapolatedBoundaryFace(const FaceInfo & fi) const override
  {
    return _vector.isExtrapolatedBoundaryFace(fi);
  }

private:
  const VectorFunctor & _vector;
  const unsigned int _component;

  ValueType evaluate(const libMesh::Elem * const & elem,
                     const unsigned int state) const override final
  {
    return _vector(elem, state)(_component);
  }

  ValueType evaluate(const ElemFromFaceArg & elem_from_face,
                     const unsigned int state) const override final
  {
    return _vector(elem_from_face, state)(_component);
  }

  ValueType evaluate(const FaceArg & face, const unsigned int state) const override final
  {
    return _vector(face, state)(_component);
  }

  ValueType evaluate(const SingleSidedFaceArg & face, const unsigned int state) const override final
  {
    return _vector(face, state)(_component);
  }

  ValueType evaluate(const ElemQpArg & elem_qp, const unsigned int state) const override final
  {
    return _vector(elem_qp, state)(_component);
  }

  ValueType evaluate(const ElemSideQpArg & elem_side_qp,
                     const unsigned int state) const override final
  {
    return _vector(elem_side_qp, state)(_component);
  }
};
