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
#include "FaceInfo.h"
#include "CellCenteredMapFunctor.h"
#include "MathFVUtils.h"
#include "libmesh/elem.h"

#include <unordered_map>
#include <utility>

namespace Moose
{
namespace FV
{
template <typename T, typename Map, typename Consumer>
void
reconstruct(CellCenteredMapFunctor<T, Map> & output_functor,
            const FunctorImpl<T> & input_functor,
            const unsigned int num_reconstructions,
            const bool two_term_expansion,
            const bool weight_with_sf,
            const std::vector<const FaceInfo *> & faces,
            const Consumer & consumer)
{
  if (!num_reconstructions)
    return;

  std::unordered_map<dof_id_type, std::pair<T, Real>> elem_to_num_denom;

  for (const auto * const face : faces)
  {
    mooseAssert(face, "This must be non-null");
    const Real weight = weight_with_sf ? face->faceArea() * face->faceCoord() : 1;
    const auto sub_pair = faceArgSubdomains(consumer, *face);
    const auto face_arg = makeCDFace(*face, sub_pair);
    auto face_value = input_functor(face_arg);
    std::pair<T, Real> * neighbor_pr = nullptr;
    if (face->neighborPtr() && face->neighborPtr() != libMesh::remote_elem)
    {
      neighbor_pr = &elem_to_num_denom[face->neighbor().id()];
      neighbor_pr->first += face_value * weight;
      neighbor_pr->second += weight;
    }
    auto & elem_pr = elem_to_num_denom[face->elem().id()];
    elem_pr.first += std::move(face_value) * weight;
    elem_pr.second += weight;

    if (two_term_expansion)
    {
      auto face_gradient = input_functor.gradient(face_arg);
      if (face->neighborPtr() && face->neighborPtr() != libMesh::remote_elem)
        neighbor_pr->first +=
            face_gradient * (face->neighborCentroid() - face->faceCentroid()) * weight;
      elem_pr.first +=
          std::move(face_gradient) * (face->elemCentroid() - face->faceCentroid()) * weight;
    }
  }

  for (const auto & pr : elem_to_num_denom)
  {
    const auto & data_pr = pr.second;
    output_functor[pr.first] = data_pr.first / data_pr.second;
  }

  reconstruct(output_functor,
              output_functor,
              num_reconstructions - 1,
              two_term_expansion,
              weight_with_sf,
              faces,
              consumer);
}
}
}
