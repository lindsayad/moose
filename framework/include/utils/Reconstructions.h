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
#include "libmesh/elem.h"

#include <unordered_map>
#include <utility>

namespace Moose
{
namespace FV
{
template <typename T, typename Map>
void
tanoReconstruction(CellCenteredMapFunctor<T, Map> & output_functor,
                   const CellCenteredMapFunctor<T, Map> & input_functor,
                   const unsigned int num_reconstructions,
                   const MooseMesh & mesh)
{
  if (!num_reconstructions)
    return;

  std::unordered_map<dof_id_type, std::pair<T, Real>> elem_to_num_denom;

  const auto & all_fi = mesh.allFaceInfo();
  for (const auto & fi : all_fi)
  {
    auto & elem_pr = elem_to_num_denom[fi.elem().id()];
    auto face_value = input_functor(fi);
    if (fi.neighborPtr() && fi.neighborPtr() != libMesh::remote_elem)
    {
      auto & neighbor_pr = elem_to_num_denom[fi.neighbor().id()];
      neighbor_pr.first += face_value;
      neighbor_pr.second += 1.;
    }
    elem_pr.first += std::move(face_value);
    elem_pr.second += 1.;
  }

  for (const auto & pr : elem_to_num_denom)
  {
    const auto & data_pr = pr.second;
    output_functor[pr.first] = data_pr.first / data_pr.second;
  }

  tanoReconstruction(output_functor, output_functor, num_reconstructions - 1, mesh);
}
}
}
