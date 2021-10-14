//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "Registry.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseUnitApp.h"
#include "AppFactory.h"
#include "Factory.h"
#include "InputParameters.h"
#include "MeshGeneratorMesh.h"
#include "MooseError.h"
#include "CastUniquePointer.h"
#include "GeneratedMeshGenerator.h"
#include "MooseFunctor.h"
#include "PolynomialFit.h"
#include "FaceInfo.h"
#include "MooseTypes.h"
#include "CellCenteredMapFunctor.h"

#include "libmesh/elem.h"
#include "libmesh/tensor_value.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/utility.h"

#include <memory>
#include <vector>
#include <memory>

TEST(TestReconstruction, theTest)
{
  const char * argv[2] = {"foo", "\0"};

  std::vector<unsigned int> num_elem = {64, 128, 256};
  std::vector<Real> errors;
  std::vector<Real> weller_errors;
  std::vector<Real> h(num_elem.size());
  for (const auto i : index_range(num_elem))
    h[i] = 1. / num_elem[i];

  for (const auto i : index_range(num_elem))
  {
    const auto nx = num_elem[i];
    auto app = AppFactory::createAppShared("MooseUnitApp", 1, (char **)argv);
    auto * factory = &app->getFactory();
    std::string mesh_type = "MeshGeneratorMesh";

    std::shared_ptr<MeshGeneratorMesh> mesh;
    {
      InputParameters params = factory->getValidParams(mesh_type);
      mesh = factory->create<MeshGeneratorMesh>(mesh_type, "moose_mesh", params);
    }

    app->actionWarehouse().mesh() = mesh;

    {
      std::unique_ptr<MeshBase> lm_mesh;
      InputParameters params = factory->getValidParams("GeneratedMeshGenerator");
      params.set<unsigned int>("nx") = nx;
      params.set<unsigned int>("ny") = nx;
      params.set<MooseEnum>("dim") = "2";
      auto mesh_gen =
          factory->create<GeneratedMeshGenerator>("GeneratedMeshGenerator", "mesh_gen", params);
      lm_mesh = mesh_gen->generate();
      mesh->setMeshBase(std::move(lm_mesh));
    }

    auto & lm_mesh = mesh->getMesh();

    std::unordered_map<dof_id_type, RealVectorValue> analytic_map;
    for (auto * const elem : lm_mesh.active_element_ptr_range())
    {
      const auto centroid = elem->vertex_average();
      analytic_map[elem->id()] = RealVectorValue(-std::sin(centroid(0)) * std::cos(centroid(1)),
                                                 std::cos(centroid(0)) * std::sin(centroid(1)));
    }

    CellCenteredMapFunctor<RealVectorValue, decltype(analytic_map)> u(*mesh,
                                                                      std::move(analytic_map));

    const auto & all_fi = mesh->allFaceInfo();
    mesh->applyCoordSysToFaceCoords(coord_type, rz_radial_coord);

    std::unordered_map<dof_id_type, RealVectorValue> up;
    std::unordered_map<dof_id_type, RealVectorValue> up_weller;
    std::unordered_map<dof_id_type, Real> sf_sfhat_sum;

    for (const auto & fi : all_fi)
    {
      const auto face =
          std::make_tuple(&fi,
                          Moose::FV::LimiterType::CentralDifference,
                          true,
                          std::make_pair(fi.elem().subdomain_id(),
                                         fi.neighborPtr() ? fi.neighbor().subdomain_id()
                                                          : Moose::INVALID_BLOCK_ID));
      const RealVectorValue uf(u(face));
      const RealTensorValue grad_uf(u.gradient(face));

      const Point surface_vector = fi.normal() * fi.faceArea() * fi.faceCoord();
      const auto elem_interpolant = uf + grad_uf * (fi.elemCentroid() - fi.faceCentroid());
      const auto sf_sfhat = fi.normal() * surface_vector;

      up[fi.elem().id()] += elem_interpolant * sf_sfhat;
      up_weller[fi.elem().id()] += uf * sf_sfhat;
      sf_sfhat_sum[fi.elem().id()] += sf_sfhat;

      if (fi.neighborPtr())
      {
        const auto neighbor_interpolant =
            uf + grad_uf * (fi.neighborCentroid() - fi.faceCentroid());
        up[fi.neighbor().id()] += neighbor_interpolant * sf_sfhat;
        up_weller[fi.neighbor().id()] += uf * sf_sfhat;
        sf_sfhat_sum[fi.neighbor().id()] += sf_sfhat;
      }
    }

    Real error = 0;
    Real weller_error = 0;
    const auto current_h = h[i];
    for (auto & pr : up)
    {
      const auto elem_id = pr.first;
      auto & up_current = pr.second;
      const auto sf_sfhat = libmesh_map_find(sf_sfhat_sum, elem_id);
      up_current /= sf_sfhat;
      auto & up_weller_current = libmesh_map_find(up_weller, elem_id);
      up_weller_current /= sf_sfhat;
      auto * elem = lm_mesh.elem_ptr(elem_id);
      const RealVectorValue analytic(u(elem));
      const auto diff = analytic - up_current;
      error += diff * diff * current_h * current_h;
      const auto weller_diff = analytic - up_weller_current;
      weller_error += weller_diff * weller_diff * current_h * current_h;
    }
    error = std::sqrt(error);
    weller_error = std::sqrt(weller_error);
    errors.push_back(error);
    weller_errors.push_back(weller_error);
  }

  for (auto error : errors)
    EXPECT_LT(error, std::numeric_limits<Real>::epsilon());

  std::for_each(h.begin(), h.end(), [](Real & h_elem) { h_elem = std::log(h_elem); });
  std::for_each(weller_errors.begin(), weller_errors.end(), [](Real & weller_error) {
    weller_error = std::log(weller_error);
  });
  PolynomialFit weller_fit(h, weller_errors, 1);
  weller_fit.generate();

  const auto & coeffs = weller_fit.getCoefficients();
  EXPECT_NEAR(coeffs[1], 2., .05);
}
