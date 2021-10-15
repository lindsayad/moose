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

void
testReconstruction(const Moose::CoordinateSystemType coord_type,
                   const unsigned int rz_radial_coord = libMesh::invalid_uint)
{
  const char * argv[2] = {"foo", "\0"};

  std::vector<unsigned int> num_elem = {64, 128, 256};
  std::vector<Real> errors;
  std::vector<Real> linear_errors;
  std::vector<Real> weller_errors;
  std::vector<Real> moukalled_errors;
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

    CellCenteredMapFunctor<RealVectorValue, std::unordered_map<dof_id_type, RealVectorValue>> u(
        *mesh, true);
    CellCenteredMapFunctor<RealVectorValue, std::unordered_map<dof_id_type, RealVectorValue>>
        u_linear(*mesh, false);
    for (auto * const elem : lm_mesh.active_element_ptr_range())
    {
      const auto centroid = elem->vertex_average();
      const auto value = RealVectorValue(-std::sin(centroid(0)) * std::cos(centroid(1)),
                                         std::cos(centroid(0)) * std::sin(centroid(1)));
      u[elem->id()] = value;
      u_linear[elem->id()] = value;
    }

    const auto & all_fi = mesh->allFaceInfo();
    mesh->applyCoordSysToFaceCoords(coord_type, rz_radial_coord);

    std::unordered_map<dof_id_type, RealVectorValue> up;
    std::unordered_map<dof_id_type, RealVectorValue> up_linear;
    std::unordered_map<dof_id_type, RealVectorValue> up_weller;
    std::unordered_map<dof_id_type, RealVectorValue> up_moukalled;
    std::unordered_map<dof_id_type, Real> sf_sfhat_sum;

    for (const auto & fi : all_fi)
    {
      const Point surface_vector = fi.normal() * fi.faceArea() * fi.faceCoord();
      const auto sf_sfhat = fi.normal() * surface_vector;
      sf_sfhat_sum[fi.elem().id()] += sf_sfhat;
      if (fi.neighborPtr())
        sf_sfhat_sum[fi.neighbor().id()] += sf_sfhat;

      auto interpolate = [&fi, &sf_sfhat](auto & functor, auto & container, const bool aguerre) {
        const RealVectorValue uf(functor(fi));
        const RealTensorValue grad_uf(functor.gradient(fi));
        auto elem_interpolant = uf;
        if (aguerre)
          elem_interpolant += grad_uf * (fi.elemCentroid() - fi.faceCentroid());
        container[fi.elem().id()] += elem_interpolant * sf_sfhat;
        if (fi.neighborPtr())
        {
          auto neighbor_interpolant = uf;
          if (aguerre)
            neighbor_interpolant += grad_uf * (fi.neighborCentroid() - fi.faceCentroid());
          container[fi.neighbor().id()] += neighbor_interpolant * sf_sfhat;
        }
      };

      auto moukalled_interpolate = [&fi](auto & functor, auto & container) {
        const RealVectorValue uf(functor(fi));
        const Point surface_vector = fi.normal() * fi.faceArea();
        auto product = (uf * fi.dCF()) * surface_vector;

        container[fi.elem().id()] += product * fi.gC() / fi.elemVolume();
        if (fi.neighborPtr())
          container[fi.neighbor().id()] +=
              std::move(product) * (1. - fi.gC()) / fi.neighborVolume();
      };

      interpolate(u, up, true);
      interpolate(u, up_weller, false);
      interpolate(u_linear, up_linear, true);
      moukalled_interpolate(u, up_moukalled);
    }

    Real error = 0;
    Real weller_error = 0;
    Real linear_error = 0;
    Real moukalled_error = 0;
    const auto current_h = h[i];
    for (auto * const elem : lm_mesh.active_element_ptr_range())
    {
      const auto elem_id = elem->id();
      const auto sf_sfhat = libmesh_map_find(sf_sfhat_sum, elem_id);
      const RealVectorValue analytic(u(elem));

      auto compute_elem_error = [elem_id, current_h, &sf_sfhat, &analytic](
                                    auto & container, auto & error, const bool apply_sfhat) {
        auto & current = libmesh_map_find(container, elem_id);
        if (apply_sfhat)
          current /= sf_sfhat;
        const auto diff = analytic - current;
        error += diff * diff * current_h * current_h;
      };

      compute_elem_error(up, error, true);
      compute_elem_error(up_weller, weller_error, true);
      compute_elem_error(up_linear, linear_error, true);
      compute_elem_error(up_moukalled, moukalled_error, false);
    }
    error = std::sqrt(error);
    weller_error = std::sqrt(weller_error);
    linear_error = std::sqrt(linear_error);
    moukalled_error = std::sqrt(moukalled_error);
    errors.push_back(error);
    weller_errors.push_back(weller_error);
    linear_errors.push_back(linear_error);
    moukalled_errors.push_back(moukalled_error);
  }

  for (auto error : errors)
    EXPECT_LT(error, std::numeric_limits<Real>::epsilon());

  std::for_each(h.begin(), h.end(), [](Real & h_elem) { h_elem = std::log(h_elem); });

  auto expect_errors = [&h, coord_type](auto & errors_arg, Real expected_error) {
    std::for_each(
        errors_arg.begin(), errors_arg.end(), [](Real & error) { error = std::log(error); });
    PolynomialFit fit(h, errors_arg, 1);
    fit.generate();

    const auto & coeffs = fit.getCoefficients();
    EXPECT_NEAR(coeffs[1], expected_error, .05);
  };

  expect_errors(weller_errors, coord_type == Moose::COORD_RZ ? 1.5 : 2);
  expect_errors(linear_errors, 2.5);
  expect_errors(moukalled_errors, 2);
}

TEST(TestReconstruction, Cartesian) { testReconstruction(Moose::COORD_XYZ); }

TEST(TestReconstruction, Cylindrical) { testReconstruction(Moose::COORD_RZ, 0); }
