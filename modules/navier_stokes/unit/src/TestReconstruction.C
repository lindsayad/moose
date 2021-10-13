#include "gtest/gtest.h"

#include "Registry.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "NavierStokesApp.h"
#include "AppFactory.h"
#include "Factory.h"
#include "InputParameters.h"
#include "MeshGeneratorMesh.h"
#include "MooseError.h"
#include "CastUniquePointer.h"
#include "GeneratedMeshGenerator.h"

#include <memory>

TEST(TestReconstruction, theTest)
{
  const char * argv[2] = {"foo", "\0"};
  auto app = AppFactory::createAppShared("NavierStokesApp", 1, (char **)argv);
  auto * factory = &app->getFactory();
  std::string mesh_type = "MeshGeneratorMesh";

  std::vector<unsigned int> num_elem = {2, 4, 8, 16, 32};

  for (const auto nx : num_elem)
  {
    std::shared_ptr<MeshGeneratorMesh> mesh;
    {
      InputParameters params = factory->getValidParams(mesh_type);
      mesh = factory->create<MeshGeneratorMesh>(mesh_type, "moose_mesh", params);
    }

    app->actionWarehouse().mesh() = mesh;
    std::unique_ptr<MeshBase> lm_mesh;

    {
      InputParameters params = factory->getValidParams("GeneratedMeshGenerator");
      params.set<unsigned int>("nx") = nx;
      params.set<unsigned int>("ny") = nx;
      params.set<MooseEnum>("dim") = "2";
      auto mesh_gen =
          factory->create<GeneratedMeshGenerator>("GeneratedMeshGenerator", "mesh_gen", params);
      lm_mesh = mesh_gen->generate();
    }

    mesh->setMeshBase(std::move(lm_mesh));

    const auto & all_fi = mesh->allFaceInfo();
  }
}
