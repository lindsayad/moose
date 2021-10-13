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

#include <memory>

class XYFunctor : public Functor<Real>
{
public:
  using typename Functor<Real>::FaceArg;
  using typename Functor<Real>::SingleSidedFaceArg;
  using typename Functor<Real>::ElemFromFaceArg;
  using typename Functor<Real>::ElemQpArg;
  using typename Functor<Real>::ElemSideQpArg;

  XYFunctor(const MooseMesh & mesh) : _mesh(mesh) {}

private:
  const MooseMesh & _mesh;

  virtual Real value(const Point & point) const = 0;

  Real evaluate(const libMesh::Elem * const & elem, unsigned int) const override final
  {
    return value(elem->vertex_average());
  }

  Real evaluate(const ElemFromFaceArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }

  Real evaluate(const FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *std::get<0>(face);
    const auto elem_value = (*this)(&fi.elem());
    if (fi.neighborPtr())
      return fi.gC() * elem_value + (1 - fi.gC()) * (*this)(fi.neighborPtr());
    else
      // Two term expansion
      return elem_value + this->gradient(&fi.elem()) * (fi->faceCentroid() - fi->elemCentroid())
  }

  Real evaluate(const SingleSidedFaceArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  Real evaluate(const ElemQpArg &, unsigned int) const override final { mooseError("Not needed"); }
  Real evaluate(const ElemSideQpArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  Real evaluate(const std::tuple<Moose::ElementType, unsigned int, SubdomainID> &,
                unsigned int) const override final
  {
    mooseError("Not needed");
  }

  VectorValue<Real> evaluateGradient(const libMesh::Elem * const & elem,
                                     unsigned int) const override final
  {
    return Moose::FV::greenGaussGradient(elem, *this, true, _mesh, Moose::XYZ);
  }

  VectorValue<Real> evaluateGradient(const ElemFromFaceArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }

  VectorValue<Real> evaluateGradient(const FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *std::get<0>(face);
    const auto elem_gradient = this->gradient(&fi.elem());
    if (fi.neighborPtr())
      return fi.gC() * elem_gradient + (1 - fi.gC()) * this->gradient(fi.neighborPtr());
    else
      // One term expansion
      return elem_gradient;
  }

  VectorValue<Real> evaluateGradient(const SingleSidedFaceArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  VectorValue<Real> evaluateGradient(const ElemQpArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  VectorValue<Real> evaluateGradient(const ElemSideQpArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  VectorValue<Real>
  evaluateGradient(const std::tuple<Moose::ElementType, unsigned int, SubdomainID> &,
                   unsigned int) const override final
  {
    mooseError("Not needed");
  }

  Real evaluateDot(const libMesh::Elem * const &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  Real evaluateDot(const ElemFromFaceArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  Real evaluateDot(const FaceArg &, unsigned int) const override final { mooseError("Not needed"); }
  Real evaluateDot(const SingleSidedFaceArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  Real evaluateDot(const ElemQpArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  Real evaluateDot(const ElemSideQpArg &, unsigned int) const override final
  {
    mooseError("Not needed");
  }
  Real evaluateDot(const std::tuple<Moose::ElementType, unsigned int, SubdomainID> &,
                   unsigned int) const override final
  {
    mooseError("Not needed");
  }
};

class UXFunctor : public XYFunctor
{
public:
  UXFunctor(const MooseMesh & mesh) : XYFunctor(mesh) {}

private:
  Real value(const Point & point) const override final
  {
    return -std::sin(point(0)) * std::cos(point(1));
  }
};

class UYFunctor : public XYFunctor
{
public:
  UXFunctor(const MooseMesh & mesh) : XYFunctor(mesh) {}

private:
  Real value(const Point & point) const override final
  {
    return std::cos(point(0)) * std::sin(point(1));
  }
};

TEST(TestReconstruction, theTest)
{
  const char * argv[2] = {"foo", "\0"};
  auto app = AppFactory::createAppShared("MooseUnitApp", 1, (char **)argv);
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

    UXFunctor ux(*mesh);
    UYFunctor uy(*mesh);
    const auto & all_fi = mesh->allFaceInfo();

    for (const auto & fi : all_fi)
    {
      const auto face =
          std::make_tuple(&fi,
                          Moose::FV::LimiterType::CentralDifference,
                          true,
                          std::make_pair(fi.elem().subdomain_id(),
                                         fi.neighborPtr() ? fi.neighbor().subdomain_id()
                                                          : Moose::INVALID_BLOCK_ID));
      RealVectorValue uf(ux(face), uy(face));
      RealTensorValue grad_uf;
      RealVectorValue grad_uxf = ux.gradient(face);
      RealVectorValue grad_uyf = uy.gradient(face);
      for (const auto i : make_range(unsigned(2)))
      {
        grad_uf(0, i) = grad_uxf(i);
        grad_uf(1, i) = grad_uyf(i);
      }

      const Point surface_vector = fi.normal() * fi.faceArea() * fi.faceCoord();
      const auto Ff = uf * surface_vector;

      const auto sf_sfhat = fi.normal() * surface_vector;
      const auto weighted_flux = it->second * sf_sfhat;

      // Now should we compute _b2 for this element?
      if (dof_map.is_evaluable(fi.elem(), _var_numbers[0]))
      {
        _b2[fi.elem().id()] += weighted_flux;
        sf_sfhat_sum[fi.elem().id()] += sf_sfhat;
      }

      // Or for the neighbor?
      if (fi.neighborPtr() && dof_map.is_evaluable(fi.neighbor(), _var_numbers[0]))
      {
        _b2[fi.neighbor().id()] += weighted_flux;
        sf_sfhat_sum[fi.neighbor().id()] += sf_sfhat;
      }
    }

    if (_standard_body_forces)
      for (const auto & pr : _b)
        _b2[pr.first] = pr.second;

    // We now no longer need to store _b so we can drop its memory
    _b.clear();

    if (_standard_body_forces)
      return;

    for (auto & pr : _b2)
      pr.second /= libmesh_map_find(sf_sfhat_sum, pr.first);
  }
}
