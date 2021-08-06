//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVRhieChowInterpolator.h"
#include "INSFVAttributes.h"
#include "GatherRCDataElementThread.h"
#include "GatherRCDataFaceThread.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "SystemBase.h"

#include "libmesh/mesh_base.h"
#include "libmesh/elem_range.h"

using namespace libMesh;

registerMooseObject("NavierStokesApp", INSFVRhieChowInterpolator);

InputParameters
INSFVRhieChowInterpolator::validParams()
{
  auto params = GeneralUserObject::validParams();
  params.set<ExecFlagEnum>("execute_on") = {EXEC_NONLINEAR, EXEC_LINEAR};
  params.addRequiredParam<VariableName>("u", "The x-component of velocity");
  params.addParam<VariableName>("v", "The y-component of velocity");
  params.addParam<VariableName>("w", "The z-component of velocity");
  return params;
}

INSFVRhieChowInterpolator::INSFVRhieChowInterpolator(const InputParameters & params)
  : GeneralUserObject(params),
    _moose_mesh(_subproblem.mesh()),
    _mesh(_moose_mesh.getMesh()),
    _u(_subproblem.getVariable(0, getParam<VariableName>("u"))),
    _v(isParamValid("v") ? &_subproblem.getVariable(0, getParam<VariableName>("v")) : nullptr),
    _w(isParamValid("w") ? &_subproblem.getVariable(0, getParam<VariableName>("w")) : nullptr)
{
  _var_numbers.insert(_u.number());
  if (_v)
    _var_numbers.insert(_v->number());
  if (_w)
    _var_numbers.insert(_w->number());
}

void
INSFVRhieChowInterpolator::addToA(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  if (elem->processor_id() != this->processor_id())
    _elements_to_pull.push_back(elem->id());

  _a[elem->id()](component) += value;
}

void
INSFVRhieChowInterpolator::addToB(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  mooseAssert(elem->processor_id() == this->processor_id(), "Sources should be local");

  _b[elem->id()](component) += value;
}

const VectorValue<ADReal> &
INSFVRhieChowInterpolator::rcCoeff(const libMesh::Elem * const elem) const
{
  auto it = _a.find(elem->id());
  mooseAssert(it != _a.end(), "Could not find the requested element with id " << elem->id());
  return it->second;
}

void
INSFVRhieChowInterpolator::initialize()
{
  _elements_to_pull.clear();

  _a.clear();
  _b.clear();
}

void
INSFVRhieChowInterpolator::execute()
{
  PARALLEL_TRY
  {
    GatherRCDataElementThread et(_fe_problem, _var_numbers);
    Threads::parallel_reduce(*_moose_mesh.getActiveLocalElementRange(), et);
  }
  PARALLEL_CATCH;

  PARALLEL_TRY
  {
    using FVRange = StoredRange<std::vector<const FaceInfo *>::const_iterator, const FaceInfo *>;
    GatherRCDataFaceThread<FVRange> fvr(_fe_problem, _var_numbers);
    FVRange faces(_fe_problem.mesh().faceInfo().begin(), _fe_problem.mesh().faceInfo().end());
    Threads::parallel_reduce(faces, fvr);
  }
  PARALLEL_CATCH;
}

void
INSFVRhieChowInterpolator::finalize()
{
}
