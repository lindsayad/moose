//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseError.h"
#include "DualReal.h"
#include "metaphysicl/dualdynamicsparsenumberarray.h"
#include "metaphysicl/metaphysicl_exceptions.h"

#ifdef MOOSE_SPARSE_AD
namespace Moose
{
inline void
derivInsert(DNDerivativeType & derivs, libMesh::dof_id_type index, Real value)
{
#ifndef NDEBUG
  try
  {
    derivs.insert(index) = value;
  }
  catch (MetaPhysicL::LogicError &)
  {
    mooseError("We should have unlimited container size");
  }
#else
  derivs.insert(index) = value;
#endif
}
}
#endif
