//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SparseOps.h"

#ifdef MOOSE_SPARSE_AD
thread_local SharedPool<DynamicSparseNumberArray<Real, libMesh::dof_id_type>> ad_derivatives_pool;
#endif
