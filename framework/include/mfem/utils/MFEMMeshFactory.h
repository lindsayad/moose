//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once
#include "libmesh/ignore_warnings.h"
#include "mfem.hpp"
#include "libmesh/restore_warnings.h"
#include "libmesh/elem.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/equation_systems.h"
#include "libmesh/face_quad4.h"
#include "libmesh/ignore_warnings.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_input.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/vtk_io.h"
#include "CubitElementInfo.h"
#include "MooseMesh.h"
#include "MFEMMesh.h"
#include <memory>
#include <map>
#include <tuple>
#include <vector>

/// Returns an MFEM mesh object corresponding to the provided
/// MooseMesh argument. If the argument is of type MFEMMesh then the
/// MFEM mesh it contains is what is returned. Otherwise the MFEM mesh
/// will be constructed based on libmesh data.
///
/// The boolean arguments only have an effect if the mesh is not
/// already an MFEMMesh. In that case, `fallback` indicates whether
/// element types which can not be represented in MFEM should be
/// represented by simpler ones which can be. `first_order` indicates
/// whether the mesh should be converted to have only first-order
/// elements.
std::shared_ptr<mfem::ParMesh> buildMFEMMesh(MooseMesh & mesh, bool fallback, bool first_order);

using IDMap = std::map<int, std::vector<int>>;

/**
 * An internal method used to create maps from each boundary ID to vectors of side IDs
 * and element IDs.
 */
std::tuple<IDMap, IDMap> buildBoundaryInfo(MooseMesh & mesh);

/**
 * Create a mapping from each boundary ID to a vector of vectors containing the global node
 * IDs of nodes that lie on the faces of elements that fall on the boundary.
 */
std::map<int, std::vector<std::vector<unsigned int>>>
buildBoundaryNodeIDs(const MooseMesh & mesh,
                     const std::vector<int> & unique_side_bound_ids,
                     const IDMap & element_ids_for_bound,
                     const IDMap & side_ids_for_bound);

/**
 * Builds two maps:
 * 1. Mapping from each block ID --> vector containing all element IDs for block.
 * 2. Mapping from each element --> vector containing all global node IDs for element.
 */
std::tuple<IDMap, IDMap> buildElementAndNodeIDs(MeshBase & libmesh,
                                                const CubitBlockInfo block_info,
                                                const std::vector<int> & unique_block_ids);

/**
 * Iterates through each block to find the elements in the block. For each
 * element in a block, it runs through the nodes in the block and adds only
 * the corner nodes to a vector. This is then sorted and only unique global
 * node IDs are retained.
 */
std::vector<int> buildUniqueCornerNodeIDs(const CubitBlockInfo & block_info,
                                          const std::vector<int> & unique_block_ids,
                                          const IDMap & element_ids_for_block_id,
                                          const IDMap & node_ids_for_element_id);

/**
 * Updates the two-way mappings between libMesh global node IDs and MFEM local DoF indices so
 * that they reflect the parallel DoF numbering of the ParMesh rather than the serial one.
 * The partitioning array (indexed by serial element index, value = owning rank) must be the
 * same array that was used to construct the ParMesh.
 */
void convertSerialDofMappingsToParallel(
    const mfem::Mesh & serial_mesh,
    const mfem::ParMesh & parallel_mesh,
    const int * partitioning,
    std::map<int, int> & libmesh_global_node_id_for_mfem_local_node_id,
    std::map<int, int> & mfem_local_node_id_for_libmesh_global_node_id);

/**
 * Assemble information on block elements .
 */
CubitBlockInfo buildCubitBlockInfo(MeshBase & libmesh,
                                   const std::vector<int> & unique_block_ids,
                                   bool fallback,
                                   bool first_order);

/**
 * Blocks/subdomains are separate subsets of the mesh that could have different
 * material properties etc. This method returns a vector containing the unique
 * IDs of each block in the mesh. This will be passed to the MFEMMesh constructor
 * which sets the attribute of each element to the ID of the block that it is a
 * part of.
 */
std::vector<int> getLibmeshBlockIDs(const MeshBase & libmesh);

/**
 * Returns a vector containing the IDs of all boundaries.
 */
std::vector<int> getSideBoundaryIDs(const MeshBase & libmesh);

/**
 * Maps from the element ID to the block ID.
 */
std::map<int, int>
getBlockIDForElementID(const std::map<int, std::vector<int>> & element_ids_for_block_id);

/**
 * Maps from the boundary ID to a vector containing the block IDs of all elements that lie on
 * the boundary.
 */
IDMap getBlockIDsForBoundaryID(const std::map<int, std::vector<int>> & element_ids_for_block_id,
                               const std::map<int, std::vector<int>> & element_ids_for_boundary_id);

#endif
