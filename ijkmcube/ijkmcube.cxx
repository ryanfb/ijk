/// \file ijkmcube.cxx
/// Marching cubes/hypercubes isosurface generation

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006,2007 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <assert.h>
#include <sstream>
#include <string>
#include <vector>

#include "ijkmcube.h"
#include "ijkoctree.h"
#include "ijktable.h"

using namespace IJK;
using namespace IJKGRID;
using namespace IJKMCUBE;
using namespace IJKTABLE;

// **************************************************
// MARCHING CUBES (HYPERCUBES)
// **************************************************

namespace {

  // error checking functions
  bool check_isotable_encoding
  (const ISOSURFACE_TABLE & isotable, 
   const ISOSURFACE_TABLE::ENCODING encoding, ERROR & error);
  bool check_nep_num_dup(const int num_dup, ERROR & error);
  bool check_nep(const ISOSURFACE_TABLE & isotable, 
		 const MERGE_DATA & merge_data, const int num_dup,
		 ERROR & error);
}

void IJKMCUBE::marching_cubes
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);
  

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::marching_cubes
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord)
// same as previous function but without reporting time
// see previous function for explanation of parameter list
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  MCUBE_INFO mcube_info;
  ISO_MERGE_DATA merge_data(dimension, axis_size);
  marching_cubes(scalar_grid, isotable, isovalue, 
		 simplex_vert, vertex_coord, merge_data, mcube_info);
}

void IJKMCUBE::marching_cubes
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
 MCUBE_INFO & mcube_info)
// same as previous functions but uses minmax to extract isosurface simplices
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_from_minmax
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, minmax);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::marching_cubes
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree, 
 MCUBE_INFO & mcube_info)
// same as previous functions but uses octree to extract isosurface simplices
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_from_octree
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, octree);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

// **************************************************
// NEP MARCHING CUBES (NEP: NEGATIVE-EQUALS-POSITIVE)
// **************************************************

void IJKMCUBE::marching_cubes_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, NEP_INFO & nep_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// nep_num_dup = nep control for number of duplicate isosurface simplices
//        0 : Do not create duplicate isosurface simplices
//        1 : Merge duplicate isosurface simplices into one simplex
//        2 : Allow duplicate isosurface simplices
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// nep_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("marching_cubes_nep");

  if (!check_nep_num_dup(nep_num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  nep_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, nep_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, NEP,
     simplex_vert, vertex_coord, merge_data, nep_info);

  // store times
  clock_t t_end = clock();
  nep_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::marching_cubes_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
 NEP_INFO & nep_info)
// same as previous function but uses minmax to extract isosurface simplices
{
  PROCEDURE_ERROR error("marching_cubes_nep");

  if (!check_nep_num_dup(nep_num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  nep_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_from_minmax_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, nep_info,
     minmax);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, NEP,
     simplex_vert, vertex_coord, merge_data, nep_info);

  // store times
  clock_t t_end = clock();
  nep_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::marching_cubes_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree,
 NEP_INFO & nep_info)
// same as previous function but uses minmax to extract isosurface simplices
{
  PROCEDURE_ERROR error("marching_cubes_nep");

  if (!check_nep_num_dup(nep_num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  nep_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_from_octree_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, nep_info,
     octree);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, NEP,
     simplex_vert, vertex_coord, merge_data, nep_info);

  // store times
  clock_t t_end = clock();
  nep_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

// **************************************************
// MARCHING CUBES INTERVAL VOLUME
// **************************************************

void IJKMCUBE::MCVol
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, 
 MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  clock_t t0 = clock();
  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  VERTEX_INDEX num_mixed_cubes;
  extract_ivol_simplices
    (scalar_grid, isotable, isovalue0, isovalue1, 
     iso_simplices, num_mixed_cubes);
  mcube_info.grid.num_mixed_cubes = num_mixed_cubes;
  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> ivol_vlist;
  merge_identical_vertices
    (dimension, axis_size, iso_simplices, ivol_vlist, simplex_vert, 
     merge_data);
  clock_t t2 = clock();

  // num_ivolv = number of interval volume vertices
  int num_ivolv = ivol_vlist.size();

  int numc = num_ivolv*dimension;
  vertex_coord.resize(numc);
  position_ivol_vertices_linear
    (scalar_grid, isovalue0, isovalue1, ivol_vlist, &(vertex_coord[0]));
  clock_t t3 = clock();

  // store times
  mcube_info.time.extract = float(t1-t0)/CLOCKS_PER_SEC;
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
  mcube_info.time.total = float(t3-t_start)/CLOCKS_PER_SEC;
}

// **************************************************
// SNAPMC: MARCHING CUBES WITH QUALITY TRIANGLES
// **************************************************

void IJKMCUBE::snapMC
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,  const SNAP_TYPE snap_value,
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, 
 SNAP_INFO & snap_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// snap_value = controls snapping. 0 = no snapping. 0.5 = max snapping
// nep_num_dup = controls duplicate isosurface patches.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// snap_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  SNAP_GRID snap_grid(scalar_grid, isovalue, snap_value, 
			      snap_info.time.snap);

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  snapMC(scalar_grid, snap_grid, isotable, isovalue,
	 nep_num_dup, simplex_vert, vertex_coord, merge_data, snap_info);

  clock_t t_end = clock();

  // store time
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::snapMC
(const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, SNAP_INFO & snap_info)
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_nep
    (snap_grid, isotable, isovalue, nep_num_dup, iso_simplices, snap_info);

  merge_and_position_vertices_snapMC
    (scalar_grid, snap_grid, isovalue, iso_simplices,
     simplex_vert, vertex_coord, merge_data, snap_info);

  // store times
  clock_t t_end = clock();
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


namespace {
  void extract_iso_simplices_nep_boundary
  (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const int num_cube_facet_vertices,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices, int & num_non_empty_cubes);
}


void IJKMCUBE::snapMC_from_list
(const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 const std::vector<VERTEX_INDEX> & snap_vlist,
 const bool extract_from_boundary,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, SNAP_INFO & snap_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("snapMC_from_list");

  if (!snap_grid.MatchesSize(scalar_grid, error)) { throw error; }
  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };
  if (!merge_data.Check(NEP, error)) { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  clock_t t0 = clock();

  extract_iso_simplices_from_list
    (snap_grid, isotable, isovalue, &snap_vlist[0], snap_vlist.size(), 
     iso_simplices, snap_info);

  GRID_SIZE_TYPE num_cube_facet_vertices = 
    compute_num_cube_facet_vertices(dimension);

  snap_info.num_non_empty_boundary_facets = 0;
  if (extract_from_boundary) {
    VERTEX_INDEX num_mixed_cubes;
    extract_iso_simplices_nep_boundary
      (snap_grid, isotable, isovalue, num_cube_facet_vertices,
       iso_simplices, num_mixed_cubes);
    snap_info.num_non_empty_boundary_facets = num_mixed_cubes;
  };

  clock_t t2 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical_vertices(dimension, axis_size,
			   iso_simplices, iso_vlist, simplex_vert, merge_data);
  clock_t t3 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = iso_vlist.size();

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  VERTEX_INDEX num_snapped;
  position_snap_vertices_linear
    (scalar_grid, snap_grid.ScalarPtrConst(), snap_grid.SnapBack(),
     isovalue, iso_vlist, &(vertex_coord[0]), num_snapped);
  snap_info.num_snapped_iso_vertices = num_snapped;

  clock_t t4 = clock();

  // store times
  snap_info.time.extract = float(t2-t0)/CLOCKS_PER_SEC;
  snap_info.time.merge = float(t3-t2)/CLOCKS_PER_SEC;
  snap_info.time.position = float(t4-t3)/CLOCKS_PER_SEC;
  snap_info.time.total = float(t4-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::snapMC
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,  const SNAP_TYPE snap_value,
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_MINMAX & minmax,
 SNAP_INFO & snap_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// snap_value = controls snapping. 0 = no snapping. 0.5 = max snapping
// nep_num_dup = controls duplicate isosurface patches.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// minmax = snapMC min and max of regions
// snap_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  SNAP_GRID snap_grid(scalar_grid, isovalue, snap_value, minmax,
		      snap_info.time.snap);

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  snapMC(scalar_grid, snap_grid, isotable, isovalue,
	 nep_num_dup, simplex_vert, vertex_coord, merge_data, 
	 minmax, snap_info);

  clock_t t_end = clock();

  // store time
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::snapMC
(const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_MINMAX & minmax, SNAP_INFO & snap_info)
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_from_minmax_nep
    (snap_grid, isotable, isovalue, nep_num_dup, iso_simplices, snap_info,
     minmax);

  merge_and_position_vertices_snapMC
    (scalar_grid, snap_grid, isovalue, iso_simplices,
     simplex_vert, vertex_coord, merge_data, snap_info);

  // store times
  clock_t t_end = clock();
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::snapMC
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,  const SNAP_TYPE snap_value,
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_OCTREE & octree,
 SNAP_INFO & snap_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// snap_value = controls snapping. 0 = no snapping. 0.5 = max snapping
// nep_num_dup = controls duplicate isosurface patches.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// octree = snapMC octree
// snap_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  SNAP_GRID snap_grid(scalar_grid, isovalue, snap_value, octree,
		      snap_info.time.snap);

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  snapMC(scalar_grid, snap_grid, isotable, isovalue,
	 nep_num_dup, simplex_vert, vertex_coord, merge_data, 
	 octree, snap_info);

  clock_t t_end = clock();

  // store time
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKMCUBE::snapMC
(const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_OCTREE & octree, SNAP_INFO & snap_info)
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_from_octree_nep
    (snap_grid, isotable, isovalue, nep_num_dup, iso_simplices, snap_info,
     octree);

  merge_and_position_vertices_snapMC
    (scalar_grid, snap_grid, isovalue, iso_simplices,
     simplex_vert, vertex_coord, merge_data, snap_info);

  // store times
  clock_t t_end = clock();
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


// **************************************************
// CLASS VERTEX INCREMENT
// **************************************************

namespace {

  class VERTEX_INCREMENT {

  protected:
    int dimension;
    VERTEX_INDEX * axis_size;
    VERTEX_INDEX * axis;
    VERTEX_INDEX * cube;
    VERTEX_INDEX num_cube_vertices;
    VERTEX_INDEX * iso;
    ISO_VERTEX_INDEX num_iso_vertices;
    int num_isov_per_gridv;

  public:
    VERTEX_INCREMENT(const AXIS_SIZE_TYPE * axis_size,
		     const ISOSURFACE_TABLE & isotable,
		     const ISOTABLE_TYPE isotable_type);
    ~VERTEX_INCREMENT();

    // get functions
    int Dimension() const { return(dimension); };
    VERTEX_INDEX NumCubeVertices() const { return(num_cube_vertices); };
    ISO_VERTEX_INDEX NumIsoVertices() const { return(num_iso_vertices); };
    int NumIsoVPerGridV() const { return(num_isov_per_gridv); };
    const AXIS_SIZE_TYPE * AxisSize() const { return(axis_size); };
    VERTEX_INDEX AxisSize(const int d) const { return(axis_size[d]); };
    VERTEX_INDEX Axis(const int d) const { return(axis[d]); };
    VERTEX_INDEX Cube(const int i) const { return(cube[i]); };
    VERTEX_INDEX Iso(const int i) const { return(iso[i]); };
  };

  VERTEX_INCREMENT::VERTEX_INCREMENT
  (const AXIS_SIZE_TYPE * axis_size, 
   const ISOSURFACE_TABLE & isotable, 
   const ISOTABLE_TYPE isotable_type)
  // constructor
  {
    PROCEDURE_ERROR error("VERTEX_INCREMENT constructor");

    this->axis_size = NULL;
    axis = NULL;
    cube = NULL;
    iso = NULL;

    assert(isotable.Polyhedron().NumVertices() == 
	   compute_num_cube_vertices(isotable.Dimension()));

    dimension = isotable.Dimension();
    num_cube_vertices = isotable.Polyhedron().NumVertices();
    num_iso_vertices = isotable.NumIsosurfaceVertices();

    this->axis_size = new VERTEX_INDEX[dimension];
    axis = new VERTEX_INDEX[dimension];
    cube = new VERTEX_INDEX[num_cube_vertices];
    iso = new VERTEX_INDEX[num_iso_vertices];

    for (int d = 0; d < dimension; d++) {
      this->axis_size[d] = axis_size[d];
    };

    compute_increment(dimension, axis_size, axis);

    compute_hypercube_vertex_increment
      (dimension, axis, num_cube_vertices, cube);

    switch(isotable_type) {
    case BINARY:
      num_isov_per_gridv = get_num_iso_vertices_per_grid_vertex(dimension);
      compute_iso_vertex_increment(isotable, cube, iso);
      break;

    case NEP:
      num_isov_per_gridv = get_num_nep_iso_vertices_per_grid_vertex(dimension);
      compute_nep_iso_vertex_increment(isotable, cube, iso);
      break;

    case IVOL:
      num_isov_per_gridv = get_num_ivol_vertices_per_grid_vertex(dimension);
      compute_ivol_vertex_increment(isotable, cube, iso);
      break;

    default:
      error.AddMessage("Programming error.  Illegal isosurface table type.");
      throw error;
      break;
    }

  }

  VERTEX_INCREMENT::~VERTEX_INCREMENT()
  // destructor
  {
    dimension = 0;
    num_cube_vertices = 0;
    num_iso_vertices = 0;
    num_isov_per_gridv = 0;

    delete [] axis_size;
    delete [] axis;
    delete [] cube;
    delete [] iso;

    axis_size = NULL;
    axis = NULL;
    cube = NULL;
    iso = NULL;
  }
    
}


// **************************************************
// CLASS FACET VERTEX INCREMENT
// **************************************************

namespace {

  class FACET_VERTEX_INCREMENT:public VERTEX_INCREMENT {

  protected:
    int orth_axis;    // axis orthogonal to the facet
    int side;         // indicate side of cube
    // side = 0. Facet orthogonal to orth_axis with min orth_axis coordinates.
    // side = 1. Facet orthogonal to orth_axis with max orth_axis coordinates.

    VERTEX_INDEX * facet_iso;
    VERTEX_INDEX * facet_vertex;
    VERTEX_INDEX * cube_vertex_index;

    void ComputeFacetIso(const ISOSURFACE_TABLE & isotable,
			 const int num_cube_vertices);

  public:
    FACET_VERTEX_INCREMENT
    (const AXIS_SIZE_TYPE * axis_size, const ISOSURFACE_TABLE & isotable,
     const ISOTABLE_TYPE isotable_type, const int orth_axis, const int side);
    ~FACET_VERTEX_INCREMENT();

    // get functions
    VERTEX_INDEX FacetIso(const int i) const { return(facet_iso[i]); };
    VERTEX_INDEX FacetVertex(const int i) const { return(facet_vertex[i]); };
    VERTEX_INDEX CubeVertexIndex(const int i) const
    { return(cube_vertex_index[i]); };
    int OrthogonalAxis() const { return(orth_axis); };
    int Side() const { return(side); };

    // disable member function Iso()
    VERTEX_INDEX Iso(const int i) const;
  };

  FACET_VERTEX_INCREMENT::FACET_VERTEX_INCREMENT
    (const AXIS_SIZE_TYPE * axis_size, const ISOSURFACE_TABLE & isotable,
     const ISOTABLE_TYPE isotable_type, const int orth_axis, const int side):
    VERTEX_INCREMENT(axis_size, isotable, isotable_type)
  {
    const int dimension = isotable.Dimension();
    PROCEDURE_ERROR error("FACET_VERTEX_INCREMENT constructor");
    const ISO_VERTEX_INDEX num_iso_vertices = isotable.NumIsosurfaceVertices();

    if (isotable_type != NEP) {
      error.AddMessage("Programming error. FACET_VERTEX_INCREMENT can only be used with NEP isotable.");
      throw error;
    }

    const GRID_SIZE_TYPE num_cube_vertices = 
      compute_num_cube_vertices(dimension);
    const GRID_SIZE_TYPE num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);

    this->orth_axis = orth_axis;
    this->side = side;

    facet_iso = new VERTEX_INDEX[num_iso_vertices];
    facet_vertex = new VERTEX_INDEX[num_facet_vertices];
    cube_vertex_index = new VERTEX_INDEX[num_facet_vertices];

    compute_facet_vertex_increment
      (isotable.Polyhedron(), orth_axis, side, axis_size,
       facet_vertex, cube_vertex_index);

    ComputeFacetIso(isotable, num_cube_vertices);
  }


  FACET_VERTEX_INCREMENT::~FACET_VERTEX_INCREMENT()
  {
    delete [] facet_iso;
    facet_iso = NULL;
    delete [] facet_vertex;
    facet_vertex = NULL;
    delete [] cube_vertex_index;
    cube_vertex_index = NULL;
  };

  void FACET_VERTEX_INCREMENT::ComputeFacetIso
  (const ISOSURFACE_TABLE & isotable, const int num_cube_vertices)
  {
    const int dimension = isotable.Dimension();
    VERTEX_INDEX cube_increment2[num_cube_vertices];

    const GRID_SIZE_TYPE num_facet_vertices =
      compute_num_cube_facet_vertices(dimension);

    for (int i = 0; i < num_cube_vertices; i++) 
      { cube_increment2[i] = 0; };

    for (int j = 0; j < num_facet_vertices; j++) {
      int iv = cube_vertex_index[j];
      cube_increment2[iv] = cube[iv];

      if (side == 0)
	{ cube_increment2[iv] -= axis[orth_axis]; };

      cube_increment2[iv] = cube_increment2[iv]*NumIsoVPerGridV()+dimension;
    }

    for (int k = 0; k < isotable.NumIsosurfaceVertices(); k++) {
      facet_iso[k] = 0;
      if (isotable.IsosurfaceVertex(k).Type() == ISOSURFACE_VERTEX::VERTEX) {
	VERTEX_INDEX iv = isotable.IsosurfaceVertex(k).Face();
	facet_iso[k] = cube_increment2[iv];
      }
    }
  }
}

// **************************************************
// CLASS NEP TABLE INDEX INCREMENT
// **************************************************

namespace {

  class NEP_TABLE_INDEX_INCREMENT {

  protected:
    int dimension;
    TABLE_INDEX * positive;
    TABLE_INDEX * equals;
    VERTEX_INDEX num_cube_vertices;
  
  public:
    NEP_TABLE_INDEX_INCREMENT(const int dimension);
    ~NEP_TABLE_INDEX_INCREMENT();

    // get functions
    int Dimension() const { return(dimension); };
    TABLE_INDEX Positive(const int i) const 
    { return(positive[i]); };
    TABLE_INDEX Equals(const int i) const 
    { return(equals[i]); };
  };

  NEP_TABLE_INDEX_INCREMENT::NEP_TABLE_INDEX_INCREMENT
  (const int dimension)
  {
    positive = NULL;
    equals = NULL;

    this->dimension = dimension;
    num_cube_vertices = compute_num_cube_vertices(dimension);
    
    positive = new TABLE_INDEX[num_cube_vertices];
    equals = new TABLE_INDEX[num_cube_vertices];

    TABLE_INDEX x = 1;
    for (VERTEX_INDEX j = 0; j < num_cube_vertices; j++) {
      equals[j] = x;
      positive[j] = 2*x;
      x = 3*x;
    }
  }

  NEP_TABLE_INDEX_INCREMENT::~NEP_TABLE_INDEX_INCREMENT()
  {
    delete [] positive;
    delete [] equals;
    positive = NULL;
    equals = NULL;
    dimension = 0;
    num_cube_vertices = 0;
  }
}

// **************************************************
// LOCAL CLASSES
// **************************************************

namespace {

  // base class for grid cubes
  class GRID_CUBE_BASE {
  protected:
    int dimension;
    VERTEX_INDEX * vertex_increment;
    VERTEX_INDEX vertex0;
    VERTEX_INDEX num_vertices;

    void FreeAll();

  public:
    GRID_CUBE_BASE(const int dimension, const AXIS_SIZE_TYPE * axis_size);
    ~GRID_CUBE_BASE() { FreeAll(); };
    
    // get functions
    int NumVertices() const { return(num_vertices); };
    int Vertex0() const { return(vertex0); };
    int Vertex(const int k) const
    { return(vertex0+vertex_increment[k]); };
  };


  GRID_CUBE_BASE::GRID_CUBE_BASE
  (const int dimension, const AXIS_SIZE_TYPE * axis_size)
  {
    VERTEX_INDEX axis_increment[dimension];

    this->dimension = dimension;
    vertex0 = 0;
    num_vertices = (1L << dimension);
    vertex_increment = new VERTEX_INDEX[num_vertices];
    if (vertex_increment == NULL) {
      FreeAll();
      return;
    };

    compute_increment(dimension, axis_size, axis_increment);
    compute_hypercube_vertex_increment
      (dimension, axis_increment, num_vertices, vertex_increment);
  }

  void GRID_CUBE_BASE::FreeAll()
  {
    delete [] vertex_increment;
    vertex_increment = NULL;
    num_vertices = 0;
    vertex0 = 0;
  }

  // grid cube with set function
  class GRID_CUBE1:public GRID_CUBE_BASE {
  public:
    GRID_CUBE1(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      GRID_CUBE_BASE(dimension, axis_size) {};
    void SetVertex0(int vertex0) { this->vertex0 = vertex0; };
  };

  // grid cube with Next() function
  class GRID_CUBE2:public GRID_CUBE_BASE {
  protected:
    AXIS_SIZE_TYPE * axis_size;
    VERTEX_INDEX * axis_increment;
    GRID_COORD_TYPE * coord;
    VERTEX_INDEX num_grid_vertices;

    void FreeAll();

  public:
    GRID_CUBE2(const int dimension, const AXIS_SIZE_TYPE * axis_size);
    ~GRID_CUBE2() { FreeAll(); };
    
    // get functions
    VERTEX_INDEX NumGridVertices() const { return(num_grid_vertices); };

    // change to next cube in grid
    void Next();
  };

  GRID_CUBE2::GRID_CUBE2
  (const int dimension, const AXIS_SIZE_TYPE * axis_size):
    GRID_CUBE_BASE(dimension, axis_size)
  {
    this->axis_size = new AXIS_SIZE_TYPE[dimension];
    axis_increment = new VERTEX_INDEX[dimension];
    coord = new GRID_COORD_TYPE[dimension];

    if (this->axis_size == NULL || axis_increment == NULL ||
	vertex_increment == NULL || coord == NULL) {
      FreeAll();
      return;
    };

    for (int d = 0; d < dimension; d++) {
      this->axis_size[d] = axis_size[d];
      this->coord[d] = 0;
    };

    compute_increment(dimension, axis_size, axis_increment);
    num_grid_vertices = compute_num_grid_vertices(dimension, axis_size);
  }

  void GRID_CUBE2::Next()
  {
    // increment vertex0 and coord
    vertex0++;
    coord[0]++;
    int d = 0;
    while (d < dimension && coord[d]+1 >= axis_size[d]) {
      coord[d] = 0;
      vertex0 = vertex0 + axis_increment[d];
      d++;
      if (d < dimension) { coord[d]++; };
    }
  }

  void GRID_CUBE2::FreeAll()
  {
    GRID_CUBE_BASE::FreeAll();

    delete [] axis_size;
    axis_size = NULL;
    delete [] coord;
    coord = NULL;
    num_grid_vertices = 0;
  }

  // region
  class REGION {
    int dimension;               // grid dimension
    int index;                   // region index
    AXIS_SIZE_TYPE * esize;      // # grid edges along region
    VERTEX_INDEX vertex0;        // first vertex in region
    GRID_COORD_TYPE * coord_first; // coordinates of first vertex in region
    GRID_COORD_TYPE * coord_last;  // coordinates of last vertex in region
    AXIS_SIZE_TYPE * axis_size;    // # voxels along grid axis
    VERTEX_INDEX * axis_increment; // increment for next vertex along axis
    bool end_of_regions;         // flag end of regions

    void FreeAll();
    void Init(const int dimension, const AXIS_SIZE_TYPE * axis_size, 
	      const AXIS_SIZE_TYPE * region_esize);

  public:
    REGION(const int dimension, const AXIS_SIZE_TYPE * axis_size, 
	   const AXIS_SIZE_TYPE * region_esize);
    REGION(const int dimension, const AXIS_SIZE_TYPE * axis_size, 
	   const AXIS_SIZE_TYPE region_esize);
    ~REGION() { FreeAll(); };

    // get functions
    int Dimension() const { return(dimension); };
    bool EndOfRegions() const { return(end_of_regions); };
    int Index() const { return(index); };
    VERTEX_INDEX Vertex0() const { return(vertex0); };
    GRID_COORD_TYPE CoordFirst(const int d) const { return(coord_first[d]); };
    GRID_COORD_TYPE CoordLast(const int d) const { return(coord_last[d]); };
    AXIS_SIZE_TYPE ESize(const int d) const { return(esize[d]); };
    AXIS_SIZE_TYPE AxisSize(const int d) const { return(axis_size[d]); };
    VERTEX_INDEX AxisIncrement(const int d) const 
    { return(axis_increment[d]); };

    // Change to next region
    void Next();
  };

  // region cube
  class REGION_CUBE: public GRID_CUBE2 {

  protected:
    REGION region;
    bool end_of_cubes;           // flag end of cubes in region

    void Init();
    void FreeAll();

  public:
    REGION_CUBE(const int dimension, const AXIS_SIZE_TYPE * axis_size, 
		const AXIS_SIZE_TYPE * num_region_edges);
    REGION_CUBE
    (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
     const AXIS_SIZE_TYPE num_region_edges);


    ~REGION_CUBE() { FreeAll(); };

    // get functions
    bool EndOfRegions() const { return(region.EndOfRegions()); };
    bool EndOfCubes() const { return(end_of_cubes); };
    int RegionIndex() const { return(region.Index()); };

    // Disable Next() function
    void Next();

    // Change to next cube in region
    void NextCube();

    // Change to next region
    void NextRegion();

  };

}

// **************************************************
// LINEAR INTERPOLATION
// **************************************************

namespace {

inline void linear_interpolate
(const int dimension, const SCALAR_TYPE s0, const GRID_COORD_TYPE * coord0,
 const SCALAR_TYPE s1, const GRID_COORD_TYPE * coord1,
 const SCALAR_TYPE isovalue, COORD_TYPE * coord2)
{
  COORD_TYPE w0, w1;
  SCALAR_TYPE s_diff = s1 - s0;
  const double EPSILON = 0.00001;
  if (s_diff > EPSILON || s_diff < -EPSILON) { 
    w0 = (s1 - isovalue) / s_diff;
    w1 = (isovalue - s0) / s_diff;
  }
  else {
    // arbitrarily set weights to 0.5
    w0 = w1 = 0.5;
  };

  for (int d = 0; d < dimension; d++)
    coord2[d] = w0*coord0[d] + w1*coord1[d];
}

inline void linear_weights
(const int dimension, const SCALAR_TYPE s0, const SCALAR_TYPE s1, 
 const SCALAR_TYPE isovalue, COORD_TYPE & w0, COORD_TYPE & w1)
{
  SCALAR_TYPE s_diff = s1 - s0;
  const double EPSILON = 0.00001;
  if (s_diff > EPSILON || s_diff < -EPSILON) { 
    w0 = (s1 - isovalue) / s_diff;
    w1 = (isovalue - s0) / s_diff;
  }
  else {
    // arbitrarily set weights to 0.5
    w0 = w1 = 0.5;
  };
}

}

// **************************************************
// MARCHING CUBES SUBROUTINES
// **************************************************

// local extraction routines
namespace {

  /// Extract isosurface simplices in cube
/// Note: Make this inline for faster execution
inline void extract_iso_simplices_in_cube
(const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 int & num_simplices)
// extract isosurface simplices in cube with primary vertex iv0
// returns list representing isosurface simplices
// scalar[] = array of scalar values
//   point (x0,x1,x2,...) has scalar value
//     scalar[x0 + x1*grid_length[0] + 
//                   x2*grid_length[0]*grid_length[1] + ...]
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iv0 = primary cube vertex (cube vertex with lowest coordinates)
// increment = vertex increments
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// num_simplices = number of simplices extracted
{
  // check whether cube intersects isosurface
  bool lt_flag = true;;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] >= isovalue) {
      lt_flag = false;
      break;
    };
  };
  if (lt_flag) {
    num_simplices = 0;
    return;
  }

  bool gt_flag = true;;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] <= isovalue) {
      gt_flag = false;
      break;
    };
  };
  if (gt_flag) {
    num_simplices = 0;
    return;
  }

  long it = 0;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] >= isovalue) {
      it = (it | (1L << j));
    };
  };

  num_simplices = isotable.NumSimplices(it);
  ISO_VERTEX_INDEX isov0 = iv0*increment.NumIsoVPerGridV();
  for (int is = 0; is < num_simplices; is++) {
    for (int j = 0; j < increment.Dimension(); j++) {
      int jv = isotable.SimplexVertex(it, is, j);
      ISO_VERTEX_INDEX isov = isov0 + increment.Iso(jv);
      iso_simplices.push_back(isov);
    };
  };
};

  /// Extract isosurface simplices in cube
/// Note: Make this inline for faster execution
inline void extract_iso_simplices_in_cube_nep
(const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
 const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 int & num_simplices)
// extract isosurface simplices in cube with primary vertex iv0
//   using isotable which differentiates between scalar values
//   less than (negative), equals or greater than (positive) the isovalue
// returns list representing isosurface simplices
// scalar[] = array of scalar values
//   point (x0,x1,x2,...) has scalar value
//     scalar[x0 + x1*grid_length[0] + 
//                   x2*grid_length[0]*grid_length[1] + ...]
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iv0 = primary cube vertex (cube vertex with lowest coordinates)
// increment[] = vertex increments
// table_index_increment[] = table index increments
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// num_simplices = number of simplices extracted
{
  // check whether cube intersects isosurface
  bool lt_flag = true;;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] >= isovalue) {
      lt_flag = false;
      break;
    };
  };
  if (lt_flag) {
    num_simplices = 0;
    return;
  }

  bool gt_flag = true;;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] <= isovalue) {
      gt_flag = false;
      break;
    };
  };
  if (gt_flag) {
    num_simplices = 0;
    return;
  }

  // compute cube index
  long it = 0;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] > isovalue) {
      it += table_index_increment.Positive(j);
    }
    else if (scalar[iv1] == isovalue) {
      it += table_index_increment.Equals(j);
    };
  }

  num_simplices = isotable.NumSimplices(it);
  ISO_VERTEX_INDEX isov0 = iv0*increment.NumIsoVPerGridV();
  for (int is = 0; is < num_simplices; is++) {
    for (int j = 0; j < increment.Dimension(); j++) {
      int jv = isotable.SimplexVertex(it, is, j);
      ISO_VERTEX_INDEX isov = isov0 + increment.Iso(jv);
      iso_simplices.push_back(isov);
    };
  };
};

  /// Extract isosurface simplices in cube
/// Note: Make this inline for faster execution
inline void extract_iso_simplices_in_cube_nep
(const SCALAR_TYPE * scalar, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
 const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 int & num_simplices, std::vector<VERTEX_INDEX> & in_facet_cube)
// extract isosurface simplices in cube with primary vertex iv0
//   using isotable which differentiates between scalar values
//   less than (negative), equals or greater than (positive) the isovalue
// returns list representing isosurface simplices
// scalar[] = array of scalar values
//   point (x0,x1,x2,...) has scalar value
//     scalar[x0 + x1*grid_length[0] + 
//                   x2*grid_length[0]*grid_length[1] + ...]
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iv0 = primary cube vertex (cube vertex with lowest coordinates)
// increment[] = vertex increments
// table_index_increment[] = table index increments
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// num_simplices = number of simplices extracted
// in_facet_cube = list of cubes whose isosurface patches lie in a grid facet
{
  // check whether cube intersects isosurface
  bool lt_flag = true;;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] >= isovalue) {
      lt_flag = false;
      break;
    };
  };
  if (lt_flag) {
    num_simplices = 0;
    return;
  }

  bool gt_flag = true;;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] <= isovalue) {
      gt_flag = false;
      break;
    };
  };
  if (gt_flag) {
    num_simplices = 0;
    return;
  }

  // compute cube index
  long it = 0;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] > isovalue) {
      it += table_index_increment.Positive(j);
    }
    else if (scalar[iv1] == isovalue) {
      it += table_index_increment.Equals(j);
    };
  }

  if (isotable.IsInFacet(it))
    {
      num_simplices = 0;
      in_facet_cube.push_back(iv0);
      return;
    }

  num_simplices = isotable.NumSimplices(it);
  ISO_VERTEX_INDEX isov0 = iv0*increment.NumIsoVPerGridV();
  for (int is = 0; is < num_simplices; is++) {
    for (int j = 0; j < increment.Dimension(); j++) {
      int jv = isotable.SimplexVertex(it, is, j);
      ISO_VERTEX_INDEX isov = isov0 + increment.Iso(jv);
      iso_simplices.push_back(isov);
    };
  };
};

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_binary
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   VERTEX_INDEX & num_mixed_cubes)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_iso_simplices_in_cube
	  (scalar, isotable, isovalue, iv, increment,
	   iso_simplices, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
    };

  }

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_nep
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   VERTEX_INDEX & num_mixed_cubes)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_iso_simplices_in_cube_nep
	  (scalar, isotable, isovalue, iv, increment, table_index_increment,
	   iso_simplices, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
    };
  }

  void extract_iso_patches_in_facets
  (const SCALAR_TYPE * scalar, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   const std::vector<VERTEX_INDEX> & in_facet_cube,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   NEP_INFO & nep_info)
  // extract isosurface patches which lie on grid facets
  {
    const int dimension = isotable.Dimension();
    const int numf = isotable.Polyhedron().NumFacets();
    VERTEX_INDEX facet_increment[numf];
    std::vector<VERTEX_INDEX> facet_list;
    std::vector<VERTEX_INDEX> sorted_facet_list;
    std::vector<VERTEX_INDEX> facet_list2;

    compute_facet_increment(isotable.Polyhedron(), increment.AxisSize(),
			    facet_increment);

    for (int icube = 0; icube < in_facet_cube.size(); icube++) {
      // compute cube index
      VERTEX_INDEX iv0 = in_facet_cube[icube];
      long it = 0;
      for (int j = 0; j < increment.NumCubeVertices(); j++) {
	VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
	if (scalar[iv1] > isovalue) {
	  it += table_index_increment.Positive(j);
	}
	else if (scalar[iv1] == isovalue) {
	  it += table_index_increment.Equals(j);
	};
      }

      VERTEX_INDEX jf = isotable.ContainingFacet(it);
      VERTEX_INDEX kf = iv0*dimension + facet_increment[jf];

      facet_list.push_back(kf);
      sorted_facet_list.push_back(kf);
    }

    nep_info.num_in_facet_cubes = facet_list.size();

    std::sort(sorted_facet_list.begin(), sorted_facet_list.end());

    nep_info.num_dup_iso_patches = 0;
    int j = 0;
    while (j < sorted_facet_list.size()) {

      if (j+1 < sorted_facet_list.size() && 
	  sorted_facet_list[j] == sorted_facet_list[j+1]) {
	// skip sorted_facet_list[j] and sorted_facet_list[j+1]
	j = j+2;
	nep_info.num_dup_iso_patches++;
      }
      else {
	facet_list2.push_back(sorted_facet_list[j]);
	j++;
      };
    }

    int num_simplices = 0;
    for (int icube = 0; icube < in_facet_cube.size(); icube++) {

      VERTEX_INDEX jf = facet_list[icube];
      if (binary_search(facet_list2.begin(), facet_list2.end(), jf)) {
	VERTEX_INDEX iv0 = in_facet_cube[icube];
	extract_iso_simplices_in_cube_nep
	  (scalar, isotable, isovalue, iv0, increment, table_index_increment,
	   iso_simplices, num_simplices);
      };

      if (num_simplices > 0) { nep_info.grid.num_mixed_cubes++; }
    }

  }

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_nep_del_dup
  // extract nep iso simplices deleting duplicate isosurface patches
  (const SCALAR_TYPE * scalar, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   NEP_INFO & nep_info)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);
    std::vector<VERTEX_INDEX> in_facet_cube;

    int num_simplices = 0;
    nep_info.grid.num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_iso_simplices_in_cube_nep
	  (scalar, isotable, isovalue, iv, increment, table_index_increment,
	   iso_simplices, num_simplices, in_facet_cube);

	if (num_simplices > 0) { nep_info.grid.num_mixed_cubes++; }
      };
    };

    // handle list of cubes whose isosurface patches lie in facets
    extract_iso_patches_in_facets
      (scalar, isotable, isovalue, increment, table_index_increment,
       in_facet_cube, iso_simplices, nep_info);
  }

  void extract_iso_simplices_in_facet
  (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const int iv0,
   const int num_cube_facet_vertices, const int orth_dir, const int side, 
   const FACET_VERTEX_INCREMENT & facet_vertex_increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices, int & num_simplices)
  {
    const int dimension = isotable.Dimension();
    const int num_cube_vertices = isotable.Polyhedron().NumVertices();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();

    num_simplices = 0;

    int num_equals = 0;
    for (int k = 0; k < num_cube_facet_vertices; k++) {
      int iv = iv0 + facet_vertex_increment.FacetVertex(k);
      if (scalar[iv] > isovalue)
	{ return; }
      else if (scalar[iv] == isovalue)
	{ num_equals++; };
    }

    // if there are not enough vertices labelled '=' to form a simplex,
    //   there is no isosurface simplex in the facet
    if (num_equals < isotable.NumVerticesPerSimplex())
      return;

    long it = 0;
    for (int k = 0; k < num_cube_facet_vertices; k++) {
      int iv = iv0 + facet_vertex_increment.FacetVertex(k);

      if (scalar[iv] == isovalue) {
	int j = facet_vertex_increment.CubeVertexIndex(k);
	it += table_index_increment.Equals(j);
      }
    }

    num_simplices = isotable.NumSimplices(it);
    ISO_VERTEX_INDEX isov0 = iv0*facet_vertex_increment.NumIsoVPerGridV();
    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
	int jv = isotable.SimplexVertex(it, is, j);
	ISO_VERTEX_INDEX isov = isov0 + facet_vertex_increment.FacetIso(jv);

	iso_simplices.push_back(isov);
      }
    }
  }


  void extract_iso_simplices_nep_boundary
  (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const int num_cube_facet_vertices,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices, int & num_non_empty_cubes)
  {
    const int dimension = isotable.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

    num_non_empty_cubes = 0;


    NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);

    VERTEX_INDEX max_num_cubes_in_facet = 0;
    for (int d = 0; d < dimension; d++) {

      GRID_SIZE_TYPE num_cubes_in_facet = 
	compute_num_cubes_in_grid_facet(dimension, axis_size, d);

      if (max_num_cubes_in_facet < num_cubes_in_facet)
	{ max_num_cubes_in_facet = num_cubes_in_facet; }
    }

    VERTEX_ARRAY vlist(max_num_cubes_in_facet);

    for (int d = 0; d < dimension; d++) {
      for (int side = 0; side < 2; side++) {

	FACET_VERTEX_INCREMENT facet_vertex_increment
	  (axis_size, isotable, NEP, d, side);

	GRID_SIZE_TYPE num_cubes_in_facet = 
	  compute_num_cubes_in_grid_facet(dimension, axis_size, d);

	get_cubes_in_grid_facet(dimension, axis_size, d, side, vlist.Ptr());

	for (int j = 0; j < num_cubes_in_facet; j++) {
	  int num_simplices = 0;
	  extract_iso_simplices_in_facet
	    (scalar_grid, isotable, isovalue, vlist[j], 
	     num_cube_facet_vertices, d, side,
	     facet_vertex_increment, nep_table_index_increment,
	     iso_simplices, num_simplices);

	  if (num_simplices > 0) { num_non_empty_cubes++; }
	}
      }
    }
  }


  template <class DATASTRUCT>
  void extract_iso_simplices_from_datastruct
  (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   MCUBE_INFO & mcube_info, DATASTRUCT & datastruct)
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

    clock_t t0 = clock();

    const GRID_SIZE_TYPE num_cubes = 
      compute_num_grid_cubes(dimension, axis_size);

    VERTEX_ARRAY vlist(num_cubes);
    VERTEX_INDEX vlist_length = 0;
    get_mixed_cubes(dimension, axis_size, datastruct, isovalue,
		    vlist.Ptr(), vlist_length);

    // initialize output
    iso_simplices.clear();

    extract_iso_simplices_from_list
      (scalar_grid, isotable, isovalue, vlist.PtrConst(), vlist_length, 
       iso_simplices, mcube_info);

    clock_t t1 = clock();
    mcube_info.time.extract = float(t1-t0)/CLOCKS_PER_SEC;

    assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
  }

  template <class DATASTRUCT>
  void extract_iso_simplices_from_datastruct_nep
  (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const int nep_num_dup,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   NEP_INFO & nep_info, DATASTRUCT & datastruct)
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

    clock_t t0 = clock();

    const GRID_SIZE_TYPE num_cubes = 
      compute_num_grid_cubes(dimension, axis_size);

    VERTEX_ARRAY vlist(num_cubes);
    VERTEX_INDEX vlist_length = 0;
    get_mixed_cubes(dimension, axis_size, datastruct, isovalue,
		    vlist.Ptr(), vlist_length);

    // initialize output
    iso_simplices.clear();

    extract_iso_simplices_from_list_nep
      (scalar_grid, isotable, isovalue, nep_num_dup,
       vlist.PtrConst(), vlist_length, true, iso_simplices, nep_info);

    clock_t t1 = clock();
    nep_info.time.extract = float(t1-t0)/CLOCKS_PER_SEC;

    assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
  }

}


/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (isotable.Encoding() != ISOSURFACE_TABLE::BINARY) {
    error.AddMessage("Programming error.  Isosurface table must be BINARY.");
    throw error;
  }

  clock_t t0 = clock();

  // get cubes in facet 0
  const GRID_SIZE_TYPE num_cubes_in_facet0 = 
    compute_num_cubes_in_grid_facet0(dimension, axis_size);
  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;

  if (isotable.Encoding() == ISOSURFACE_TABLE::BINARY) {

    VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

    extract_iso_simplices_binary
      (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
       increment, iso_simplices, num_mixed_cubes);
  }
  else if (isotable.Encoding() == ISOSURFACE_TABLE::BASE3) {

    NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);
    VERTEX_INCREMENT increment(axis_size, isotable, NEP);

    extract_iso_simplices_nep
      (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
       increment, nep_table_index_increment, iso_simplices, num_mixed_cubes);
  }
  else {
    error.AddMessage("Illegal isotable encoding.");
    throw error;
  }

  mcube_info.grid.num_mixed_cubes = num_mixed_cubes;

  clock_t t1 = clock();
  mcube_info.time.extract = float(t1-t0)/CLOCKS_PER_SEC;

  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
}

/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, NEP_INFO & nep_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BASE3);

  clock_t t0 = clock();

  // get cubes in facet 0
  const GRID_SIZE_TYPE num_cubes_in_facet0 = 
    compute_num_cubes_in_grid_facet0(dimension, axis_size);
  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  nep_info.grid.num_mixed_cubes = 0;

  NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);
  VERTEX_INCREMENT increment(axis_size, isotable, NEP);

  switch (nep_num_dup) {

  case 0:
    extract_iso_simplices_nep_del_dup
      (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
       increment, nep_table_index_increment, iso_simplices, nep_info);
    break;

  case 1:
    // NOT YET IMPLEMENTED
    // DEFAULTS TO case 2

  case 2:
    {
      VERTEX_INDEX num_mixed_cubes;
      extract_iso_simplices_nep
	(scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
	 increment, nep_table_index_increment, iso_simplices, num_mixed_cubes);


      nep_info.grid.num_mixed_cubes = num_mixed_cubes;

      GRID_SIZE_TYPE num_cube_facet_vertices = 
	compute_num_cube_facet_vertices(dimension);
      extract_iso_simplices_nep_boundary
	(scalar_grid, isotable, isovalue, num_cube_facet_vertices,
	 iso_simplices, num_mixed_cubes);
      nep_info.num_non_empty_boundary_facets = num_mixed_cubes;
      break;
    }

  default:
    error.AddMessage("Programming error.  Illegal value ", nep_num_dup,
		     " for nep_num_dup.");
    error.AddMessage("   Value should be 0, 1 or 2.");
    throw error;
  }

  clock_t t1 = clock();
  nep_info.time.extract = float(t1-t0)/CLOCKS_PER_SEC;
}

/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_from_list
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const VERTEX_INDEX * vlist,
 const VERTEX_INDEX num_cubes, std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// vlist[] = list of primary cube vertices 
//   (vertices with lowest coord in cubes)
// num_cubes = number of cubes
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// mcube_info = mcube information
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices_from_list");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;
  int num_simplices = 0;

  if (isotable.Encoding() == ISOSURFACE_TABLE::BINARY) {

    VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

    for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

      extract_iso_simplices_in_cube
	(scalar, isotable, isovalue, vlist[j], increment,
	 iso_simplices, num_simplices);

      if (num_simplices > 0) { num_mixed_cubes++; }
    };

  }
  else if (isotable.Encoding() == ISOSURFACE_TABLE::BASE3) {

    NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);
    VERTEX_INCREMENT increment(axis_size, isotable, NEP);

    for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

      extract_iso_simplices_in_cube_nep
	(scalar, isotable, isovalue, vlist[j], increment, 
	 nep_table_index_increment, iso_simplices, num_simplices);

      if (num_simplices > 0) { num_mixed_cubes++; }
    };

  }
  else {
    error.AddMessage("Illegal isotable encoding.");
    throw error;
  }

  mcube_info.grid.num_mixed_cubes = num_mixed_cubes;
}

/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_from_list_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int num_dup,
 const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
 const bool extract_from_boundary,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, NEP_INFO & nep_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// vlist[] = list of primary cube vertices 
//   (vertices with lowest coord in cubes)
// num_cubes = number of cubes
// extract_from_boundary = extract simplices lying in grid boundary if true
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// nep_info = nep information
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices_from_list_nep");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_nep_num_dup(num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;
  int num_simplices = 0;

  NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);
  VERTEX_INCREMENT increment(axis_size, isotable, NEP);

  switch(num_dup) {

  case 0:
    {
      std::vector<VERTEX_INDEX> in_facet_cube;

      for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

	extract_iso_simplices_in_cube_nep
	  (scalar, isotable, isovalue, vlist[j], increment, 
	   nep_table_index_increment, iso_simplices, num_simplices,
	   in_facet_cube);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
      nep_info.grid.num_mixed_cubes = num_mixed_cubes;

      // handle list of cubes whose isosurface patches lie in facets
      extract_iso_patches_in_facets
	(scalar, isotable, isovalue, increment, nep_table_index_increment,
	 in_facet_cube, iso_simplices, nep_info);

      break;
    }

  case 1:
    // NOT YET IMPLEMENTED
    // DEFAULTS TO case 2

  case 2:
    {
      for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

	extract_iso_simplices_in_cube_nep
	  (scalar, isotable, isovalue, vlist[j], increment, 
	   nep_table_index_increment, iso_simplices, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
      nep_info.grid.num_mixed_cubes = num_mixed_cubes;
      break;
    }
  }

  if (extract_from_boundary) {
    GRID_SIZE_TYPE num_cube_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);

    VERTEX_INDEX num_mixed_cubes;
    extract_iso_simplices_nep_boundary
      (scalar_grid, isotable, isovalue, num_cube_facet_vertices,
       iso_simplices, num_mixed_cubes);
    nep_info.num_non_empty_boundary_facets = num_mixed_cubes;
  }

}


void IJKMCUBE::extract_iso_simplices_from_octree
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info,
 const IJKOCTREE::OCTREE & octree)
{
  extract_iso_simplices_from_datastruct
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, octree);
}

void IJKMCUBE::extract_iso_simplices_from_minmax
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info,
 const IJKMCUBE::MINMAX_REGIONS & minmax)
{
  extract_iso_simplices_from_datastruct
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, minmax);
}

void IJKMCUBE::extract_iso_simplices_from_octree_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, NEP_INFO & nep_info,
 const IJKOCTREE::OCTREE & octree)
{
  extract_iso_simplices_from_datastruct_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, 
     nep_info, octree);
}

void IJKMCUBE::extract_iso_simplices_from_minmax_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, NEP_INFO & nep_info,
 const MINMAX_REGIONS & minmax)
{
  extract_iso_simplices_from_datastruct_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, 
     nep_info, minmax);
}

void IJKMCUBE::merge_and_position_vertices
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 const ISOTABLE_TYPE isotable_type,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
// call merge_identical_vertices and then position_vertices_linear
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// isotable_type = type of isosurface table
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  if (isotable_type == NEP || isotable_type == IVOL) {
    if (merge_data.NumObjPerVertex() < 1) {
      error.AddMessage("Incompatible isosurface lookup table and class MERGE_DATA.");
      error.AddMessage("Number of objects per vertex in MERGE_DATA should be at least 1.");
      throw error;
    };
  };

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical_vertices(dimension, axis_size, iso_simplices, 
			   iso_vlist, simplex_vert, merge_data);
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = iso_vlist.size();

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  position_iso_vertices_linear
    (scalar_grid, isovalue, isotable_type, iso_vlist, &(vertex_coord[0]));

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}

void IJKMCUBE::merge_and_position_vertices_snapMC
(const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const SCALAR_TYPE isovalue, 
 const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, SNAP_INFO & snap_info)
// call merge_identical_vertices and then position_vertices_linear
// scalar_grid = scalar grid data
// snap_grid = snap grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// snap_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices_snapMC");

  if (!merge_data.Check(NEP, error)) { throw error; };

  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical_vertices(dimension, axis_size,
			   iso_simplices, iso_vlist, simplex_vert, merge_data);
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = iso_vlist.size();

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  VERTEX_INDEX num_snapped;
  position_snap_vertices_linear
    (scalar_grid, snap_grid.ScalarPtrConst(), snap_grid.SnapBack(),
     isovalue, iso_vlist, &(vertex_coord[0]), num_snapped);
  snap_info.num_snapped_iso_vertices = num_snapped;

  clock_t t3 = clock();

  // store times
  snap_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  snap_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}


void IJKMCUBE::merge_identical_vertices
(const int dimension, const AXIS_SIZE_TYPE * grid_length,
 const std::vector<ISO_VERTEX_INDEX> & vlist0,
 std::vector<ISO_VERTEX_INDEX> & vlist1, 
 std::vector<ISO_VERTEX_INDEX> & vlist0_map, MERGE_DATA & merge_data)
// merge identical vertices in list of isosurface vertices
// return vlist1 and vlist0_map
// dimension = volume dimension
// vlist0[] = input list of isosurface vertices
//   Precondition: Vertex indices are in range [0,merged_data.NumEdges()-1]
// vlist1[] = output list of isosurface vertices
// vlist0_map[] = vlist0[iv] is mapped to vlist1[vlist0_map[iv]]
// merge_data = internal merge data structure
{
  PROCEDURE_ERROR error("merge_identical_vertices");

  if (!merge_data.Check(error)) { throw error; }

  // initialize
  merge_data.ClearList();
  vlist1.clear();
  vlist0_map.clear();

  ISO_VERTEX_INDEX v;
  for (int i = 0; i < vlist0.size(); i++) {
    v = vlist0[i];

    if (!merge_data.InList(v)) {
      vlist1.push_back(v);
      merge_data.Insert(v);
    };
    vlist0_map.push_back(merge_data.ListLoc(v));
  };
}

namespace {

void position_iso_vertices_linear_binary
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
// compute position of isosurface vertices using linear interpolation
// scalar_grid = scalar grid data
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  const int numv = vlist.size();
  GRID_COORD_TYPE coord0[dimension];
  GRID_COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const int num_isov_per_gridv = 
    get_num_iso_vertices_per_grid_vertex(dimension);


  for (int i = 0; i < numv; i++) {
    int d = vlist[i]%num_isov_per_gridv;
    VERTEX_INDEX v0 = vlist[i]/num_isov_per_gridv;
    VERTEX_INDEX v1 = v0+axis_increment[d];

    SCALAR_TYPE s0 = scalar[v0];
    SCALAR_TYPE s1 = scalar[v1];

    compute_coord(v0, dimension, axis_size, coord0);
    compute_coord(v1, dimension, axis_size, coord1);

    linear_interpolate(dimension, s0, coord0, s1, coord1, isovalue, coord2);
    for (int d = 0; d < dimension; d++)
      coord[i*dimension+d] = coord2[d];
  }
}

void position_iso_vertices_linear_nep
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
// compute position of isosurface vertices using linear interpolation
//   vertices can have positive, negative or equals values
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  const int numv = vlist.size();
  GRID_COORD_TYPE coord0[dimension];
  GRID_COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];
  const int VERTEX_OFFSET = dimension;

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const int num_isov_per_gridv = 
    get_num_nep_iso_vertices_per_grid_vertex(dimension);

  for (int i = 0; i < numv; i++) {
    int k = vlist[i]%num_isov_per_gridv;
    VERTEX_INDEX v0 = vlist[i]/num_isov_per_gridv;

    if (k == VERTEX_OFFSET) {
      // isosurface vertex lies on grid vertex v0
      compute_coord(v0, dimension, axis_size, coord0);
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord0[d];
    }
    else {
      // isosurface vertex lies on grid edge
      VERTEX_INDEX v1 = v0+axis_increment[k];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0);
      compute_coord(v1, dimension, axis_size, coord1);

      linear_interpolate(dimension, s0, coord0, s1, coord1, isovalue, 
			 coord2);
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord2[d];
    }
  }

}

}

void IJKMCUBE::position_iso_vertices_linear
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
// compute position of isosurface vertices using linear interpolation
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_linear");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binary(scalar_grid, isovalue, vlist, coord);
    break;

  case NEP:
    position_iso_vertices_linear_nep(scalar_grid, isovalue, vlist, coord);
    break;

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }
}

// **************************************************
// INTERVAL VOLUME SUBROUTINES
// **************************************************

namespace {

/// Extract interval volume simplices in cube
/// Note: Make this inline for faster execution
inline void extract_ivol_simplices_in_cube
(const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE lower_isovalue, const SCALAR_TYPE upper_isovalue,
 const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
 std::vector<ISO_VERTEX_INDEX> & ivol_simplices, 
 int & num_simplices)
// extract isosurface simplices in cube with primary vertex iv0
// returns list representing isosurface simplices
// scalar[] = array of scalar values
//   point (x0,x1,x2,...) has scalar value
//     scalar[x0 + x1*grid_length[0] + 
//                   x2*grid_length[0]*grid_length[1] + ...]
// isotable = hypercube isosurface table for given dimension
// iv0 = primary cube vertex (cube vertex with lowest coordinates)
// increment[] = vertex increments
// ivol_simplices[] = vector of interval volume simplex vertices
//   ivol_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// num_simplices = number of simplices in interval volume patch
{
  long it = 0;
  long k = 1;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    int iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] >= upper_isovalue) {
      it += 2*k;
    }
    else if (scalar[iv1] >= lower_isovalue) {
      it += k;
    };
    k = 3*k;
  };

  num_simplices = isotable.NumSimplices(it);
  ISO_VERTEX_INDEX ivol_v0 = iv0*increment.NumIsoVPerGridV();

  for (int is = 0; is < num_simplices; is++) {
    for (int jv = 0; jv < isotable.NumVerticesPerSimplex(); jv++) {
      int isov = isotable.SimplexVertex(it, is, jv);
      ISO_VERTEX_INDEX ivol_v = ivol_v0 + increment.Iso(isov);
      ivol_simplices.push_back(ivol_v);
    };
  };
};

}

/// Extract interval volume simplices
void IJKMCUBE::extract_ivol_simplices
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
 VERTEX_INDEX & num_mixed_cubes)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// isotable = hypercube isosurface table for given dimension
// scalar_grid[] = array of scalar values
//   point (x0,x1,x2,...) has scalar value
//     scalar_grid[x0 + x1*grid_length[0] + 
//                   x2*grid_length[0]*grid_length[1] + ...]
// grid_length[i] = grid dimension i
//   # grid points = 
//      grid_length[0]*grid_length[1]*grid_length[2]*...
// isovalue = isosurface scalar value
// ivol_simplices[] = inteval volume simplex vertices
//   ivol_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of interval volume simplex is
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  const VERTEX_INDEX num_ivolv = isotable.NumIsosurfaceVertices();

  assert(isotable.Encoding() == ISOSURFACE_TABLE::BASE3);
  assert(isotable.Dimension() == isotable.SimplexDimension());
  assert(scalar_grid.Dimension() == isotable.Dimension());

  // initialize output
  ivol_simplices.clear();
  num_mixed_cubes = 0;
  int num_simplices = 0;

  SCALAR_TYPE lower_isovalue = isovalue0;
  SCALAR_TYPE upper_isovalue = isovalue1;
  if (lower_isovalue > upper_isovalue) {
    std::swap(lower_isovalue, upper_isovalue);
  };

  // get cubes in facet 0
  const GRID_SIZE_TYPE num_cubes_in_facet0 = 
    compute_num_cubes_in_grid_facet0(dimension, axis_size);
  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  VERTEX_INCREMENT increment(axis_size, isotable, IVOL);

  for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
    for (VERTEX_INDEX iv = facet_vlist[j];
	 iv < facet_vlist[j] + axis_size[0]-1; iv++) {

      extract_ivol_simplices_in_cube
	(scalar, isotable, lower_isovalue, upper_isovalue, iv,
	 increment, ivol_simplices, num_simplices);

      if (num_simplices > 0) { num_mixed_cubes++; }
    };
  }
}

/// Extract interval volume simplices
void IJKMCUBE::extract_ivol_simplices_from_list
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
 std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
 VERTEX_INDEX & num_mixed_cubes)
// extract isosurface mesh
// returns list representing isosurface simplices
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// vlist[] = list of primary cube vertices 
//   (vertices with lowest coord in cubes)
// num_cubes = number of cubes
// ivol_simplices[] = inteval volume simplex vertices
//   ivol_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of interval volume simplex is
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  const VERTEX_INDEX num_ivolv = isotable.NumIsosurfaceVertices();

  assert(isotable.Encoding() == ISOSURFACE_TABLE::BASE3);
  assert(isotable.Dimension() == isotable.SimplexDimension());
  assert(scalar_grid.Dimension() == isotable.Dimension());

  // initialize output
  ivol_simplices.clear();
  num_mixed_cubes = 0;
  int num_simplices = 0;

  SCALAR_TYPE lower_isovalue = isovalue0;
  SCALAR_TYPE upper_isovalue = isovalue1;
  if (lower_isovalue > upper_isovalue) {
    std::swap(lower_isovalue, upper_isovalue);
  };

  VERTEX_INCREMENT increment(axis_size, isotable, IVOL);

  for (VERTEX_INDEX j = 0; j < num_cubes; j++) {
    extract_ivol_simplices_in_cube
      (scalar, isotable, lower_isovalue, upper_isovalue, vlist[j],
       increment, ivol_simplices, num_simplices);

    if (num_simplices > 0) { num_mixed_cubes++; }
  }

}

void IJKMCUBE::position_ivol_vertices_linear
(const MC_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
// compute position of interval volume vertices using linear interpolation
// scalar_grid = scalar grid data
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  const int numv = vlist.size();
  GRID_COORD_TYPE coord0[dimension];
  GRID_COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];
  const int VERTEX_OFFSET = 2*dimension;

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const int num_ivolv_per_gridv = 
    get_num_ivol_vertices_per_grid_vertex(dimension);


  for (int i = 0; i < numv; i++) {
    int k = vlist[i]%num_ivolv_per_gridv;
    VERTEX_INDEX v0 = vlist[i]/num_ivolv_per_gridv;

    if (k == VERTEX_OFFSET) {
      // isosurface vertex lies on grid vertex v0
      compute_coord(v0, dimension, axis_size, coord0);
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord0[d];
    }
    else {
      // isosurface vertex lies on grid edge
      int d = k/2;
      VERTEX_INDEX v1 = v0+axis_increment[d];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0);
      compute_coord(v1, dimension, axis_size, coord1);

      if (k % 2 == 0) {
	// linear interpolate using isovalue0
	linear_interpolate(dimension, s0, coord0, s1, coord1, isovalue0, 
			   coord2);
      }
      else {
	// linear interpolate using isovalue1
	linear_interpolate(dimension, s0, coord0, s1, coord1, isovalue1, 
			   coord2);
      }
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord2[d];
    }
  }

}

// **************************************************
// SNAPMC SUBROUTINES
// **************************************************

namespace {
  inline void snap_vertex
  (const VERTEX_INDEX v0, const VERTEX_INDEX v1, const SNAP_TYPE dist,
   const SCALAR_TYPE * scalar, const SCALAR_TYPE isovalue, 
   SCALAR_TYPE * scalar_snap, VERTEX_INDEX * snap_back, 
   SNAP_TYPE * snap_dist)
  {
    if (scalar_snap[v0] != isovalue) {
      scalar_snap[v0] = isovalue;
      snap_back[v0] = v1;
      snap_dist[v0] = dist;
    }
    else if (snap_dist[v0] > dist) {
      scalar_snap[v0] = isovalue;
      snap_back[v0] = v1;
      snap_dist[v0] = dist;
    }
    else if (snap_dist[v0] < dist) {
      return;
    }
    else {
      // snap_dist[v0] == dist
      // break distance tie by using lower vertex index for snap_back[v0]
      if (v1 < snap_back[v0]) {
	snap_back[v0] = v1;
      }
    }
  }

  inline void if_close_then_snap
  (const int dimension,
   const VERTEX_INDEX v0, const SCALAR_TYPE s0,
   const VERTEX_INDEX v1,const SCALAR_TYPE s1,
   const SCALAR_TYPE * scalar, const SCALAR_TYPE isovalue,
   const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap, 
   VERTEX_INDEX * snap_back, SNAP_TYPE * snap_dist)
  {
    COORD_TYPE w0, w1;

    if ((s0 < isovalue && s1 > isovalue) ||
	(s0 > isovalue && s1 < isovalue)) {
      // isosurface intersects edge

      linear_weights(dimension, s0, s1, isovalue, w0, w1);

      if (w1 < snap_value) {

	// snap isosurface vertex to v0.  Distance = w1.
	snap_vertex(v0, v1, w1, scalar, isovalue,
		    scalar_snap, snap_back, snap_dist);
      }
      else if (w0 < snap_value) {

	// snap isosurface vertex to v1.  Distance = w1.
	snap_vertex(v1, v0, w0, scalar, isovalue,
		    scalar_snap, snap_back, snap_dist);
      }
    }
  }

  void snap_scalar_values_in_upper_facet
  (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
   const SNAP_TYPE snap_value, int orth_dir,
   SCALAR_TYPE * scalar_snap,
   VERTEX_INDEX * snap_back, SNAP_TYPE * snap_dist)
  // snap scalar values in upper facet
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    VERTEX_INDEX axis_increment[scalar_grid.Dimension()];

    if (scalar_grid.Dimension() < 2) { return; };
    if (scalar_grid.AxisSize(orth_dir) < 1) { return; };

    compute_increment(scalar_grid, axis_increment);

    const GRID_SIZE_TYPE num_facet_vertices =
      compute_num_vertices_in_grid_facet(dimension, axis_size, orth_dir);

    const VERTEX_INDEX lower2upper_increment =
      axis_increment[orth_dir]*(scalar_grid.AxisSize(orth_dir)-1);

    VERTEX_INDEX v0, v1;
    SCALAR_TYPE s0, s1;
    COORD_TYPE w0, w1;
    for (int d = 0; d < dimension; d++) {

      if (d != orth_dir) {
	const GRID_SIZE_TYPE num_ridge_vertices =
	  compute_num_vertices_in_grid_ridge
	  (dimension, axis_size, orth_dir, d);

	VERTEX_ARRAY ridge_vlist(num_ridge_vertices);
	get_vertices_in_grid_ridge(dimension, axis_size, orth_dir, d, 
				   ridge_vlist.Ptr());


	for (int i = 0; i < num_ridge_vertices; i++) {
	  // convert ridge on lower facet to ridge on upper facet
	  ridge_vlist[i] += lower2upper_increment;
	}

	for (VERTEX_INDEX iv = 0; iv < num_ridge_vertices; iv++) {

	  v0 = ridge_vlist[iv];
	  s0 = scalar[v0];

	  for (int k = 0; k < axis_size[d]-1; k++) {

	    v1 = v0 + axis_increment[d];
	    s1 = scalar[v1];

	    if_close_then_snap(dimension, v0, s0, v1, s1, scalar, isovalue,
			       snap_value, scalar_snap, snap_back, snap_dist);

	    // set v0 to v1
	    v0 = v1;
	    s0 = s1;
	  }
	}

      }
    }
  }

}


void IJKMCUBE::snap_scalar_values
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
 VERTEX_INDEX * snap_back)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  PROCEDURE_ERROR error("snape_scalar_values");

  assert(scalar != NULL && snap_back != NULL);

  if (snap_value < 0.0 || snap_value > 0.5) {
    error.AddMessage("Snap value must be in range [0.0, 0.5].");
    throw error;
  }

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();

  // copy scalar_grid to scalar_snap
  std::copy(scalar, scalar+num_grid_vertices, scalar_snap);

  if (dimension < 1) { return; };

  const int d_last = dimension-1;
  const GRID_SIZE_TYPE num_cubes_in_facet =
    compute_num_cubes_in_grid_facet(dimension, axis_size, d_last);

  VERTEX_ARRAY facet_cube_list(num_cubes_in_facet);
  get_cubes_in_grid_facet(dimension, axis_size, d_last, 0, 
			  facet_cube_list.Ptr());

  ARRAY<SNAP_TYPE> snap_dist(num_grid_vertices);

  VERTEX_INDEX v0, v1;
  SCALAR_TYPE s0, s1;
  COORD_TYPE w0, w1;

  for (int k = 0; k+1 < axis_size[d_last]; k++) {

    VERTEX_INDEX facet_increment = k*axis_increment[d_last];
    for (VERTEX_INDEX icube = 0; icube < num_cubes_in_facet; icube++) {
      v0 = facet_cube_list[icube] + facet_increment;
      s0 = scalar[v0];
      for (int d = 0; d < dimension; d++) {
	v1 = v0 + axis_increment[d];
	s1 = scalar[v1];

	if_close_then_snap(dimension, v0, s0, v1, s1, scalar, isovalue,
			   snap_value, scalar_snap, snap_back, 
			   snap_dist.Ptr());
      }
    }
  }

  // snap scalar value in upper facets
  for (int d = 0; d < dimension; d++) {
    snap_scalar_values_in_upper_facet
      (scalar_grid, isovalue, snap_value, d,
       scalar_snap, snap_back, snap_dist.Ptr());
  }
}


void IJKMCUBE::snap_scalar_values
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
 VERTEX_INDEX * snap_back, const IJKMCUBE::MINMAX_REGIONS & minmax)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  PROCEDURE_ERROR error("snape_scalar_values");

  assert(scalar != NULL && snap_back != NULL);

  if (snap_value < 0.0 || snap_value > 0.5) {
    error.AddMessage("Snap value must be in range [0.0, 0.5].");
    throw error;
  }

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();

  // copy scalar_grid to scalar_snap
  std::copy(scalar, scalar+num_grid_vertices, scalar_snap);

  const GRID_SIZE_TYPE num_cubes = 
    compute_num_grid_cubes(dimension, axis_size);

  ARRAY<SNAP_TYPE> snap_dist(num_grid_vertices);

  REGION_CUBE cube(dimension, axis_size, minmax.RegionEdgeLength());
  VERTEX_INDEX v0, v1;
  SCALAR_TYPE s0, s1;
  while (!cube.EndOfRegions()) {
    VERTEX_INDEX iregion = cube.RegionIndex();
    if (minmax.Min(iregion) < isovalue &&
	minmax.Max(iregion) >= isovalue) {

      while (!cube.EndOfCubes()) {

	v0 = cube.Vertex0();
	s0 = scalar[v0];
	for (int d = 0; d < dimension; d++) {
	  v1 = v0 + axis_increment[d];
	  s1 = scalar[v1];

	  if_close_then_snap
	    (dimension, v0, s0, v1, s1, scalar, isovalue,
	     snap_value, scalar_snap, snap_back, snap_dist.Ptr());
	}
	cube.NextCube();
      };
    };
    cube.NextRegion();
  }

  // snap scalar value in upper facets
  for (int d = 0; d < dimension; d++) {
    snap_scalar_values_in_upper_facet
      (scalar_grid, isovalue, snap_value, d, scalar_snap, snap_back, 
       snap_dist.Ptr());
  }
}

void IJKMCUBE::snap_scalar_values
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
 VERTEX_INDEX * snap_back, const SNAP_OCTREE & octree)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  PROCEDURE_ERROR error("snape_scalar_values");

  assert(scalar != NULL && snap_back != NULL);

  if (snap_value < 0.0 || snap_value > 0.5) {
    error.AddMessage("Snap value must be in range [0.0, 0.5].");
    throw error;
  }

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();

  // copy scalar_grid to scalar_snap
  std::copy(scalar, scalar+num_grid_vertices, scalar_snap);

  const GRID_SIZE_TYPE num_cubes = 
    compute_num_grid_cubes(dimension, axis_size);

  ARRAY<SNAP_TYPE> snap_dist(num_grid_vertices);

  int num_levels = octree.NumLevels();
  if (num_levels == 0) return;
  int ileaf_level = num_levels-1;

  // Depth first search through octree
  IJKOCTREE::OCTREE_STACK stack(num_levels);

  IJKOCTREE::CONST_NODE_PTR root = octree.Root();
  stack.PushRoot(root);

  VERTEX_INDEX v0, v1;
  SCALAR_TYPE s0, s1;
  while (!stack.IsEmpty()) {

    IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
    if (node->MinValue() < isovalue && isovalue <= node->MaxValue()) {

      if (stack.TopIsLeaf()) {
	// leaf node
	IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
	for (int i = 0; i < node->NumChildren(); i++) {

	  v0 = octree.ComputeChildVertex0(ileaf_level, node, i);
	  s0 = scalar[v0];
	  for (int d = 0; d < dimension; d++) {
	    v1 = v0 + axis_increment[d];
	    s1 = scalar[v1];

	    if_close_then_snap
	      (dimension, v0, s0, v1, s1, scalar, isovalue,
	       snap_value, scalar_snap, snap_back, snap_dist.Ptr());
	  }

	}

	// pop stack
	stack.Pop();
      }
      else {
	// internal node
	if (stack.NoNextChild()) {
	  // pop stack
	  stack.Pop();
	}
	else {
	  // push child_node onto stack
	  stack.PushChild();
	}
      }
    }
    else {
      // pop stack
      stack.Pop();
    }
  }

  // snap scalar value in upper facets
  for (int d = 0; d < dimension; d++) {
    snap_scalar_values_in_upper_facet
      (scalar_grid, isovalue, snap_value, d,
       scalar_snap, snap_back, snap_dist.Ptr());
  }
}

void IJKMCUBE::position_snap_vertices_linear
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE * scalar_snap, 
 const VERTEX_INDEX * snap_back, const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist, 
 COORD_TYPE * coord, VERTEX_INDEX & num_snapped)
// compute position of isosurface vertices using linear interpolation
//   vertices can have positive, negative or equals values
// scalar_grid = scalar grid data
// scalar_snap[] = array of snapped values
// snap_back[i] = move i'th isosurface vertex toward snap_back[i]
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
// num_snapped = number of snapped isosurface vertices
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  const int numv = vlist.size();
  GRID_COORD_TYPE coord0[dimension];
  GRID_COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];
  const int VERTEX_OFFSET = dimension;

  num_snapped = 0;

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const int num_isov_per_gridv = 
    get_num_nep_iso_vertices_per_grid_vertex(dimension);

  for (int i = 0; i < numv; i++) {
    int k = vlist[i]%num_isov_per_gridv;
    VERTEX_INDEX v0 = vlist[i]/num_isov_per_gridv;

    if (k == VERTEX_OFFSET) {
      if (scalar_snap[v0] == scalar[v0]) {
	// isosurface vertex lies on grid vertex v0

	compute_coord(v0, dimension, axis_size, coord0);
	for (int d = 0; d < dimension; d++)
	  coord[i*dimension+d] = coord0[d];
      }
      else {
	// isosurface vertex was snapped to grid vertex v0
	// move isosurface vertex to interior of edge [v0, v1]

	VERTEX_INDEX v1 = snap_back[v0];

	SCALAR_TYPE s0 = scalar[v0];
	SCALAR_TYPE s1 = scalar[v1];

	compute_coord(v0, dimension, axis_size, coord0);
	compute_coord(v1, dimension, axis_size, coord1);

	linear_interpolate(dimension, s0, coord0, s1, coord1, isovalue, 
			   coord2);
	for (int d = 0; d < dimension; d++)
	  coord[i*dimension+d] = coord2[d];

	num_snapped++;
      }
    }
    else {
      // isosurface vertex lies on grid edge
      VERTEX_INDEX v1 = v0+axis_increment[k];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0);
      compute_coord(v1, dimension, axis_size, coord1);

      linear_interpolate(dimension, s0, coord0, s1, coord1, isovalue, 
			 coord2);
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord2[d];
    }
  }

}


// **************************************************
// UTILITY FUNCTIONS
// **************************************************

/// Return index of vertex with specified coord
/// Redefine compute_vertex_index template as function
/// returning type VERTEX_INDEX
VERTEX_INDEX IJKMCUBE::compute_vertex_index
(const GRID_COORD_TYPE * coord, const int dimension, 
 const AXIS_SIZE_TYPE * axis_size)
{
  VERTEX_INDEX iv = 
    IJKGRID::compute_vertex_index<VERTEX_INDEX>(coord, dimension, axis_size);
  return(iv);
}

void IJKMCUBE::compute_hypercube_vertex_increment
(const int dimension, const VERTEX_INDEX * axis_increment, 
 const VERTEX_INDEX num_hypercube_vertices, VERTEX_INDEX * increment)
  // compute increment to add to vertex 0 to get vertex i of hypercube
  // increment[i] = increment for vertex i
  // Precondition: array increment is allocated with size at least
  //               num_hypercube_vertices
{
  assert(axis_increment != NULL);
  assert(increment != NULL);

  for (long j = 0; j < num_hypercube_vertices; j++) {
    increment[j] = 0;
    long j0 = j;
    for (int d = 0; d < dimension; d++) {
      if ((j0 % 2) == 1) {
	increment[j] = increment[j] + axis_increment[d];
      };
      j0 = j0/2;
    };
  }
}

void IJKMCUBE::compute_iso_vertex_increment
(const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
 VERTEX_INDEX * increment)
  // compute increment to add to primary vertex to get index 
  //   of isosurface vertex i of hypercube
  // vertex_increment[k] = vertex increment for hypercube vertex k
  // increment[i] = increment for isosurface vertex i
  //              = index of edge containing vertex i
  // Precondition: array increment is allocated with size at least
  //               isotable.NumIsosurfaceVertices()
{
  PROCEDURE_ERROR error("IJKMCUBE::compute_iso_vertex_increment");

  assert(vertex_increment != NULL);
  assert(increment != NULL);
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BINARY);
  assert(isotable.NumIsosurfaceVertices() == isotable.Polyhedron().NumEdges());

  const ISO_VERTEX_INDEX num_iso_vertices = 
    isotable.NumIsosurfaceVertices();
  const int dimension = isotable.Dimension();
  const int num_isov_per_gridv = 
    get_num_iso_vertices_per_grid_vertex(dimension);

  for (long j = 0; j < num_iso_vertices; j++) {

    if (isotable.IsosurfaceVertex(j).Type() != ISOSURFACE_VERTEX::EDGE) {
      error.AddMessage("Illegal isosurface vertex ", j, ".");
      error.AddMessage("Isosurface vertex type must be EDGE.");
      throw error;
    };

    IJKTABLE::EDGE_INDEX ie = isotable.IsosurfaceVertex(j).Face();
    int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
    int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

    int d = 0;
    while (d < dimension && 
	   isotable.Polyhedron().VertexCoord(iv0, d) ==
	   isotable.Polyhedron().VertexCoord(iv1, d)) { d++; };

    // error checking
    if (d == dimension) {
      error.AddMessage("Illegal cube edge ", ie, ".");
      error.AddMessage("Endpoints ", iv0, " and ", iv1, 
		       " have the exact same coordinates.");
      throw error;
    };

    for (int d2 = 0; d2 < dimension; d2++) {
      if (d2 != d && 
	  isotable.Polyhedron().VertexCoord(iv0, d2) !=
	  isotable.Polyhedron().VertexCoord(iv1, d2)) { 
	error.AddMessage("Illegal cube edge ", ie, ".");
	error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			 " differ in more than one coordinate.");
	throw error;
      };
    }

    if (isotable.Polyhedron().VertexCoord(iv0, d) >
	isotable.Polyhedron().VertexCoord(iv1, d)) { std::swap(iv0, iv1); };

    increment[j] = vertex_increment[iv0]*num_isov_per_gridv + d;
  }
}

void IJKMCUBE::compute_nep_iso_vertex_increment
(const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
 VERTEX_INDEX * increment)
  // compute increment to add to primary vertex to get index 
  //   of isosurface vertex i of hypercube
  // vertex_increment[k] = vertex increment for hypercube vertex k
  // increment[i] = increment for isosurface vertex i
  //              = index of edge containing vertex i
  // Precondition: array increment is allocated with size at least
  //               isotable.NumIsosurfaceVertices()
{
  PROCEDURE_ERROR error("IJKMCUBE::compute_nep_iso_vertex_increment");

  assert(vertex_increment != NULL);
  assert(increment != NULL);
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BASE3);
  assert(isotable.NumIsosurfaceVertices() == 
	 isotable.Polyhedron().NumEdges()+isotable.Polyhedron().NumVertices());

  const ISO_VERTEX_INDEX num_iso_vertices = 
    isotable.NumIsosurfaceVertices();
  const int dimension = isotable.Dimension();
  const int num_isov_per_gridv = 
    get_num_nep_iso_vertices_per_grid_vertex(dimension);

  for (long j = 0; j < num_iso_vertices; j++) {

    if (isotable.IsosurfaceVertex(j).Type() == ISOSURFACE_VERTEX::EDGE) {

      IJKTABLE::EDGE_INDEX ie = isotable.IsosurfaceVertex(j).Face();
      int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
      int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

      int d = 0;
      while (d < dimension && 
	     isotable.Polyhedron().VertexCoord(iv0, d) ==
	     isotable.Polyhedron().VertexCoord(iv1, d)) { d++; };

      // error checking
      if (d == dimension) {
	error.AddMessage("Illegal cube edge ", ie, ".");
	error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			 " have the exact same coordinates.");
	throw error;
      };

      for (int d2 = 0; d2 < dimension; d2++) {
	if (d2 != d && 
	    isotable.Polyhedron().VertexCoord(iv0, d2) !=
	    isotable.Polyhedron().VertexCoord(iv1, d2)) { 
	  error.AddMessage("Illegal cube edge ", ie, ".");
	  error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			   " differ in more than one coordinate.");
	  throw error;
	};
      }

      if (isotable.Polyhedron().VertexCoord(iv0, d) >
	  isotable.Polyhedron().VertexCoord(iv1, d)) { std::swap(iv0, iv1); };

      increment[j] = vertex_increment[iv0]*num_isov_per_gridv + d;
    }
    else if (isotable.IsosurfaceVertex(j).Type() == 
	     ISOSURFACE_VERTEX::VERTEX) {
      VERTEX_INDEX iv = isotable.IsosurfaceVertex(j).Face();
      increment[j] = vertex_increment[iv]*num_isov_per_gridv + dimension;
    }
    else {
      error.AddMessage("Illegal isosurface vertex ", j, ".");
      error.AddMessage("Isosurface vertex type must be EDGE or VERTEX.");
      throw error;
    }
  }
}

void IJKMCUBE::compute_ivol_vertex_increment
(const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
 ISO_VERTEX_INDEX * increment)
  // compute increment to add to primary vertex to get index 
  //   of interval volume vertex i of hypercube
  // vertex_increment[k] = vertex increment for hypercube vertex k
  // increment[i] = increment for isosurface vertex i
  //              = index of edge containing vertex i
  // Precondition: array increment is allocated with size at least
  //               isotable.NumIsosurfaceVertices()
{
  PROCEDURE_ERROR error("IJKMCUBE::compute_ivol_vertex_increment");

  assert(vertex_increment != NULL);
  assert(increment != NULL);
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BASE3);
  assert(isotable.NumIsosurfaceVertices() == 
	 2*isotable.Polyhedron().NumEdges()+
	 isotable.Polyhedron().NumVertices());

  const ISO_VERTEX_INDEX num_ivol_vertices = 
    isotable.NumIsosurfaceVertices();
  const int dimension = isotable.Dimension();
  const int num_ivolv_per_gridv = 
    get_num_ivol_vertices_per_grid_vertex(dimension);

  for (long j = 0; j < num_ivol_vertices; j++) {

    if (isotable.IsosurfaceVertex(j).Type() == ISOSURFACE_VERTEX::EDGE) {

      IJKTABLE::EDGE_INDEX ie = isotable.IsosurfaceVertex(j).Face();
      int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
      int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

      int d = 0;
      while (d < dimension && 
	     isotable.Polyhedron().VertexCoord(iv0, d) ==
	     isotable.Polyhedron().VertexCoord(iv1, d)) { d++; };

      // error checking
      if (d == dimension) {
	error.AddMessage("Illegal cube edge ", ie, ".");
	error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			 " have the exact same coordinates.");
	throw error;
      };

      for (int d2 = 0; d2 < dimension; d2++) {
	if (d2 != d && 
	    isotable.Polyhedron().VertexCoord(iv0, d2) !=
	    isotable.Polyhedron().VertexCoord(iv1, d2)) { 
	  error.AddMessage("Illegal cube edge ", ie, ".");
	  error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			   " differ in more than one coordinate.");
	  throw error;
	};
      }

      if (isotable.Polyhedron().VertexCoord(iv0, d) >
	  isotable.Polyhedron().VertexCoord(iv1, d)) { std::swap(iv0, iv1); };

      if (isotable.IsosurfaceVertex(j).IsLabelSet()) {
	if (isotable.IsosurfaceVertex(j).Label() == "0") {
	  increment[j] = vertex_increment[iv0]*num_ivolv_per_gridv + 2*d;
	}
	else if (isotable.IsosurfaceVertex(j).Label() == "1") {
	  increment[j] = vertex_increment[iv0]*num_ivolv_per_gridv + 2*d+1;
	}
	else {
	  error.AddMessage("Interval volume vertex ", j, " has improper label.");
	  error.AddMessage("Interval volume vertices on grid edges must have labels 0 or 1.");
	  throw error;
	};
      }
      else {
	error.AddMessage("Interval volume vertex ", j, " has no label.");
	error.AddMessage("Interval volume vertices on grid edges must have labels 0 or 1.");
	throw error;
      };

    }
    else if (isotable.IsosurfaceVertex(j).Type() == 
	     ISOSURFACE_VERTEX::VERTEX) {
      VERTEX_INDEX iv = isotable.IsosurfaceVertex(j).Face();
      increment[j] = vertex_increment[iv]*num_ivolv_per_gridv + 2*dimension;
    } else {
      error.AddMessage("Interval volume vertex ", j, " has improper type.");
      error.AddMessage("Interval volume vertices must lie on grid vertices or grid edges.");
      throw error;
    };

  }

}

void IJKMCUBE::compute_hypercube_edge_increment
(const int dimension, const ISOSURFACE_TABLE & isotable, 
 const VERTEX_INDEX * hcube_vertex_increment, EDGE_INDEX * increment)
// compute increment from first edge in hypercube
// increment[i] = increment of i'th edge
//   Precondition: array increment is allocated with size at least
//                 isotable.Polyhedron().NumEdges()
{
  const int nume = isotable.Polyhedron().NumEdges();
  PROCEDURE_ERROR error("ijkmcube");

  assert(increment != NULL);

  for (int ie = 0; ie < nume; ie++) {
    int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
    int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

    if (iv0 > iv1) {
      std::swap(iv0, iv1);
    };

    bool edge_found = false;
    for (int d = 0; d < dimension; d++) {
      int jv1 = iv0 + (1L << d);
      if (jv1 == iv1) {
	increment[ie] = hcube_vertex_increment[iv0]*dimension+d;
	edge_found = true;
	break;
      }
    }

    if (!edge_found) {
      error.AddMessage("Error computing hypercube edge increments.");
      error.AddMessage
	("Hypercube vertices may not be in lexicographic order.");
      throw error;
    };
  }

}

void IJKMCUBE::compute_facet_increment
(const ISOSURFACE_TABLE_POLYHEDRON & poly, const AXIS_SIZE_TYPE * axis_size,
 VERTEX_INDEX * facet_increment)
  // compute increment to add to index of current voxel to get
  //   facet index
  // Precondition: Polyhedron is a cube.
  // Precondition: array increment is allocated with size 
  //               at least poly.NumFacets();
{
  const int dimension = poly.Dimension();
  const int numf = poly.NumFacets();
  const int numv = poly.NumVertices();
  VERTEX_INDEX axis_increment[dimension];
  int orthogonal_axis[numf];
  int axis_intercept[numf];
  int num_orthogonal[dimension];  // number of facets orthogonal to dimension
  int facet0[dimension];
  int facet1[dimension];


  PROCEDURE_ERROR error("compute_facet_increment");
  error.AddMessage("Programming error.");

  assert(numf == 2*dimension);
  assert(numv > 0 && numv == (1L << dimension));

  compute_increment(dimension, axis_size, axis_increment);

  // find orthogonal axis
  for (int jf = 0; jf < poly.NumFacets(); jf++) {
    bool is_orthogonal_axis_set = false;
    for (int d = 0; d < dimension; d++) {
      int coord;
      bool coord_is_set = false;
      bool is_orthogonal = true;
      for (int iv = 0; iv < numv; iv++) {
	if (poly.FacetVertexFlag(jf, iv)) {
	  if (coord_is_set) {
	    if (poly.VertexCoord(iv, d) != coord) {
	      is_orthogonal = false;
	      break;
	    }
	  }
	  else {
	    coord = poly.VertexCoord(iv, d); 
	    coord_is_set = true;
	  }
	}
      }

      if (coord_is_set && is_orthogonal) {
	orthogonal_axis[jf] = d;
	axis_intercept[jf] = coord;
	is_orthogonal_axis_set = true;
	break;
      }

      if (!coord_is_set) {
	error.AddMessage("No vertices lie on polyhedron facet ", jf, ".");
	throw error;
      }
    }

    if (!is_orthogonal_axis_set) {
      error.AddMessage("Unable to find axis orthogonal to polyhedron facet ", 
		       jf, ".");
      throw error;
    }
  }

  // initialize num_orthogonal[d] to 0
  for (int d = 0; d < dimension; d++) {
    num_orthogonal[d] = 0;
  }

  // set facet0 and facet1
  for (int jf = 0; jf < poly.NumFacets(); jf++) {
    int d = orthogonal_axis[jf];
    if (num_orthogonal[d] == 0) {
      facet0[d] = jf;
    }
    else if (num_orthogonal[d] == 1) {
      facet1[d] = jf;
    }
    else {
      error.AddMessage("More than two facets are orthogonal to axis ", 
		       d, ".");
      throw error;
    }

    num_orthogonal[d]++;
  };

  for (int d = 0; d < dimension; d++) {

    if (num_orthogonal[d] != 2) {
      error.AddMessage("Axis ", d, " does not have two orthogonal facets.");
      throw error;
    }

    int jf0 = facet0[d];
    int jf1 = facet1[d];
    if (axis_intercept[jf0] > axis_intercept[jf1])
      { std::swap(jf0, jf1); };

    facet_increment[jf0] = d;
    facet_increment[jf1] = dimension*axis_increment[d] + d;
  }

}

namespace {

  void get_min_max_coord
  (const ISOSURFACE_TABLE_POLYHEDRON & poly, const int d,
   int & minc, int & maxc)
  {
    minc = maxc = 0;
    if (poly.NumVertices() < 1) { return; };

    minc = poly.VertexCoord(0, d);
    maxc = poly.VertexCoord(0, d);
    for (int iv = 1; iv < poly.NumVertices(); iv++)
      {
	int coord = poly.VertexCoord(iv, d);
	if (coord < minc) { minc = coord; };
	if (coord > maxc) { maxc = coord; };
      }
  }

}

void IJKMCUBE::compute_facet_vertex_increment
(const ISOSURFACE_TABLE_POLYHEDRON & poly, const int orth_dir, const int side,
 const AXIS_SIZE_TYPE * axis_size, 
 VERTEX_INDEX * facet_vertex_increment, VERTEX_INDEX * cube_vertex_index)
// compute increment to add to primary facet vertex
//   to get other facet vertices
// poly = isosurface table polyhedron.
//   Precondition: Polyhedron is a cube.
// orth_dir = direction orthogonal to facet. 
// side = 0. Facet orthogonal to orth_dir with min orth_dir coordinates.
// side = 1. Facet orthogonal to orth_dir with max orth_dir coordinates.
// facet_vertex_increment[k]: increment to add to primary vertex
//               to compute k'th facet vertex
// Precondition: facet_vertex increment is allocated with size 
//               num_facet_vertices
// cube_vertex_index[k]: Index of k'th facet vertex in cube
//   Precondition: cube_vertex_index is preallocated with size 
//               num_facet_vertices
{
  const int dimension = poly.Dimension();
  VERTEX_INDEX axis_increment[dimension];
  int minc[dimension];
  int maxc[dimension];
  PROCEDURE_ERROR error("compute_facet_vertex_increment");

  if (dimension == 0) { return; };

  assert(0 <= orth_dir && orth_dir < dimension);

  const GRID_SIZE_TYPE num_cube_facet_vertices = 
    compute_num_cube_facet_vertices(dimension);
  compute_increment(dimension, axis_size, axis_increment);

  for (int d = 0; d < dimension; d++) 
    { get_min_max_coord(poly, d, minc[d], maxc[d]); }

  int facet_coord = minc[orth_dir];
  if (side == 0)
    { facet_coord = maxc[orth_dir]; };

  int k = 0;
  for (int iv = 0; iv < poly.NumVertices(); iv++) {
    int coord0 = poly.VertexCoord(iv, orth_dir);
    if (coord0 == facet_coord) {
      facet_vertex_increment[k] = 0;
      cube_vertex_index[k] = iv;
      for (int d = 0; d < dimension; d++) {
	if (d != orth_dir) {
	  int coord1 = poly.VertexCoord(iv, d);
	  if (coord1 == maxc[d])
	    { facet_vertex_increment[k] += axis_increment[d]; }
	}
      }
      k++;
    }
  }

  if (k != num_cube_facet_vertices) {
    error.AddMessage("Programming error.  Failed to process correct number of facet vertices.");
    error.AddMessage("  Processed ", k, " vertices.  Number of cube facet vertices = ", num_cube_facet_vertices, ".");
    throw error;
  }
}

void IJKMCUBE::get_mixed_cubes
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MINMAX_REGIONS & minmax,  const SCALAR_TYPE isovalue, 
 VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length)
// get mixed cubes (cubes with positive and negative vertices)
// cubes are represented by primary vertices (lowest coord vertices in cubes)
// minmax = minmax data structures
// vlist[k] = list of primary cube vertices
//    Precondition: vlist[] is preallocated to length at least 
//                  number of cubes
// vlist_length = number of vertices in vlist
{
  assert(axis_size != NULL && vlist != NULL);

  const GRID_SIZE_TYPE num_cubes = 
    compute_num_grid_cubes(dimension, axis_size);

  vlist_length = 0;
  REGION_CUBE cube(dimension, axis_size, minmax.RegionEdgeLength());
  while (!cube.EndOfRegions()) {
    int iregion = cube.RegionIndex();
    if (minmax.Min(iregion) < isovalue &&
	minmax.Max(iregion) >= isovalue) {
      while (!cube.EndOfCubes()) {
	vlist[vlist_length] = cube.Vertex0();
	vlist_length++;
	cube.NextCube();
      };
    };
    cube.NextRegion();
  }

  if (num_cubes < vlist_length) {
    PROCEDURE_ERROR error("get_mixed_cubes");
    error.AddMessage("Programming error.");
    error.AddMessage
      ("Procedure get_mixed_cubes returns more cubes from minmax than exist in grid.");
  }
}

void IJKMCUBE::get_mixed_cubes
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const IJKOCTREE::OCTREE & octree,  const SCALAR_TYPE isovalue, 
 VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length)
// get mixed cubes (cubes with positive and negative vertices)
// cubes are represented by primary vertices (lowest coord vertices in cubes)
// minmax = minmax data structures
// vlist[k] = list of primary cube vertices
//    Precondition: vlist[] is preallocated to length at least 
//                  number of cubes
// vlist_length = number of vertices in vlist
{
  assert(axis_size != NULL && vlist != NULL);

  const GRID_SIZE_TYPE num_cubes = 
    compute_num_grid_cubes(dimension, axis_size);

  vlist_length = 0;

  int num_levels = octree.NumLevels();
  if (num_levels == 0) return;
  int ileaf_level = num_levels-1;

  // Depth first search through octree
  IJKOCTREE::OCTREE_STACK stack(num_levels);

  IJKOCTREE::CONST_NODE_PTR root = octree.Root();
  stack.PushRoot(root);

  while (!stack.IsEmpty()) {

    IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
    if (node->MinValue() < isovalue && isovalue <= node->MaxValue()) {

      if (stack.TopIsLeaf()) {
	// leaf node
	IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
	for (int i = 0; i < node->NumChildren(); i++) {
	  VERTEX_INDEX vertex0 = 
	    octree.ComputeChildVertex0(ileaf_level, node, i);
	  vlist[vlist_length] = vertex0;
	  vlist_length++;
	}

	// pop stack
	stack.Pop();
      }
      else {
	// internal node
	if (stack.NoNextChild()) {
	  // pop stack
	  stack.Pop();
	}
	else {
	  // push child_node onto stack
	  stack.PushChild();
	}
      }
    }
    else {
      // pop stack
      stack.Pop();
    }
  }


  if (num_cubes < vlist_length) {
    PROCEDURE_ERROR error("get_mixed_cubes");
    error.AddMessage("Programming error.");
    error.AddMessage
      ("Procedure get_mixed_cubes returns more cubes from minmax than exist in grid.");
  }
}

void IJKMCUBE::get_nonempty_snap_cubes
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, std::vector<VERTEX_INDEX> & snap_vlist)
// get all possible nonempty cubes under all possible snap values
// snap_vlist may contain some extra cubes
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  SCALAR_TYPE * scalar_snap = NULL;  // snapped scalar values
  VERTEX_INDEX * snap_back = NULL;   // snap grid vertex iv toward 
                                     //   grid vertex snapto[iv]
  const SNAP_TYPE max_snap = 0.5;

  PROCEDURE_ERROR error("get_nonempty_snap_cubes");

  snap_vlist.clear();

  // create grid of snapped scalar values
  const GRID_SIZE_TYPE num_grid_vertices = 
    compute_num_grid_vertices(dimension, axis_size);
  scalar_snap = new SCALAR_TYPE[num_grid_vertices];
  snap_back = new VERTEX_INDEX[num_grid_vertices];

  // NOTE: SHOULD REPLACE THIS BY A FAST ROUTINE WHICH SIMPLY
  //  SETS TO SCALAR VALUE ANY VERTEX WITHIN 0.5 OF AN ISOSURFACE VERTEX
  snap_scalar_values
    (scalar_grid, isovalue, max_snap, scalar_snap, snap_back);

  VERTEX_INCREMENT increment(axis_size, isotable, NEP);

  GRID_SIZE_TYPE num_cubes_in_facet0 = 
    compute_num_cubes_in_grid_facet0(dimension, axis_size);

  VERTEX_INDEX * facet_vlist = new VERTEX_INDEX[num_cubes_in_facet0];
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist);

  for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
    for (VERTEX_INDEX iv0 = facet_vlist[j]; 
	 iv0 < facet_vlist[j] + axis_size[0]-1; iv0++) {

      bool found_equals = false;
      for (int k = 0; k < increment.NumCubeVertices(); k++) {
	VERTEX_INDEX iv1 = iv0 + increment.Cube(k);
	if (scalar_snap[iv1] == isovalue) {
	  found_equals = true;
	  break;
	};
      };

      if (found_equals) {
	snap_vlist.push_back(iv0);
      }
      else {
	bool found_negative = false;
	for (int k = 0; k < increment.NumCubeVertices(); k++) {
	  VERTEX_INDEX iv1 = iv0 + increment.Cube(k);
	  if (scalar_snap[iv1] < isovalue) {
	    found_negative = true;
	    break;
	  };
	};

	bool found_positive = false;
	for (int k = 0; k < increment.NumCubeVertices(); k++) {
	  VERTEX_INDEX iv1 = iv0 + increment.Cube(k);
	  if (scalar_snap[iv1] >= isovalue) {
	    found_positive = true;
	    break;
	  };
	};

	if (found_positive && found_negative)
	  { snap_vlist.push_back(iv0); }
      }
    }
  }

  delete [] facet_vlist;
  delete [] scalar_snap;
  delete [] snap_back;
}


void IJKMCUBE::set_in_facet_iso_patches(NEP_ISOSURFACE_TABLE & isotable)
{
  const int numf = isotable.Polyhedron().NumFacets();
  const int numv_per_simplex = isotable.NumVerticesPerSimplex();

  isotable.SetIsInFacet(false);

  for (TABLE_INDEX it = 0; it < isotable.NumTableEntries(); it++) {
    int nums = isotable.NumSimplices(it);
    if (nums == 0)   // no isosurface patch for table entry it
      { continue; }

    for (int jf = 0; jf < numf; jf++) {
      bool in_facet = true;
      for (int is = 0; is < nums && in_facet; is++) {
	for (int k = 0; k < numv_per_simplex && in_facet; k++) {
	  ISOSURFACE_VERTEX_INDEX isov = isotable.SimplexVertex(it, is, k);
	  if (isotable.IsosurfaceVertex(isov).Type() ==
	      ISOSURFACE_VERTEX::VERTEX) {
	    int iv = isotable.IsosurfaceVertex(isov).Face();
	    if (!isotable.Polyhedron().FacetVertexFlag(jf, iv)) 
	      { in_facet = false; };
	  }
	  else {
	    in_facet = false;
	  }
	}
      }

      if (in_facet) {
	isotable.SetContainingFacet(it, jf);
	break;
      }
    }
  }

}

std::string IJKMCUBE::get_isotable_filename
(const ISOTABLE_TYPE isotable_type, const int dimension)
// return isotable file name base on isotable_type and dimension
{
  PROCEDURE_ERROR error("get_isotable_filename");
  std::string isotable_filename;

  switch(isotable_type){

  case BINARY:
    isotable_filename = "iso";
    break;

  case NEP:
    isotable_filename = "iso.nep";
    break;

  case IVOL:
    isotable_filename = "ivol";
    break;

  default:
    error.AddMessage("Illegal isosurface table type.");
    throw error;
  };

  isotable_filename += ".cube";

  std::ostringstream dimension_stringstream;
  dimension_stringstream << dimension;

  isotable_filename += "." + dimension_stringstream.str() + "D.xit";

  return(isotable_filename);
}


// **************************************************
// SNAP_GRID_BASE
// **************************************************

bool SNAP_GRID_BASE::MatchesSize
(const MC_GRID_BASE & scalar_grid, ERROR & error) const
{
  const AXIS_SIZE_TYPE * axis_size = AxisSize();
  const AXIS_SIZE_TYPE * scalar_axis_size = scalar_grid.AxisSize();

  if (Dimension() != scalar_grid.Dimension()) {
    error.AddMessage("Snap grid dimension (", Dimension(), 
		     ") does not equal scalar grid dimension(",
		     scalar_grid.Dimension(), ").");
    return(false);
  }

  for (int d = 0; d < Dimension(); d++) {
    if (axis_size[d] != scalar_axis_size[d]) {
      error.AddMessage
	("Snap grid axis ", d , " size does not equal scalar grid axis ", 
	 d, " size.");
      error.AddMessage
	("Snap grid axis size: ", axis_size[d], 
	 ".  Scalar grid axis size: ", scalar_axis_size[d], ".");
      return(false);
    }
  }

  return(true);
}

// class for creating snap grid information
// Note: In contrast with MC_GRID_BASE and SNAP_GRID_BASE,
//       the constructor creates memory to store the grid information
SNAP_GRID::SNAP_GRID
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  Create(scalar_grid, isovalue, snap_value);
}

SNAP_GRID::SNAP_GRID
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, float & creation_time):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  clock_t t0 = clock();
  Create(scalar_grid, isovalue, snap_value);
  clock_t t1 = clock();
  creation_time = float(t1-t0)/CLOCKS_PER_SEC;
}

SNAP_GRID::SNAP_GRID
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value,  const IJKMCUBE::MINMAX_REGIONS & minmax, 
 float & creation_time):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  clock_t t0 = clock();
  Create(scalar_grid, isovalue, snap_value, minmax);
  clock_t t1 = clock();
  creation_time = float(t1-t0)/CLOCKS_PER_SEC;
}

SNAP_GRID::SNAP_GRID
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value,  const SNAP_OCTREE & octree, 
 float & creation_time):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  clock_t t0 = clock();
  Create(scalar_grid, isovalue, snap_value, octree);
  clock_t t1 = clock();
  creation_time = float(t1-t0)/CLOCKS_PER_SEC;
}

SNAP_GRID::~SNAP_GRID()
// destructor
{
  delete [] scalar;
  scalar = NULL;
  delete [] snap_back;
  snap_back = NULL;
}

void SNAP_GRID::Create
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value)
{
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  this->scalar = new SCALAR_TYPE[num_grid_vertices];
  this->snap_back = new VERTEX_INDEX[num_grid_vertices];
  snap_scalar_values
    (scalar_grid, isovalue, snap_value, this->scalar, this->snap_back);
}

void SNAP_GRID::Create
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, const IJKMCUBE::MINMAX_REGIONS & minmax)
{
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  this->scalar = new SCALAR_TYPE[num_grid_vertices];
  this->snap_back = new VERTEX_INDEX[num_grid_vertices];
  snap_scalar_values
    (scalar_grid, isovalue, snap_value, this->scalar, this->snap_back,
     minmax);
}

void SNAP_GRID::Create
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, const SNAP_OCTREE & octree)
{
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  this->scalar = new SCALAR_TYPE[num_grid_vertices];
  this->snap_back = new VERTEX_INDEX[num_grid_vertices];
  snap_scalar_values
    (scalar_grid, isovalue, snap_value, this->scalar, this->snap_back,
     octree);
}


// **************************************************
// MCUBE TIME
// **************************************************

IJKMCUBE::MCUBE_TIME::MCUBE_TIME()
// constructor
{
  Clear();
}

void IJKMCUBE::MCUBE_TIME::Clear()
{
  preprocessing = 0.0;
  snap = 0.0;
  extract = 0.0;
  merge = 0.0;
  position = 0.0;
  total = 0.0;
}

void IJKMCUBE::MCUBE_TIME::Add(const MCUBE_TIME & mcube_time)
{
  preprocessing += mcube_time.preprocessing;
  snap += mcube_time.snap;
  extract += mcube_time.extract;
  merge += mcube_time.merge;
  position += mcube_time.position;
  total += mcube_time.total;
}

// **************************************************
// INFO CLASSES
// **************************************************

IJKMCUBE::GRID_INFO::GRID_INFO()
// constructor
{
  Clear();
}

void IJKMCUBE::GRID_INFO::Clear()
{
  num_cubes = 0;
  num_mixed_cubes = 0;
  num_bipolar_edges = 0;
}

IJKMCUBE::MCUBE_INFO::MCUBE_INFO()
{
  Clear();
}

void IJKMCUBE::MCUBE_INFO::Clear()
{
  time.Clear();
  grid.Clear();
}

void IJKMCUBE::NEP_INFO::Clear()
{
  MCUBE_INFO::Clear();
  num_in_facet_cubes = 0;
  num_dup_iso_patches = 0;
  num_non_empty_boundary_facets = 0;
}

void IJKMCUBE::SNAP_INFO::Clear()
{
  NEP_INFO::Clear();
  num_snapped_iso_vertices = 0;
}

// **************************************************
// MERGE DATA
// **************************************************

void IJKMCUBE::MERGE_DATA::Init
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MERGE_INDEX num_obj_per_vertex, const MERGE_INDEX num_obj_per_edge)
{
  this->num_obj_per_vertex = num_obj_per_vertex;
  this->num_obj_per_edge = num_obj_per_edge;
  num_vertices = compute_num_grid_vertices(dimension, axis_size);
  num_edges = dimension*num_vertices;
  vertex_id0 = num_obj_per_edge*num_edges;
  MERGE_INDEX num_obj = 
    num_obj_per_vertex*num_vertices + num_obj_per_edge*num_edges;
  INTEGER_LIST<MERGE_INDEX>::Init(num_obj);
}

bool IJKMCUBE::MERGE_DATA::Check(ERROR & error) const
{
  if (MaxInt()+1 < 
      NumObjPerVertex()*NumVertices() + NumObjPerEdge()*NumEdges()) {
    error.AddMessage("Not enough allocated memory.");
    return(false);
  };

  return(true);
}

bool IJKMCUBE::MERGE_DATA::Check
(const ISOTABLE_TYPE isotable_type, ERROR & error) const
{
  if (isotable_type == NEP) {
    if (NumObjPerVertex() < 1) {
      error.AddMessage("Programming error:  MERGE_DATA has too few object per vertex.");
      error.AddMessage
	("  NEP isosurface lookup table may position isosurface vertices ");
      error.AddMessage("  at grid vertices.");
      error.AddMessage
	("  Number of objects per vertex in MERGE_DATA should be at least 1.");

      return(false);
    }
  }

  return(true);
}


// **************************************************
// NEP_ISOSURFACE_TABLE
// **************************************************

NEP_ISOSURFACE_TABLE::NEP_ISOSURFACE_TABLE(const int d):ISOSURFACE_TABLE(d)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(d, d-1);
}

NEP_ISOSURFACE_TABLE::NEP_ISOSURFACE_TABLE
(const int dimension, const int simplex_dimension):
  ISOSURFACE_TABLE(dimension, simplex_dimension)
// constructor
// dimension = dimension of space containing isosurface.  Should be 2, 3 or 4.
// simplex_dimension = dimension of isosurface simplices
//        simplex_dimension equals dimension-1 for isosurfaces.
//        simplex_dimension equals dimension for interval volumes
{
  Init(dimension, simplex_dimension);
}


void NEP_ISOSURFACE_TABLE::Init
(const int dimension, const int simplex_dimension)
// initialization routine.  Called by constructors.
{
  ISOSURFACE_TABLE::Init(dimension, simplex_dimension);

  is_in_facet = NULL;
  containing_facet = NULL;
}

NEP_ISOSURFACE_TABLE::~NEP_ISOSURFACE_TABLE()
{
  FreeAll();
}

void NEP_ISOSURFACE_TABLE::SetIsInFacet(const bool flag)
// set is_in_facet[k] = flag, for all table entries k
{
  for (IJKTABLE::TABLE_INDEX it = 0; it < NumTableEntries(); it++)
    { is_in_facet[it] = flag; };
}

void NEP_ISOSURFACE_TABLE::SetNumTableEntries(const int num_table_entries)
{
  const char * procname = "NEP_ISOSURFACE_TABLE::SetNumTableEntries";

  ISOSURFACE_TABLE::SetNumTableEntries(num_table_entries);

  is_in_facet = new bool[num_table_entries];
  containing_facet = new IJKTABLE::FACET_INDEX[num_table_entries];

  if (is_in_facet == NULL || containing_facet == NULL)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for NEP isosurface table.");
}

void NEP_ISOSURFACE_TABLE::FreeAll()
{
  ISOSURFACE_TABLE::FreeAll();

  delete [] is_in_facet;
  is_in_facet = NULL;
  delete [] containing_facet;
  containing_facet = NULL;
}

// **************************************************
// MINMAX DATA
// **************************************************

namespace {

  REGION::REGION
  (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
   const AXIS_SIZE_TYPE * region_esize)
  {
    Init(dimension, axis_size, region_esize);
  }

  REGION::REGION
  (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
   const AXIS_SIZE_TYPE region_esize)
  {
    AXIS_SIZE_TYPE num_region_edges[dimension];

    for (int d = 0; d < dimension; d++) {
      num_region_edges[d] = region_esize;
    }

    Init(dimension, axis_size, num_region_edges);
  }

  void REGION::Init
  (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
   const AXIS_SIZE_TYPE * region_esize)
  {
    this->dimension = dimension;
    end_of_regions = false;
    index = 0;
    vertex0 = 0;

    this->axis_size = new AXIS_SIZE_TYPE[dimension];
    this->esize = new AXIS_SIZE_TYPE[dimension];
    axis_increment = new VERTEX_INDEX[dimension];
    coord_first = new GRID_COORD_TYPE[dimension];
    coord_last = new GRID_COORD_TYPE[dimension];

    for (int d = 0; d < dimension; d++) {
      this->axis_size[d] = axis_size[d];
      this->esize[d] = region_esize[d];
      coord_first[d] = 0;
      coord_last[d] = std::min(coord_first[d]+GRID_COORD_TYPE(esize[d]),
			       GRID_COORD_TYPE(axis_size[d]-1));
    }

    compute_increment(dimension, axis_size, axis_increment);
  }

  void REGION::Next()
  {
    coord_first[0] = coord_first[0] + esize[0];
    vertex0 = vertex0 + esize[0];
    coord_last[0] = std::min(coord_first[0]+GRID_COORD_TYPE(esize[0]),
			     GRID_COORD_TYPE(axis_size[0]-1));
    int d = 0;
    while (d < dimension && coord_first[d]+1 >= axis_size[d]) {
      vertex0 = vertex0 - axis_increment[d]*coord_first[d];
      coord_first[d] = 0;
      coord_last[d] = std::min(esize[d], axis_size[d]-1);
      d++;
      if (d < dimension) { 
	coord_first[d] = coord_first[d] + esize[d];
	coord_last[d] = std::min(coord_first[d]+GRID_COORD_TYPE(esize[d]), 
				 GRID_COORD_TYPE(axis_size[d]-1));
	vertex0 = vertex0 + axis_increment[d]*esize[d];
      }
      else {
	end_of_regions = true;
      }
    };
    index++;
  }

  void REGION::FreeAll() 
  {
    end_of_regions = true;
    index = 0;
    vertex0 = 0;

    delete [] axis_size;
    axis_size = NULL;
    delete [] esize;
    esize = NULL;
    delete [] axis_increment;
    axis_increment = NULL;
    delete [] coord_first;
    coord_first = NULL;
    delete [] coord_last;
    coord_last = NULL;
  }

  REGION_CUBE::REGION_CUBE
  (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
   const AXIS_SIZE_TYPE * num_region_edges) :
    GRID_CUBE2(dimension, axis_size),
    region(dimension, axis_size, num_region_edges)
  {
    Init();
  }

  REGION_CUBE::REGION_CUBE
  (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
   const AXIS_SIZE_TYPE num_region_edges) :
    GRID_CUBE2(dimension, axis_size),
    region(dimension, axis_size, num_region_edges)
  {
    Init();
  }

  void REGION_CUBE::Init()
  {
    end_of_cubes = false;
  }

  void REGION_CUBE::FreeAll() 
  {
    end_of_cubes = true;
  }

  void REGION_CUBE::NextCube()
  // increment vertex0 and coord0 to next cube in region
  // set end_of_cubes to true if no next cube
  {
    vertex0++;
    coord[0]++;
    int d = 0;
    while (d < dimension && 
	   (coord[d]+1 >= axis_size[d] || coord[d] >= region.CoordLast(d))) {
      vertex0 = vertex0 - (coord[d]-region.CoordFirst(d))*axis_increment[d];
      coord[d] = region.CoordFirst(d);
      d++;
      coord[d]++; 
      if (d < dimension) {
	coord[d] = coord[d]++;
	vertex0 = vertex0 + axis_increment[d];
      }
      else {
	end_of_cubes = true;
      }
    }
  
  }

  void REGION_CUBE::NextRegion()
  {
    region.Next();

    // reset vertex0 and coord[] to region.Vertex0() and region.CoordFirst()
    end_of_cubes = false;
    vertex0 = region.Vertex0();
    for (int d = 0; d < dimension; d++) {
      coord[d] = region.CoordFirst(d);
    }
  }

}

// **************************************************
// CLASS SNAP_MINMAX
// **************************************************

void IJKMCUBE::SNAP_MINMAX::ComputeMinMax
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const SCALAR_TYPE * scalar, const AXIS_SIZE_TYPE region_edge_length)
{
  // Offset each region by one cube to include all snapped cubes
  const AXIS_SIZE_TYPE offset_length = 1;
  MINMAX_REGIONS::ComputeMinMax(dimension, axis_size, scalar,
				region_edge_length, offset_length);
}

void IJKMCUBE::SNAP_MINMAX::ComputeMinMax
(const MC_GRID_BASE & scalar_grid, AXIS_SIZE_TYPE region_edge_length)
{
  ComputeMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
		scalar_grid.ScalarPtrConst(), region_edge_length);
}

// **************************************************
// CLASS SNAP_OCTREE
// **************************************************

void IJKMCUBE::SNAP_OCTREE::SetMinMax(const SCALAR_TYPE * scalar)
{
  if (NumLevels() == 0) { return; };

  SNAP_MINMAX minmax;

  minmax.ComputeMinMax(Dimension(), AxisSize(), scalar, 2);

  OCTREE::SetMinMax(minmax);
}

// **************************************************
// ERROR CHECKING
// **************************************************

namespace {

  bool check_isotable_encoding
  (const ISOSURFACE_TABLE & isotable, 
   const ISOSURFACE_TABLE::ENCODING encoding, ERROR & error)
  {
    if (isotable.Encoding() != encoding) {
      error.AddMessage("Illegal isosurface table encoding.");
      if (encoding == ISOSURFACE_TABLE::BINARY)
	{ error.AddMessage("  BINARY encoding required."); }
      else if (encoding == ISOSURFACE_TABLE::BASE3)
	{ error.AddMessage("  BASE3 encoding required."); }

      return(false);
    }

    return(true);
  }

  bool check_nep_num_dup(const int nep_num_dup, ERROR & error)
  {
    if (nep_num_dup < 0 || nep_num_dup > 2) {
      error.AddMessage("Programming error.  Illegal value ", nep_num_dup,
		       " for nep_num_dup.");
      error.AddMessage("   Value should be 0, 1 or 2.");

      return(false);
    }

    return(true);
  }

  bool check_nep(const ISOSURFACE_TABLE & isotable, 
		 const MERGE_DATA & merge_data,
		 const int num_dup, ERROR & error)
  {
    if (!check_nep_num_dup(num_dup, error)) { throw error; }
    if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
      { throw error; };
    if (!merge_data.Check(NEP, error)) { throw error; };
  }
}
