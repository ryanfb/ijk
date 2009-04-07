/// \file ijkmcube.h
/// Generate isosurface in arbitrary dimensions

/*
  IJK: Isosurface Jeneration Code
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

#ifndef _IJKMCUBE__
#define _IJKMCUBE_

#include <string>

#include "ijk.txx"
#include "ijkgrid.txx"

#include "ijktable.h"
#include "ijkoctree.h"

using IJKTABLE::ISOSURFACE_TABLE;
using IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON;
using IJKTABLE::TABLE_INDEX;

namespace IJKMCUBE {

  typedef float SCALAR_TYPE;
  typedef float COORD_TYPE;
  typedef float SNAP_TYPE;
  typedef int GRID_COORD_TYPE;
  typedef int AXIS_SIZE_TYPE;
  typedef int VERTEX_INDEX;
  typedef int MERGE_INDEX;
  typedef int ISO_VERTEX_INDEX;

  // vertex and edge indices must have the same type
  typedef VERTEX_INDEX EDGE_INDEX;

  // interval volume vertex types
  typedef enum { LOWER, MIDDLE, UPPER} INTERVAL_VOLUME_VERTEX_TYPE;

  // isotable type
  typedef enum { BINARY, NEP, IVOL } ISOTABLE_TYPE;

  typedef IJK::ARRAY<VERTEX_INDEX> VERTEX_ARRAY;

  typedef IJKGRID::SCALAR_GRID_BASE
    <int, AXIS_SIZE_TYPE, VERTEX_INDEX, SCALAR_TYPE> MC_GRID_BASE;
  typedef IJKGRID::SCALAR_GRID_WRAPPER
    <int, AXIS_SIZE_TYPE, VERTEX_INDEX, SCALAR_TYPE> MC_GRID_WRAPPER;
  typedef IJKGRID::SCALAR_GRID
    <int, AXIS_SIZE_TYPE, VERTEX_INDEX, SCALAR_TYPE> MC_GRID;
  typedef IJKGRID::MINMAX_REGIONS
    <int, AXIS_SIZE_TYPE, VERTEX_INDEX, SCALAR_TYPE> MINMAX_REGIONS;


// **************************************************
// DATA STRUCTURES
// **************************************************

  class MCUBE_INFO;
  class NEP_INFO;
  class SNAP_INFO;
  class MERGE_DATA;
  class NEP_ISOSURFACE_TABLE;
  class SNAP_GRID_BASE;
  class SNAP_GRID;
  class SNAP_MINMAX;
  class SNAP_OCTREE;


// **************************************************
// MARCHING CUBES (HYPERCUBES)
// **************************************************

 void marching_cubes
   (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA & merge_data, MCUBE_INFO & mcube_info);

  void marching_cubes
   (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord);

  void marching_cubes
   (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
    MCUBE_INFO & mcube_info);

  void marching_cubes
   (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree, 
    MCUBE_INFO & mcube_info);

// **************************************************
// NEP MARCHING CUBES (NEP: NEGATIVE-EQUALS-POSITIVE)
// **************************************************

  void marching_cubes_nep
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, NEP_INFO & nep_info);

  void marching_cubes_nep
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
     NEP_INFO & nep_info);

  void marching_cubes_nep
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree,
     NEP_INFO & nep_info);

// **************************************************
// MCVOL: MARCHING CUBES INTERVAL VOLUME
// **************************************************

  void MCVol
   (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord, MERGE_DATA & merge_data, 
    MCUBE_INFO & mcube_info);

// **************************************************
// SNAPMC: MARCHING CUBES WITH QUALITY TRIANGLES
// **************************************************

  void snapMC
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const SNAP_TYPE snap_value,
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

  void snapMC
    (const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

  void snapMC_from_list
    (const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     const std::vector<VERTEX_INDEX> & snap_vlist,
     const bool extract_from_boundary,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

  void snapMC
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const SNAP_TYPE snap_value,
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_MINMAX & minmax,
     SNAP_INFO & snap_info);

  void snapMC
    (const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_MINMAX & minmax,
     SNAP_INFO & snap_info);

  void snapMC
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const SNAP_TYPE snap_value,
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_OCTREE & octree,
     SNAP_INFO & snap_info);

  void snapMC
    (const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_OCTREE & octree,
     SNAP_INFO & snap_info);


// **************************************************
// MARCHING CUBES SUBROUTINES
// **************************************************

  void extract_iso_simplices
    (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  void extract_iso_simplices_from_list
    (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, 
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  void extract_iso_simplices_from_octree
    (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info,  const IJKOCTREE::OCTREE & octree);

  void extract_iso_simplices_from_minmax
    (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info, const MINMAX_REGIONS & minmax);

  void merge_identical_vertices
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const std::vector<ISO_VERTEX_INDEX> & vlist0,
     std::vector<ISO_VERTEX_INDEX> & vlist1, 
     std::vector<ISO_VERTEX_INDEX> & vlist0_map, MERGE_DATA & merge_data);

  void position_iso_vertices_linear
    (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  void merge_and_position_vertices
    (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     const ISOTABLE_TYPE isotable_type,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, MCUBE_INFO & mcube_info);

// **************************************************
// NEP MARCHING CUBES SUBROUTINES
// **************************************************

  void extract_iso_simplices_nep
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, NEP_INFO & nep_info);

  void extract_iso_simplices_from_octree_nep
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, NEP_INFO & nep_info,
     const IJKOCTREE::OCTREE & octree);

  void extract_iso_simplices_from_minmax_nep
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, NEP_INFO & nep_info,
     const MINMAX_REGIONS & minmax);

  void extract_iso_simplices_from_list_nep
    (const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     const bool extract_from_boundary,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     NEP_INFO & nep_info);

// **************************************************
// INTERVAL VOLUME SUBROUTINES
// **************************************************

  void extract_ivol_simplices
    (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
     VERTEX_INDEX & num_mixed_cubes);

  void extract_ivol_simplices_from_list
    (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
     VERTEX_INDEX & num_mixed_cubes);

  void position_ivol_vertices_linear
    (const MC_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

// **************************************************
// SNAPMC SUBROUTINES
// **************************************************

  void snap_scalar_values
    (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const SNAP_TYPE snap_value0, SCALAR_TYPE * scalar_snap,
     VERTEX_INDEX * snaptoward);

  void snap_scalar_values
    (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
     VERTEX_INDEX * snap_back, const IJKMCUBE::MINMAX_REGIONS & minmax);

  void snap_scalar_values
    (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
     VERTEX_INDEX * snap_back, const SNAP_OCTREE & octree);

  void position_snap_vertices_linear
    (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE * snap_grid, 
     const VERTEX_INDEX * snap_toward, const SCALAR_TYPE isovalue,
     const std::vector<ISO_VERTEX_INDEX> & vlist, 
     COORD_TYPE * coord, VERTEX_INDEX & num_snapped);

  void merge_and_position_vertices_snapMC
    (const MC_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const SCALAR_TYPE isovalue, 
     const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

// **************************************************
// UTILITY FUNCTIONS
// **************************************************

  /// Return index of vertex with specified coord
  /// Redefine compute_vertex_index template as function
  /// returning type VERTEX_INDEX
  VERTEX_INDEX compute_vertex_index
    (const GRID_COORD_TYPE * coord, const int dimension, 
     const AXIS_SIZE_TYPE * axis_size);

  // get grid cubes intersected by the isosurface
  void get_mixed_cubes
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const MINMAX_REGIONS & minmax,  const SCALAR_TYPE isovalue, 
     VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length);
  void get_mixed_cubes
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const IJKOCTREE::OCTREE & octree,  const SCALAR_TYPE isovalue, 
     VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length);
  void get_nonempty_snap_cubes
    (const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, std::vector<VERTEX_INDEX> & snap_vlist);

  inline int get_num_iso_vertices_per_grid_vertex(const int dimension)
    { return(dimension); };
  inline int get_num_nep_iso_vertices_per_grid_vertex(const int dimension)
    { return(dimension+1); };
  inline int get_num_ivol_vertices_per_grid_vertex(const int dimension)
    { return(2*dimension+1); };

  void compute_hypercube_vertex_increment
    (const int dimension, const VERTEX_INDEX * axis_increment, 
     const VERTEX_INDEX num_hypercube_vertices, VERTEX_INDEX * increment);
  void compute_iso_vertex_increment
    (const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * vertex_increment, ISO_VERTEX_INDEX * increment);
  void compute_nep_iso_vertex_increment
    (const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * vertex_increment, ISO_VERTEX_INDEX * increment);
  void compute_ivol_vertex_increment
    (const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * vertex_increment, ISO_VERTEX_INDEX * increment);
  void compute_hypercube_edge_increment
    (const int dimension, const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * hcube_vertex_increment, EDGE_INDEX * increment);
  void compute_facet_increment
    (const ISOSURFACE_TABLE_POLYHEDRON & poly, 
     const AXIS_SIZE_TYPE * axis_size, VERTEX_INDEX * increment);
  void compute_facet_vertex_increment
    (const ISOSURFACE_TABLE_POLYHEDRON & poly, 
     const int orth_dir, const int side,
     const AXIS_SIZE_TYPE * axis_size, VERTEX_INDEX * facet_vertex_increment,
     VERTEX_INDEX * cube_vertex_index);

  void set_in_facet_iso_patches(NEP_ISOSURFACE_TABLE & isotable);

  std::string get_isotable_filename
    (const ISOTABLE_TYPE type, const int dimension);

// **************************************************
// SNAP_GRID_BASE
// **************************************************

  /// class storing snap grid information
  /// scalar[] (inherited from MC_GRID_BASE) stores snapped scalar values
  /// Note: SNAP_GRID_BASE does not provide any way to allocate memory
  ///    for scalar[] or snap_back[]

  class SNAP_GRID_BASE:public MC_GRID_BASE {
    
  protected:
    VERTEX_INDEX * snap_back;
    // grid vertex (x0,x1,x2,...) snaps back to isosurface vertex
    //     snap_toward[x0 + x1*axis_size[0] + 
    //                   x2*axis_size[0]*axis_size[1] + ...]

  public:
    SNAP_GRID_BASE(const int dimension, const AXIS_SIZE_TYPE * axis_size):
    MC_GRID_BASE(dimension, axis_size)
      { this->snap_back = NULL; };
    ~SNAP_GRID_BASE(){ this->snap_back = NULL; };
    // Note: constructor and destructor do not allocate or free any memory

    const VERTEX_INDEX * SnapBack() const { return(snap_back); };

    bool MatchesSize(const MC_GRID_BASE & scalar_grid,
		     IJK::ERROR & error_msg) const;
  };

// class for creating snap grid information
// Note: In contrast with MC_GRID_BASE and SNAP_GRID_BASE,
//       this class creates and destroys memory to store the grid information
  class SNAP_GRID:public SNAP_GRID_BASE {

  protected:
    void Create(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
		const SNAP_TYPE snap_value);
    void Create
      (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const MINMAX_REGIONS & minmax);
    void Create
      (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const SNAP_OCTREE & octree);

  public:
    SNAP_GRID
      (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value);
    SNAP_GRID
      (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, float & creation_time);
    SNAP_GRID
      (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const MINMAX_REGIONS & minmax,
       float & creation_time);
    SNAP_GRID
      (const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const SNAP_OCTREE & octree,
       float & creation_time);
    ~SNAP_GRID();
    // Note: constructor and destructor ALLOCATE and FREE memory
  };


// **************************************************
// MCUBE TIME
// **************************************************

  class MCUBE_TIME {
    // Uses system clock() function to determine time.  
    // System clock() function should return cpu time, 
    //   but may return wall time.

  public:
    // all times are in seconds
    float preprocessing;  
      // time to create data structure for faster isosurface extraction
    float snap;       // time to snap grid scalar values 
    float extract;    // time to extract isosurface mesh
    float merge;      // time to merge identical vertices
    float position;   // time to position isosurface vertices
    float total;      // extract_time+merge_time+position_time

    MCUBE_TIME();
    void Clear();
    void Add(const MCUBE_TIME & mcube_time);
  };

// **************************************************
// GRID INFO
// **************************************************

  // Grid information
  class GRID_INFO {
  public:
    VERTEX_INDEX num_cubes;         // number of grid cubes

    VERTEX_INDEX num_mixed_cubes;   // number of mixed cubes
    //   (number of cubes with some vertex value less than the isovalue
    //    and some vertex value greater than or equal to the isovalue)

    VERTEX_INDEX num_bipolar_edges; // number of bipolar edges
    //   (number of edges with some vertex value less than the isovalue
    //    and some vertex value greater than or equal to the isovalue)

    GRID_INFO();
    void Clear();     // clear all data
  };

// **************************************************
// MCUBE INFO
// **************************************************

  class MCUBE_INFO {

  public:
    MCUBE_TIME time;
    GRID_INFO grid;

    MCUBE_INFO();
    void Clear();     // clear all data
  };

  class NEP_INFO: public MCUBE_INFO {
  public:
    
    VERTEX_INDEX num_in_facet_cubes;
    // number of cubes whose isosurface patch lies in a cube facet

    VERTEX_INDEX num_dup_iso_patches;
    // number of duplicate isosurface patches
    // each pair of duplicate patches is counted only once

    VERTEX_INDEX num_non_empty_boundary_facets;
    // number of facets containing an isosurface patch

    NEP_INFO() { Clear(); };
    void Clear();     // clear all data
  };

  class SNAP_INFO:public NEP_INFO {
  public:
    VERTEX_INDEX num_snapped_iso_vertices;
    // number of snapped isosurface vertices

    SNAP_INFO() { Clear(); };
    void Clear();     // clear all data
  };

// **************************************************
// INTEGER_LIST
// **************************************************

// list of non-negative integers
  template < class INTEGER > class INTEGER_LIST {

  protected:
    typedef std::vector<INTEGER> INTEGER_VECTOR;
    typedef typename INTEGER_VECTOR::iterator INTEGER_VECTOR_ITERATOR;

    INTEGER max_int;             // maximum integer
    INTEGER_VECTOR list;
    // NOTE: DO WE REALLY NEED TO KEEP THE LIST LOCATIONS????
    INTEGER * list_loc;          // list location
    bool * in_list_flag;         // true if integer is in list

    void Init(const INTEGER max_int)
    // Precondition: max_int is non-negative
    {
      assert(max_int >= 0);

      this->max_int = max_int;
      list_loc = new INTEGER[max_int+1];
      in_list_flag = new bool[max_int+1];

      // initialize in_list_flag[] to false
      for (INTEGER i = 0; i < max_int+1; i++) {
	in_list_flag[i] = false;
      }
    };

    INTEGER_LIST()
      // constructor for inherited classes
      {
	this->max_int = 0;
	list_loc = NULL;
	in_list_flag = NULL;
      };

  public:
    INTEGER_LIST(const INTEGER max_int) { Init(max_int); };
    // Note:  INTEGER_LIST allocates and initializes memory of size max_int
    ~INTEGER_LIST() 
      {
	delete [] list_loc;
	list_loc = NULL;
	delete [] in_list_flag;
	in_list_flag = NULL;
	max_int = 0;
      };

    // get functions
    INTEGER MaxInt() const { return(max_int); };
    bool InList(const INTEGER i) const { return(in_list_flag[i]); };
    INTEGER ListLoc(const INTEGER i) const { return(list_loc[i]); };
    INTEGER ListLength() const { return(list.size()); };
    INTEGER List(const int i) const { return(list[i]); };

    // set functions
    INTEGER Insert(const INTEGER i)
    {
      if (!InList(i)) {
	int loc = list.size();
	list.push_back(i);
	list_loc[i] = loc;
	in_list_flag[i] = true;
	return(loc);
      }
      else {
	return(ListLoc(i));
      }
    };

    void ClearList()
    {
      // set in_list_flag[] to false
      for (INTEGER_VECTOR_ITERATOR iter = list.begin(); 
	   iter != list.end(); iter++) {
	int i = *iter;
	in_list_flag[i] = false;
      }

      list.clear();
    };

  };

// **************************************************
// MERGE DATA
// **************************************************

/// Internal data structure for merge_identical_vertices 
  class MERGE_DATA: public INTEGER_LIST<MERGE_INDEX> {

  protected:
    MERGE_INDEX num_edges;             // number of edges
    MERGE_INDEX num_vertices;          // number of vertices
    MERGE_INDEX num_obj_per_vertex;    // number of objects per vertex
    MERGE_INDEX num_obj_per_edge;      // number of objects per edge
    MERGE_INDEX vertex_id0;            // first vertex identifier

    void Init(const int dimension, const AXIS_SIZE_TYPE * axis_size,
	      const MERGE_INDEX num_obj_per_vertex,
	      const MERGE_INDEX num_obj_per_edge);

  public:
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size)
      { Init(dimension, axis_size, 0, 1); };
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size,
	       const MERGE_INDEX num_obj_per_vertex, 
	       const MERGE_INDEX num_obj_per_edge)
      { Init(dimension, axis_size, num_obj_per_vertex, num_obj_per_edge); };

    // get functions
    MERGE_INDEX NumEdges() const { return(num_edges); };
    MERGE_INDEX NumVertices() const { return(num_vertices); };
    MERGE_INDEX NumObjPerVertex() const { return(num_obj_per_vertex); };
    MERGE_INDEX NumObjPerEdge() const { return(num_obj_per_edge); };
    MERGE_INDEX VertexIdentifier(const MERGE_INDEX iv) const 
    { return(vertex_id0 + iv); };
    MERGE_INDEX VertexIdentifier
      (const MERGE_INDEX iv, const MERGE_INDEX j) const 
    { return(vertex_id0 + iv + j * num_vertices); };
    MERGE_INDEX EdgeIdentifier(const MERGE_INDEX ie) const 
    { return(ie); };
    MERGE_INDEX EdgeIdentifier
      (const MERGE_INDEX ie, const MERGE_INDEX j) const 
    { return(ie + j * num_edges); };

    // check function
    bool Check(IJK::ERROR & error) const;
    bool Check(const ISOTABLE_TYPE type, IJK::ERROR & error) const;
  };

  class IVOL_MERGE_DATA: public MERGE_DATA {
  public:
    IVOL_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 1, 2) {};
  };

  class ISO_MERGE_DATA: public MERGE_DATA {
  public:
    ISO_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 0, 1) {};
  };

  class NEP_ISO_MERGE_DATA: public MERGE_DATA {
  public:
    NEP_ISO_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 1, 1) {};
  };


// **************************************************
// NEP_ISOSURFACE_TABLE
// **************************************************

// NEP isosurface lookup table
// includes information about isosurface patches
// contained on polyhedron facets
  class NEP_ISOSURFACE_TABLE:public ISOSURFACE_TABLE {

  protected:
    bool * is_in_facet;   // is_in_facet[k] = true if iso patch is in facet
    IJKTABLE::FACET_INDEX * containing_facet;
                          // containing_facetr[k] = facet containing iso patch

    // initialization routine
    void Init(const int dimension, const int simplex_dimension);

  public:
    // constructors
    NEP_ISOSURFACE_TABLE(const int d);
    NEP_ISOSURFACE_TABLE(const int dimension, const int simplex_dimension);

    ~NEP_ISOSURFACE_TABLE();                // destructor

    // get functions
    bool IsInFacet(const TABLE_INDEX it) const
    { return(is_in_facet[it]); };
    IJKTABLE::FACET_INDEX ContainingFacet(const TABLE_INDEX it) const
      { return(containing_facet[it]); };

    // set functions

    // set is_in_facet[k] = flag, for all table entries k
    void SetIsInFacet(const bool flag); 

    void SetIsInFacet(const TABLE_INDEX it, const bool flag)
    { is_in_facet[it] = flag; };
    void SetContainingFacet(const TABLE_INDEX it, 
			    const IJKTABLE::FACET_INDEX jf)
    { is_in_facet[it] = true; containing_facet[it] = jf; };

    virtual void SetNumTableEntries(const int num_table_entries);

    // free memory
    virtual void FreeAll();                     // free all memory
  };

// **************************************************
// CLASS SNAP_MINMAX
// **************************************************

  /// MINMAX_REGIONS for snapMC
  class SNAP_MINMAX:public MINMAX_REGIONS {

  public:
    SNAP_MINMAX(){};

    // Compute min and max of each region
    // Offset each region by one cube to include all snapped cubes
    void ComputeMinMax
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const SCALAR_TYPE * scalar, const AXIS_SIZE_TYPE region_edge_length);

    void ComputeMinMax
    (const MC_GRID_BASE & scalar_grid, AXIS_SIZE_TYPE region_edge_length);
  };

// **************************************************
// CLASS SNAP_OCTREE
// **************************************************

  /// OCTREE for snapMC
  class SNAP_OCTREE:public IJKOCTREE::OCTREE {

  public:
    SNAP_OCTREE(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      IJKOCTREE::OCTREE(dimension, axis_size) {};

    // Offset each leaf node by one cube to include all snapped cubes
    void SetMinMax(const SCALAR_TYPE * scalar);
  };

};

#endif
