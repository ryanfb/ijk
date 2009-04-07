/// \file ijkgrid.txx
/// ijk templates defining regular grid classes
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008 Rephael Wenger

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

#ifndef _IJKGRID_
#define _IJKGRID_

#include "ijk.txx"

namespace IJKGRID {

  // **************************************************
  // TYPE DEFINITIONS
  // **************************************************

  typedef long GRID_SIZE_TYPE;

  // **************************************************
  // TEMPLATE CLASS GRID
  // **************************************************

  template <class DTYPE, class ATYPE, class VTYPE> class GRID {

  protected:
    DTYPE dimension;     // grid dimension
    ATYPE * axis_size;   // axis_size[i] = # grid points along axis i
    GRID_SIZE_TYPE num_vertices;  // number of grid vertices

    void Init(const DTYPE dimension, const ATYPE * axis_size);
    void Set(const DTYPE dimension, const ATYPE * axis_size);
    void FreeAll();

  public:
    GRID(const DTYPE dimension, const ATYPE * axis_size);
    GRID() { Init(0, NULL); };
    ~GRID();
    GRID(const GRID & grid);                      // copy constructor
    const GRID & operator = (const GRID & right); // copy assignment

    // get functions
    DTYPE Dimension() const { return(dimension); }
    const ATYPE * AxisSize() const { return(axis_size); }
    ATYPE AxisSize(const DTYPE i) const { return(axis_size[i]); }
    GRID_SIZE_TYPE NumVertices() const { return(num_vertices); }

    // compute functions
    GRID_SIZE_TYPE ComputeNumVertices
    (const VTYPE iv0, const VTYPE iv1) const;
    GRID_SIZE_TYPE ComputeNumCubes() const;
    template <class GTYPE>
    VTYPE ComputeVertexIndex(const GTYPE * coord) const;
    template <class GTYPE>
    void ComputeCoord(const VTYPE iv, GTYPE * coord) const;

    // compare
    bool CompareSize(const DTYPE dimension, const ATYPE * axis_size);

    // check function
    bool Check(const DTYPE dimension, const ATYPE * axis_size,
	       IJK::ERROR & error);
  };

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID_BASE
  // **************************************************

  /// Class storing scalar grid information
  /// Does not allocate or destroy any scalar information
  /// Note: SCALAR_GRID_BASE does not provide any way to allocate memory
  ///    for scalar or set scalar to point to existing memory.
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  class SCALAR_GRID_BASE:public GRID<DTYPE,ATYPE,VTYPE> {

  protected:
    STYPE * scalar;
    // point (x0,x1,x2,...) has scalar value
    //     scalar_grid[x0 + x1*axis_size[0] + 
    //                   x2*axis_size[0]*axis_size[1] + ...]

  public:
    SCALAR_GRID_BASE(const DTYPE dimension, const ATYPE * axis_size):
      GRID<DTYPE,ATYPE,VTYPE>(dimension,axis_size){ scalar = NULL; };
    // Note: constructor and destructor do not allocate or free scalar memory
    ~SCALAR_GRID_BASE() { scalar = NULL; };

    // copy constructor and assignment: NOT IMPLEMENTED
    SCALAR_GRID_BASE(const SCALAR_GRID_BASE & scalar_grid); 
    const SCALAR_GRID_BASE & operator = (const SCALAR_GRID_BASE & right);

    // set functions
    void Set(const VTYPE iv, const STYPE s) { scalar[iv] = s; };
    void Replace
    (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid2,
     const VTYPE iv0, const VTYPE iv1);

    // get functions
    STYPE * ScalarPtr() { return(scalar); };
    const STYPE * ScalarPtrConst() const { return(scalar); };
    const STYPE * End() const { return(scalar+this->NumVertices()); };
    STYPE Scalar(const VTYPE iv) const { return(scalar[iv]); };
  };

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID
  // **************************************************

  /// Class storing scalar grid information
  /// Allocates and destroys array scalar[]
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  class SCALAR_GRID:public SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> {

  protected:
    void FreeAll();

  public:
    SCALAR_GRID();
    SCALAR_GRID(const DTYPE dimension, const ATYPE * axis_size);
    ~SCALAR_GRID() { FreeAll(); };

    // copy constructor and assignment: NOT IMPLEMENTED
    SCALAR_GRID(const SCALAR_GRID & scalar_grid); 
    const SCALAR_GRID & operator = (const SCALAR_GRID & right);

    // set functions
    void SetSize(const DTYPE dimension, const ATYPE * axis_size);
    void Subsample
    (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
     const ATYPE resolution);
  };

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID_WRAPPER
  // **************************************************

  /// scalar grid wrapper for scalar array
  /// Does not allocate or destroy any scalar information
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  class SCALAR_GRID_WRAPPER:public SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> {

  public:
    SCALAR_GRID_WRAPPER(const DTYPE dimension, const ATYPE * axis_size,
			STYPE * scalar);
    ~SCALAR_GRID_WRAPPER() { this->scalar = NULL; };
    // Note: constructor and destructor do not allocate or free scalar memory
  };

  // **************************************************
  // TEMPLATE CLASS MINMAX_BASE
  // **************************************************

  // base class for representing min & max of scalar grid regions
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  class MINMAX_BASE:public GRID<DTYPE,ATYPE,VTYPE> {

  protected:
    STYPE * scalar_min;
    STYPE * scalar_max;

    void FreeAll();

  public:
    MINMAX_BASE();
    MINMAX_BASE(const DTYPE dimension, const ATYPE * axis_size);
    ~MINMAX_BASE(){ FreeAll(); };

    // copy constructor and assignment: NOT IMPLEMENTED
    MINMAX_BASE(const MINMAX_BASE & minmax_grid); 
    const MINMAX_BASE & operator = (const MINMAX_BASE & right);

    // set functions
    void SetSize(const DTYPE dimension, const ATYPE * axis_size);

    // copy min & max values
    void Copy(const STYPE * scalar_min, const STYPE * scalar_max);

    // get functions
    const STYPE * Min() const { return(scalar_min); };
    const STYPE * Max() const { return(scalar_max); };
    STYPE Min(const VTYPE iv) const { return(scalar_min[iv]); };
    STYPE Max(const VTYPE iv) const { return(scalar_max[iv]); };
  };

  // **************************************************
  // TEMPLATE CLASS MINMAX_GRID
  // **************************************************

  // class for computing min & max of scalar grid
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  class MINMAX_GRID:public MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE> {

  public:
    MINMAX_GRID(){};
    MINMAX_GRID(const DTYPE dimension, const ATYPE * axis_size):
      MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>(dimension, axis_size) {};
    ~MINMAX_GRID(){};

    // copy constructor and assignment: NOT IMPLEMENTED
    MINMAX_GRID(const MINMAX_GRID & minmax_grid); 
    const MINMAX_GRID & operator = (const MINMAX_GRID & right);

    // Compute min and max of each region and store in primary vertices
    void OverwriteWithMinMax(const ATYPE region_edge_length);

    void OverwriteWithMinMax(const ATYPE region_edge_length,
			     const ATYPE offset_edge_length);

    // Project scalar grid onto minmax grid
    // Compute min and max of projection
    void ProjectMinMax
    (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
     const ATYPE num_region_edges, const ATYPE facet_index, 
     const ATYPE axis_increment);

    void ProjectMinMax
    (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
     const ATYPE istart, const ATYPE iend, const ATYPE ionto,
     const ATYPE axis_increment);

    void ProjectMinMax
    (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
     const ATYPE num_region_edges, const ATYPE facet_index, 
     const ATYPE axis_increment);
  };

  // **************************************************
  // TEMPLATE CLASS MINMAX_REGIONS
  // **************************************************

  // class for computing min & max of scalar grid regions
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  class MINMAX_REGIONS:public MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE> {

  protected:
    ATYPE region_edge_length;

  public:
    MINMAX_REGIONS();
    MINMAX_REGIONS(const DTYPE dimension, const ATYPE * axis_size,
		   const ATYPE region_edge_length);
    ~MINMAX_REGIONS();

    // copy constructor and assignment: NOT IMPLEMENTED
    MINMAX_REGIONS(const MINMAX_REGIONS & minmax_regions); 
    const MINMAX_REGIONS & operator = (const MINMAX_REGIONS & right);

    // get functions
    ATYPE RegionEdgeLength() const { return(region_edge_length); };
    GRID_SIZE_TYPE NumRegions() const 
    { return(MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>::NumVertices()); };

    // remove function NumVertices. (Use NumRegions() instead.)
    GRID_SIZE_TYPE NumVertices() const;

    // Compute min and max of each region
    // Creates a grid of regions and stores min and max for each region
    void ComputeMinMax
    (const DTYPE dimension, const ATYPE * axis_size,
     const STYPE * scalar, const ATYPE region_edge_length);

    void ComputeMinMax
    (const DTYPE dimension, const ATYPE * axis_size,
     const STYPE * scalar, const ATYPE region_edge_length,
     const ATYPE offset_edge_length);

    void ComputeMinMax
    (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid, 
     const ATYPE region_edge_length);

    void ComputeMinMax
    (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid, 
     const ATYPE region_edge_length, const ATYPE offset_edge_length);
  };

  // **************************************************
  // TEMPLATE FUNCTIONS: COUNTING AND INDEXING
  // **************************************************

  /// Return coordinate of specified vertex
  template <class VTYPE, class DTYPE, class ATYPE, class GTYPE>
  void compute_coord(const VTYPE iv, const DTYPE dimension,
		     const ATYPE * axis_size, GTYPE * coord)
  // compute coordinates of grid vertex iv
  // Precondition: axis_size[d] > 0 for all d = 0,...,dimension-1
  {
    VTYPE k = iv;
    for (DTYPE d = 0; d < dimension; d++) {
      coord[d] = k % axis_size[d];
      k = k / axis_size[d];
    };
  }

  /// Return index of vertex with specified coord
  template <class VTYPE, class GTYPE, class DTYPE, class ATYPE>
  VTYPE compute_vertex_index
  (const GTYPE * coord, const DTYPE dimension, const ATYPE * axis_size)
  {
    VTYPE iv = 0;
    VTYPE inc = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      iv += inc*coord[d];
      inc = inc*axis_size[d];
    }

    return(iv);
  }

  /// Return number of grid vertices
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    GRID_SIZE_TYPE num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++)
      { num_vertices = num_vertices * axis_size[d]; }
    return(num_vertices);
  }

  /// Return number of grid vertices between two vertices
  template <class VTYPE, class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1)
  // get number of grid vertices between iv0 and iv1
  // Precondition: 0 <= iv0 <= iv1 < total_num_grid_vertices
  {
    VTYPE coord0[dimension];
    VTYPE coord1[dimension];
    IJK::PROCEDURE_ERROR error("compute_num_grid_vertices");

    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    GRID_SIZE_TYPE num_grid_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] > coord1[d]) {
	error.AddMessage("Programming error in calculating ", d,
			 "'th coordinate.");
	error.AddMessage("  coord0 (", coord0[d], 
			 ") > coord1 (", coord1[d], ").");
	throw error;
      }

      num_grid_vertices = num_grid_vertices * (coord1[d]-coord0[d]+1);
    }

    return(num_grid_vertices);
  }

  /// Return number of grid cubes
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    GRID_SIZE_TYPE num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] < 2) { num_cubes = 0; }
      else { num_cubes = num_cubes * (axis_size[d]-1); }
    }
    return(num_cubes);
  }

  /// Return number of cube vertices
  template <class DTYPE> 
  long compute_num_cube_vertices(const DTYPE dimension)
  { return(1L << dimension); };

  /// Return number of cube facet vertices
  template <class DTYPE> 
  long compute_num_cube_facet_vertices(const DTYPE dimension)
  // Precondition: dimension > 0
  { return(1L << dimension-1); };

  /// Return number of vertices in a region all of whose edges
  ///   have length region_edge_length
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_region_vertices
  (const DTYPE dimension, const ATYPE region_edge_length)
   // region_edge_length = length of each region edge measured in grid edges
  {
    GRID_SIZE_TYPE num_region_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      num_region_vertices = num_region_vertices * (region_edge_length+1);
    }
    return(num_region_vertices);
  }

  /// Compute number of regions along a single axis
  template <class ATYPE>
  ATYPE compute_num_regions_along_axis
  (const ATYPE axis_size, const ATYPE num_region_edges)
  // axis_size = axis size of given axis
  // num_region_edges = number of region edges along given axis
  {
    ATYPE num_regions_along_axis = 
      long(axis_size+num_region_edges-2)/long(num_region_edges);
    return(num_regions_along_axis);
  }

  /// Return total number of regions
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_regions
  (const DTYPE dimension, const ATYPE * axis_size, const ATYPE num_region_edges)
  {
    IJK::PROCEDURE_ERROR error("compute_num_regions");
    if (num_region_edges <= 0) {
      error.AddMessage("Programming error. Number of region edges must be positive.");
      throw error;
    }

    GRID_SIZE_TYPE num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 1) { return(0); };

      ATYPE num_regions_along_axis = 
	long(axis_size[d]+num_region_edges-2)/long(num_region_edges);
      num_regions = num_regions * num_regions_along_axis; 
    }
    return(num_regions);
  }

  /// Return total number of regions
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE * num_region_edges, const DTYPE orth_dir)
  // num_region_edges[d] = number of region edges in direction d
  // orth_dir = direction orthogonal to facet
  {
    GRID_SIZE_TYPE num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {

      if (d != orth_dir) {
	if (axis_size[d] <= 1) { return(0); };

	ATYPE num_regions_along_axis = 
	  compute_num_regions_along_axis(axis_size[d], num_region_edges[d]);
	num_regions = num_regions * num_regions_along_axis; 
      }
    }
    return(num_regions);
  }

  /// Compute 1D subsample size
  template <class ATYPE>
  GRID_SIZE_TYPE compute_subsample_size
  (const ATYPE axis_size, const ATYPE resolution)
  {
    GRID_SIZE_TYPE subsample_size =
      long(axis_size+resolution-1)/long(resolution);
    return(subsample_size);
  }

  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_subsample_size
  (const DTYPE dimension, const ATYPE * axis_size, const ATYPE * resolution)
  // resolution[d] = resolution along axis d
  {
    IJK::PROCEDURE_ERROR error("compute_subsample_size");

    GRID_SIZE_TYPE num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] < 1) { return(0); };
      if (resolution[d] < 1) {
	error.AddMessage("Subsample resolution must be positive.");
	throw error;
      }

      ATYPE num_subsampled_along_axis =
	compute_subsample_size(axis_size[d], resolution[d]);
      num_vertices = num_vertices*num_subsampled_along_axis;
    }
    return(num_vertices);
  }

  /// Check range of vertices
  template <class VTYPE, class DTYPE, class ATYPE>
  bool check_range
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1, IJK::ERROR & error)
  // return true if 0 <= iv0 <= iv1 < total_num_grid_vertices
  {
    const VTYPE total_num_grid_vertices =
      compute_num_grid_vertices<VTYPE>(dimension, axis_size);

    if (iv0 > iv1) {
      error.AddMessage("Illegal vertex range. Vertex index ", iv0, 
		       " is greater than vertex index ", iv1, ".");
      return(false);
    }

    if (0 > iv0 || iv1 >= total_num_grid_vertices) {
      error.AddMessage("Illegal vertex indices: ", iv0, ", ", iv1, ".");
      error.AddMessage("Vertex indices should be in range: [0,",
		       total_num_grid_vertices, "].");
      return(false);
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTING INCREMENTS
  // **************************************************

  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_increment
  (const DTYPE dimension, const ATYPE * axis_size, ITYPE * increment)
  // compute increment to add to index of current vertex to get
  //   next vertex along each axis
  // Precondition: array increment is allocated with size at least dimension
  {
    IJK::PROCEDURE_ERROR error("compute_increment");

    if (dimension <= 0) { return; };

    if (axis_size == NULL || increment == NULL) {
      error.AddMessage("Programming error. axis_size == NULL or increment == NULL.");
      throw error;
    }

    increment[0] = 1;
    for (DTYPE d = 1; d < dimension; d++)
      { increment[d] = increment[d-1]*axis_size[d-1]; }
  }

  template <class DTYPE, class ATYPE, class VTYPE, class ITYPE>
  void compute_increment
  (const GRID<DTYPE,ATYPE,VTYPE> & grid, ITYPE * increment)
  // compute increment to add to index of current vertex to get
  //   next vertex along each axis
  // Precondition: array increment is allocated with size at least dimension
  {
    compute_increment(grid.Dimension(), grid.AxisSize(), increment);
  }

  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_region_vertex_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, ITYPE * increment)
  // compute increment to add to index of primary vertex to get
  //   k'th vertex in region
  // Precondition: array increment is allocated with size at least 
  //   number of region vertices
  {
    const DTYPE d_last = dimension-1;
    ITYPE axis_increment[dimension];
    IJK::PROCEDURE_ERROR error("compute_region_vertex_increment");

    if (dimension <= 0) { return; };

    if (axis_size == NULL || increment == NULL) {
      error.AddMessage("Programming error. axis_size == NULL or increment == NULL.");
      throw error;
    }

    compute_increment(dimension, axis_size, axis_increment);

    const GRID_SIZE_TYPE num_region_vertices = 
      compute_num_region_vertices(dimension, region_edge_length);

    // initialize iv
    increment[0] = 0;
    ITYPE prev_num_vert = 1;

    // Process axes 0,1,2,..., dimension-1
    ITYPE * vcur_ptr = increment+prev_num_vert;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE v_inc = axis_increment[d];
      for (ATYPE j = 1; j <= region_edge_length; j++) {
	for (ITYPE * vprev_ptr = increment; 
	     vprev_ptr != increment+prev_num_vert; vprev_ptr++) {
	  *(vcur_ptr) = v_inc + *(vprev_ptr);
	  vcur_ptr++;
	};
	v_inc = v_inc + axis_increment[d];
      }
      prev_num_vert = prev_num_vert*(region_edge_length+1);
    }

    if (prev_num_vert != num_region_vertices) {
      error.AddMessage("Programming error.  Computed increment for ", 
		       prev_num_vert, " vertices.");
      error.AddMessage("Number of vertices in region is ", 
		       num_region_vertices, ".");
      throw error;
    }
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: GETTING VERTICES
  // **************************************************

  // WARNING: Procedure get_grid_vertices() has not been tested.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const VTYPE iv1, VTYPE * vlist)
  // get list of grid vertices between iv0 and iv1
  // Precondition: 0 <= iv0 <= iv1 < total_num_grid_vertices
  // Precondition: vlist[] is preallocated to size at least num_grid_vertices
  {
    VTYPE coord0[dimension];
    VTYPE coord1[dimension];
    VTYPE coord[dimension];
    VTYPE axis_increment[dimension];
    IJK::PROCEDURE_ERROR error("get_grid_vertices");

    if (dimension <= 0) { return; };

    if (axis_size == NULL) {
      error.AddMessage("Programming error. NULL axis size.");
      throw error;
    }

    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    VTYPE num_grid_vertices = 
      compute_num_grid_vertices(dimension, axis_size, iv0, iv1);

    if (num_grid_vertices < 1) return;

    if (vlist == NULL) {
      error.AddMessage("Programming error. NULL vlist.");
      throw error;
    }

    compute_increment(dimension, axis_size, axis_increment);
    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    // initialize prev_num_cubes
    vlist[0] = iv0;
    VTYPE prev_num_vert = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vert;
    for (DTYPE d = 0; d < dimension; d++) {
      VTYPE v_inc = axis_increment[d];
      for (VTYPE c = coord0[d]+1; c <= coord1[d]; c++) {
	for (VTYPE * vprev_ptr = vlist; 
	     vprev_ptr != vlist+prev_num_vert; vprev_ptr++) {
	  *(vcur_ptr) = v_inc + *(vprev_ptr);
	  vcur_ptr++;
	};
	v_inc = v_inc + axis_increment[d];
      }
      prev_num_vert = prev_num_vert*(coord1[d]-coord0[d]+1);
    }

    if (prev_num_vert != num_grid_vertices) {
      error.AddMessage("Programming error.  Added ", prev_num_vert, 
		       " grid vertices to list.");
      error.AddMessage("Number of grid vertices between vertex ", iv0,
		       " and vertex ", iv1, " is ", num_grid_vertices, ".");
      throw error;
    }
  }

  // WARNING: Procedure get_primary_cube_vertices() has not been tested.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_primary_cube_vertices
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  // get list of primary cube vertices (lowest coord vertices in cubes)
  // vlist[k] = list of primary cube vertices
  //    Precondition: vlist[] is preallocated to length at least 
  //                  number of cubes
  {
    IJK::PROCEDURE_ERROR error("get_primary_cube_vertices");
    VTYPE axis_increment[dimension];

    if (dimension <= 0) { return; };

    if (axis_size == NULL) {
      error.AddMessage("Programming error. NULL axis size.");
      throw error;
    }

    const VTYPE num_cubes = 
      compute_num_grid_cubes(dimension, axis_size);
    if (num_cubes < 1) return;
    // Note: axis_size[d] >= 1 for all d. dimension >= 1.

    if (vlist == NULL) {
      error.AddMessage("Programming error. NULL vlist.");
      throw error;
    }

    compute_increment(dimension, axis_size, axis_increment);

    // initialize prev_num_cubes
    VTYPE prev_num_cubes = 0;

    // Process first axis separately for faster execution time
    // (Makes little difference unless first axis is very long)
    if (dimension > 0) {
      if (axis_size[0] < 1) {
	error.AddMessage("Programming error. axis_size[0] < 1.");
	throw error;
      }
      if (axis_increment[0] != 1) {
	error.AddMessage("Programming error. axis_increment[0] != 1.");
	throw error;
      }

      for (VTYPE i = 0; i < axis_size[0]-1; i++) {
	vlist[i] = i;
      }
      prev_num_cubes = axis_size[0]-1;
    }

    // Process axes 1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_cubes;
    for (int d = 1; d < dimension; d++) {
      if (axis_increment[d] < 1) {
	error.AddMessage("Programming error. axis_increment[", d,
			 "] < 1.");
	throw error;
      }

      VTYPE iv0 = axis_increment[d];
      for (VTYPE i = 1; i < axis_size[d]-1; i++) {
	for (VTYPE * vprev_ptr = vlist; 
	     vprev_ptr != vlist+prev_num_cubes; vprev_ptr++) {
	  *(vcur_ptr) = iv0 + *(vprev_ptr);
	  vcur_ptr++;
	};
	iv0 = iv0 + axis_increment[d];
      }
      prev_num_cubes = prev_num_cubes*(axis_size[d]-1);
    }

    if (prev_num_cubes != num_cubes) {
      error.AddMessage("Programming error.  Added ", prev_num_cubes, 
		       " primary cube vertices to list.");
      error.AddMessage("Number of primary cube vertices should be ", 
		       num_cubes, ".");
      throw error;
    }
  }


  template <class DTYPE, class ATYPE, class VTYPE>
  void subsample_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const ATYPE * resolution,
   VTYPE * vlist)
  // resolution[d] = subsample resolution along axis d
  // Precondition: vlist[] is preallocated to size at least num_grid_vertices
  {
    IJK::PROCEDURE_ERROR error("subsample_grid_vertices");

    const GRID_SIZE_TYPE subsample_size =
      compute_subsample_size(dimension, axis_size, resolution);

    if (subsample_size < 1) { return; };
    // Note: axis_size[d] >= 1 for all d != facet_index

    if (vlist == NULL) {
      error.AddMessage("Programming error. NULL vlist.");
      throw error;
    }

    VTYPE axis_increment[dimension];
    compute_increment(dimension, axis_size, axis_increment);

    // initialize prev_num_vertices
    vlist[0] = 0;
    VTYPE prev_num_vertices = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 1) {
	error.AddMessage("Programming error. axis_size[", d, "] <= 1.");
	throw error;
      }

      ATYPE num_subsampled_along_axis =
	(axis_size[d]+resolution[d]-1)/resolution[d];

      VTYPE iv0 = resolution[d]*axis_increment[d];
      for (VTYPE i = 1; i < num_subsampled_along_axis; i++) {
	for (VTYPE * vprev_ptr = vlist; 
	     vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
	  *(vcur_ptr) = iv0 + *(vprev_ptr);
	  vcur_ptr++;
	};
	iv0 = iv0 + resolution[d]*axis_increment[d];
      }
      prev_num_vertices = prev_num_vertices*(num_subsampled_along_axis);
    }

    if (prev_num_vertices != subsample_size) {
      error.AddMessage("Programming error.  Added ", prev_num_vertices, 
		       " to subsample list.");
      error.AddMessage("Number of subsampled vertices should be ", 
		       subsample_size, ".");
      throw error;
    }

  }
  
  // ********************************************************
  // TEMPLATE FUNCTIONS: FACET VERTICES, CUBES AND REGIONS
  // ********************************************************

  /// Return number of vertices in specified grid facet
  /// Specify grid facet by the direction orthogonal to the facet.
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_vertices_in_grid_facet
    (const DTYPE dimension, const ATYPE * axis_size,
     const DTYPE orth_dir)
  {
    GRID_SIZE_TYPE num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {
	num_vertices = num_vertices * axis_size[d];
      };
    }
    return(num_vertices);
  }

  /// Return number of vertices in specified grid ridge
  /// Specify grid facet by the directions orthogonal to ridge
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_vertices_in_grid_ridge
    (const DTYPE dimension, const ATYPE * axis_size,
     const DTYPE orth_dir0, const DTYPE orth_dir1)
  {
    GRID_SIZE_TYPE num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir0 && d != orth_dir1) {
	num_vertices = num_vertices * axis_size[d];
      };
    }
    return(num_vertices);
  }

  /// Return maximum number of vertices over all grid facets
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_max_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    GRID_SIZE_TYPE max_num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      GRID_SIZE_TYPE num_face_vertices = 
	compute_num_vertices_in_grid_facet(dimension, axis_size, d);
      if (num_face_vertices > max_num_vertices)
	max_num_vertices = num_face_vertices;
    };

    return(max_num_vertices);
  }

  /// Return number of cubes in specified grid facet
  /// Facet is lower facet orthogonal to the specified direction
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_cubes_in_grid_facet
    (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir)
  {
    GRID_SIZE_TYPE num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 1) { return(0); };
      if (d != orth_dir) {
	num_cubes = num_cubes*(axis_size[d]-1);
      };
    }
    return(num_cubes);
  }

  /// Return number of cubes in grid facet orthogonal to axis 0
  template <class DTYPE, class ATYPE>
  GRID_SIZE_TYPE compute_num_cubes_in_grid_facet0
    (const DTYPE dimension, const ATYPE * axis_size)
  {
    return(compute_num_cubes_in_grid_facet(dimension, axis_size, 0));
  }

  /// Get vertices in specified grid facet
  /// Facet is lower facet orthogonal to the specified direction
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_vertices_in_grid_facet
    (const DTYPE dimension, const ATYPE * axis_size,
     const DTYPE orth_dir, VTYPE * vlist)
  // vlist[] = list of vertices in grid_face
  //    Precondition: vlist[] is preallocated to length at least 
  //                  number of vertices in grid facet
  {
    IJK::PROCEDURE_ERROR error("get_vertices_in_grid_facet");

    const GRID_SIZE_TYPE num_facet_vertices =
      compute_num_vertices_in_grid_facet(dimension, axis_size, orth_dir);

    if (num_facet_vertices < 1) { return; };
    // Note: axis_size[d] >= 1 for all d != facet_index

    if (vlist == NULL) {
      error.AddMessage("Programming error. NULL vlist.");
      throw error;
    }

    VTYPE axis_increment[dimension];
    compute_increment(dimension, axis_size, axis_increment);

    // initialize prev_num_vertices
    vlist[0] = 0;
    VTYPE prev_num_vertices = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {
	if (axis_size[d] < 1) {
	  error.AddMessage("Programming error. axis_size[", d, "] < 1.");
	  throw error;
	}

	VTYPE iv0 = axis_increment[d];
	for (VTYPE i = 1; i < axis_size[d]; i++) {
	  for (VTYPE * vprev_ptr = vlist; 
	       vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
	    *(vcur_ptr) = iv0 + *(vprev_ptr);
	    vcur_ptr++;
	  };
	  iv0 = iv0 + axis_increment[d];
	}
	prev_num_vertices = prev_num_vertices*(axis_size[d]);
      }
    }

    if (prev_num_vertices != num_facet_vertices) {
      error.AddMessage("Programming error.  Added ", prev_num_vertices, 
		       " vertices in grid face to list.");
      error.AddMessage("Number of vertices in grid face should be ", 
		       num_facet_vertices, ".");
      throw error;
    }
  }

  /// Get vertices in specified grid ridge
  /// Ridge is lower ridge orthogonal to the specified directions
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_vertices_in_grid_ridge
    (const DTYPE dimension, const ATYPE * axis_size,
     const DTYPE orth_dir0, const DTYPE orth_dir1, VTYPE * vlist)
  // vlist[] = list of vertices in grid ridge
  //    Precondition: vlist[] is preallocated to length at least 
  //                  number of vertices in grid ridge
  {
    IJK::PROCEDURE_ERROR error("get_vertices_in_grid_ridge");

    const GRID_SIZE_TYPE num_ridge_vertices =
      compute_num_vertices_in_grid_ridge
      (dimension, axis_size, orth_dir0, orth_dir1);

    if (num_ridge_vertices < 1) { return; };
    // Note: axis_size[d] >= 1 for all d != orthogonal directions

    if (vlist == NULL) {
      error.AddMessage("Programming error. NULL vlist.");
      throw error;
    }

    VTYPE axis_increment[dimension];
    compute_increment(dimension, axis_size, axis_increment);

    // initialize prev_num_vertices
    vlist[0] = 0;
    VTYPE prev_num_vertices = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir0 && d != orth_dir1) {
	if (axis_size[d] < 1) {
	  error.AddMessage("Programming error. axis_size[", d, "] < 1.");
	  throw error;
	}

	VTYPE iv0 = axis_increment[d];
	for (VTYPE i = 1; i < axis_size[d]; i++) {
	  for (VTYPE * vprev_ptr = vlist; 
	       vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
	    *(vcur_ptr) = iv0 + *(vprev_ptr);
	    vcur_ptr++;
	  };
	  iv0 = iv0 + axis_increment[d];
	}
	prev_num_vertices = prev_num_vertices*(axis_size[d]);
      }
    }

    if (prev_num_vertices != num_ridge_vertices) {
      error.AddMessage("Programming error.  Added ", prev_num_vertices, 
		       " vertices in grid ridge to list.");
      error.AddMessage("Number of vertices in grid ridge should be ", 
		       num_ridge_vertices, ".");
      throw error;
    }
  }

  /// Get primary vertices of cubes in specified grid facet
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, const bool side, VTYPE * vlist)
  {
    IJK::PROCEDURE_ERROR error("get_cubes_in_grid_facet");

    GRID_SIZE_TYPE num_cubes =
      compute_num_cubes_in_grid_facet(dimension, axis_size, orth_dir);
    if (num_cubes < 1) return;
    // Note: axis_size[d] >= 1 for all d >= 1. dimension >= 1.

    if (vlist == NULL) {
      error.AddMessage("Programming error. NULL vlist.");
      throw error;
    }

    VTYPE axis_increment[dimension];
    compute_increment(dimension, axis_size, axis_increment);

    // initialize prev_num_cubes
    if (!side)
      { vlist[0] = 0; }
    else 
      { vlist[0] = axis_increment[orth_dir]*(axis_size[orth_dir]-1); };
    VTYPE prev_num_cubes = 1;

    // Process axes 0,1,2,..., dimension-1 (except for orth_dir)
    VTYPE * vcur_ptr = vlist + prev_num_cubes;
    for (DTYPE d = 0; d < dimension; d++) {

      if (d != orth_dir) {
	if (axis_size[d] < 1) {
	  error.AddMessage("Programming error. axis_size[", d, "] < 1.");
	  throw error;
	}

	VTYPE v_inc = axis_increment[d];
	for (VTYPE i = 1; i < axis_size[d]-1; i++) {
	  for (VTYPE * vprev_ptr = vlist; 
	       vprev_ptr != vlist+prev_num_cubes; vprev_ptr++) {
	    *(vcur_ptr) = v_inc + *(vprev_ptr);
	    vcur_ptr++;
	  };
	  v_inc = v_inc + axis_increment[d];
	}

	prev_num_cubes = prev_num_cubes*(axis_size[d]-1);
      }
    }

    if (prev_num_cubes != num_cubes) {
      error.AddMessage("Programming error.  Added ", prev_num_cubes, 
		       " primary cube vertices to list.");
      error.AddMessage("Number of primary cube vertices in face 0 should be ", 
		       num_cubes, ".");
      throw error;
    }

  }

  template <class DTYPE, class ATYPE, class VTYPE>
  void get_cubes_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_cubes_in_grid_facet(dimension, axis_size, 0, false, vlist);
  }

  /// Get regions of specified grid facet
  /// Facet is lower facet orthogonal to the specified direction
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_regions_in_grid_facet
    (const DTYPE dimension, const ATYPE * axis_size,
     const ATYPE * num_region_edges, const DTYPE orth_dir, VTYPE * vlist)
  // vlist[] = list of primary vertices of regions in grid_face
  //    Precondition: vlist[] is preallocated to length at least 
  //                  number of regions in grid facet
  {
    IJK::PROCEDURE_ERROR error("get_regions_in_grid_facet");

    const GRID_SIZE_TYPE num_regions =
      compute_num_regions_in_grid_facet
      (dimension, axis_size, num_region_edges, orth_dir);

    if (num_regions < 1) { return; };
    // Note: axis_size[d] >= 1 for all d != facet_index

    if (vlist == NULL) {
      error.AddMessage("Programming error. NULL vlist.");
      throw error;
    }

    VTYPE axis_increment[dimension];
    compute_increment(dimension, axis_size, axis_increment);

    // initialize prev_num_vertices
    vlist[0] = 0;
    VTYPE prev_num_vertices = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {
	if (axis_size[d] <= 1) {
	  error.AddMessage("Programming error. axis_size[", d, "] <= 1.");
	  throw error;
	}

	ATYPE num_regions_along_axis =
	  compute_num_regions_along_axis(axis_size[d], num_region_edges[d]);

	VTYPE iv0 = num_region_edges[d]*axis_increment[d];
	for (VTYPE i = 1; i < num_regions_along_axis; i++) {
	  for (VTYPE * vprev_ptr = vlist; 
	       vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
	    *(vcur_ptr) = iv0 + *(vprev_ptr);
	    vcur_ptr++;
	  };
	  iv0 = iv0 + num_region_edges[d]*axis_increment[d];
	}
	prev_num_vertices = prev_num_vertices*(num_regions_along_axis);
      }
    }

    if (prev_num_vertices != num_regions) {
      error.AddMessage("Programming error.  Added ", prev_num_vertices, 
		       " regions in grid face to list.");
      error.AddMessage("Number of regions in grid face should be ", 
		       num_regions, ".");
      throw error;
    }
  }

  
  // ********************************************************
  // TEMPLATE FUNCTIONS: COMPUTE CUBE/REGION MIN AND MAX
  // ********************************************************

  /// Compute minimum scalar value of each cube
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void compute_cube_min
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
   STYPE * cube_min)
  // cube_min[iv] = min scalar value of cube with primary vertex iv
  //   Precondition: cube_min[] is preallocated to size 
  //                 at least num_grid_vertices
  {
    const DTYPE dimension = scalar_grid.Dimension();
    const ATYPE * axis_size = scalar_grid.AxisSize();
    VTYPE axis_increment[dimension];
    VTYPE * facet_vlist = NULL;

    const GRID_SIZE_TYPE num_grid_vertices = 
      compute_num_grid_vertices(dimension, axis_size);
    if (num_grid_vertices == 0) { return; };

    // copy scalar_grid to cube_min
    std::copy(scalar_grid.ScalarPtrConst(), scalar_grid.End(),
	      cube_min);

    // compute axis_increment
    compute_increment(scalar_grid, axis_increment);

    try {
      const GRID_SIZE_TYPE max_num_grid_facet_vertices =
	compute_max_num_grid_facet_vertices(dimension, axis_size);

      facet_vlist = new VTYPE[max_num_grid_facet_vertices];

      for (DTYPE d = 0; d < dimension; d++) {

	const GRID_SIZE_TYPE num_grid_facet_vertices =
	  compute_num_vertices_in_grid_facet(dimension, axis_size, d);

	get_vertices_in_grid_facet(dimension, axis_size, d, facet_vlist);

	for (VTYPE iv = 0; iv < num_grid_facet_vertices; iv++) {

	  VTYPE v0 = facet_vlist[iv];
	  STYPE min0 = cube_min[v0];

	  for (ATYPE k = 1; k < axis_size[d]; k++) {

	    VTYPE v1 = v0 + axis_increment[d];
	    STYPE min1 = cube_min[v1];
	    if (min0 > min1) { cube_min[v0] = min1; }

	    v0 = v1;
	    min0 = min1;
	  }
	}
      }
    }
    catch(...) {
      delete [] facet_vlist;
      facet_vlist = NULL;
      throw;
    }

    delete [] facet_vlist;
    facet_vlist = NULL;
  }

  /// Compute maximum scalar value of each cube
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void compute_cube_max
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
   STYPE * cube_max)
  // cube_max[iv] = max scalar value of cube with primary vertex iv
  //   Precondition: cube_max[] is preallocated to size 
  //                 at least num_grid_vertices
  {
    const DTYPE dimension = scalar_grid.Dimension();
    const ATYPE * axis_size = scalar_grid.AxisSize();
    VTYPE axis_increment[dimension];
    VTYPE * facet_vlist = NULL;

    const GRID_SIZE_TYPE num_grid_vertices = 
      compute_num_grid_vertices(dimension, axis_size);
    if (num_grid_vertices == 0) { return; };

    // copy scalar_grid to cube_min
    std::copy(scalar_grid.ScalarPtrConst(), scalar_grid.End(),
	      cube_max);

    // compute axis_increment
    compute_increment(scalar_grid, axis_increment);

    try {
      const GRID_SIZE_TYPE max_num_grid_facet_vertices =
	compute_max_num_grid_facet_vertices(dimension, axis_size);

      facet_vlist = new VTYPE[max_num_grid_facet_vertices];

      for (DTYPE d = 0; d < dimension; d++) {

	const GRID_SIZE_TYPE num_grid_facet_vertices =
	  compute_num_vertices_in_grid_facet(dimension, axis_size, d);

	get_vertices_in_grid_facet(dimension, axis_size, d, facet_vlist);

	for (VTYPE iv = 0; iv < num_grid_facet_vertices; iv++) {

	  VTYPE v0 = facet_vlist[iv];
	  STYPE max0 = cube_max[v0];

	  for (ATYPE k = 1; k < axis_size[d]; k++) {

	    VTYPE v1 = v0 + axis_increment[d];
	    STYPE max1 = cube_max[v1];
	    if (max0 < max1) { cube_max[v0] = max1; }

	    v0 = v1;
	    max0 = max1;
	  }
	}
      }
    }
    catch(...) {
      delete [] facet_vlist;
      facet_vlist = NULL;
      throw;
    }

    delete [] facet_vlist;
    facet_vlist = NULL;
  }

  /// Compute minimum and maximum scalar value of each region
  template <class DTYPE, class ATYPE, class STYPE>
  void compute_region_minmax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE num_region_edges, STYPE * region_min, STYPE * region_max)
  // num_region_edges = number of grid edges per region edge
  // region_min[j] = min scalar value of region j
  // region_max[j] = max scalar value of region j
  //   Precondition: region_min[] and region_max[] are preallocated to size 
  //                 at least number of regions
  {
    IJK::PROCEDURE_ERROR error("compute_region_minmax");
    const GRID_SIZE_TYPE num_vertices = 
      compute_num_grid_vertices(dimension, axis_size);
    GRID_SIZE_TYPE axis_increment[dimension];
    GRID_SIZE_TYPE vprimary[dimension];
    std::vector<STYPE> facet_min;
    std::vector<STYPE> facet_max;

    if (num_region_edges <= 0) {
      error("Programming error.  Number of region edges must be positive.");
      throw error;
    }

    const GRID_SIZE_TYPE num_regions = 
      compute_num_regions(dimension, axis_size, num_region_edges);

    if (num_regions < 1 || num_vertices < 1) { return; };

    region_min[0] = scalar[0];
    region_max[0] = scalar[0];
    
    if (dimension < 1) { return;  };

    MINMAX_GRID<DTYPE,ATYPE,GRID_SIZE_TYPE,STYPE> 
      minmax_grid(dimension-1, axis_size);

    compute_increment(dimension, axis_size, axis_increment);

    const DTYPE d_last = dimension-1;
    const GRID_SIZE_TYPE num_facet_vertices =
      compute_num_vertices_in_grid_facet(dimension, axis_size, d_last);

    ATYPE num_regions_along_axis = 
      compute_num_regions_along_axis(axis_size[d_last], num_region_edges);
    GRID_SIZE_TYPE num_regions_in_facet =
      compute_num_regions(dimension-1, axis_size, num_region_edges);

    GRID_SIZE_TYPE iregion = 0;
    if (dimension > 1) {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, num_region_edges, j*num_region_edges,
	   axis_increment[d_last]);

	minmax_grid.OverwriteWithMinMax(num_region_edges);

	GRID_SIZE_TYPE facet_increment = (j*num_region_edges)*axis_increment[d_last];
	for (DTYPE d = 0; d < dimension; d++) 
	  { vprimary[d] = facet_increment; }

	for (int jregion = 0; jregion < num_regions_in_facet; jregion++) {

	  GRID_SIZE_TYPE iv = vprimary[0]-facet_increment;
	  region_min[iregion] = minmax_grid.Min(iv);
	  region_max[iregion] = minmax_grid.Max(iv);

	  iregion++;

	  DTYPE d = 0;
	  while (d+1 < dimension && 
		 (vprimary[d] + (num_region_edges+1)*axis_increment[d] >=
		  vprimary[d+1] + axis_increment[d+1])) {
	    d++;
	  }

	  // go to next subspace of dimension d
	  vprimary[d] = vprimary[d] + num_region_edges*axis_increment[d];

	  for (DTYPE d2 = 0; d2 < d; d2++) 
	    { vprimary[d2] = vprimary[d]; };
	}
      }
    }
    else {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, num_region_edges, j*num_region_edges,
	   axis_increment[d_last]);

	region_min[iregion] = minmax_grid.Min(0);
	region_max[iregion] = minmax_grid.Max(0);
	iregion++;
      }
    }

    if (iregion != num_regions) {
      error.AddMessage("Programming error.  Computed min value for ", iregion,
		       " regions.");
      error.AddMessage("Number of regions is ", num_regions, ".");
      throw error;
    }
  }

  /// Compute minimum and maximum scalar value of each region
  template <class DTYPE, class ATYPE, class STYPE>
  void compute_region_minmax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE num_region_edges, const ATYPE offset_edge_length,
   STYPE * region_min, STYPE * region_max)
  // num_region_edges = number of grid edges per region edge
  // offset_edge_length = offset edge length
  // region_min[j] = min scalar value of region j
  // region_max[j] = max scalar value of region j
  //   Precondition: region_min[] and region_max[] are preallocated to size 
  //                 at least number of regions
  {
    IJK::PROCEDURE_ERROR error("compute_region_minmax");
    const GRID_SIZE_TYPE num_vertices = 
      compute_num_grid_vertices(dimension, axis_size);
    GRID_SIZE_TYPE axis_increment[dimension];
    GRID_SIZE_TYPE vprimary[dimension];
    std::vector<STYPE> facet_min;
    std::vector<STYPE> facet_max;

    if (num_region_edges <= 0) {
      error("Programming error.  Number of region edges must be positive.");
      throw error;
    }

    const GRID_SIZE_TYPE num_regions = 
      compute_num_regions(dimension, axis_size, num_region_edges);

    if (num_regions < 1 || num_vertices < 1) { return; };

    region_min[0] = scalar[0];
    region_max[0] = scalar[0];
    
    if (dimension < 1) { return;  };

    MINMAX_GRID<DTYPE,ATYPE,GRID_SIZE_TYPE,STYPE> 
      minmax_grid(dimension-1, axis_size);

    compute_increment(dimension, axis_size, axis_increment);

    const DTYPE d_last = dimension-1;
    const GRID_SIZE_TYPE num_facet_vertices =
      compute_num_vertices_in_grid_facet(dimension, axis_size, d_last);

    ATYPE num_regions_along_axis = 
      compute_num_regions_along_axis(axis_size[d_last], num_region_edges);
    GRID_SIZE_TYPE num_regions_in_facet =
      compute_num_regions(dimension-1, axis_size, num_region_edges);

    GRID_SIZE_TYPE iregion = 0;
    if (dimension > 1) {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	ATYPE istart = 0;
	if (j > 0) { istart = j*num_region_edges-offset_edge_length; };
	ATYPE iend = (j+1)*num_region_edges+offset_edge_length+1;
	if (iend > axis_size[d_last]) { iend = axis_size[d_last]; };

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, istart, iend, j*num_region_edges,
	   axis_increment[d_last]);

	minmax_grid.OverwriteWithMinMax(num_region_edges, offset_edge_length);

	GRID_SIZE_TYPE facet_increment = (j*num_region_edges)*axis_increment[d_last];
	for (DTYPE d = 0; d < dimension; d++) 
	  { vprimary[d] = facet_increment; }

	for (int jregion = 0; jregion < num_regions_in_facet; jregion++) {

	  GRID_SIZE_TYPE iv = vprimary[0]-facet_increment;
	  region_min[iregion] = minmax_grid.Min(iv);
	  region_max[iregion] = minmax_grid.Max(iv);

	  iregion++;

	  DTYPE d = 0;
	  while (d+1 < dimension && 
		 (vprimary[d] + (num_region_edges+1)*axis_increment[d] >=
		  vprimary[d+1] + axis_increment[d+1])) {
	    d++;
	  }

	  // go to next subspace of dimension d
	  vprimary[d] = vprimary[d] + num_region_edges*axis_increment[d];

	  for (DTYPE d2 = 0; d2 < d; d2++) 
	    { vprimary[d2] = vprimary[d]; };
	}
      }
    }
    else {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	ATYPE istart = 0;
	if (j > 0) { istart = j*num_region_edges-offset_edge_length; };
	ATYPE iend = (j+1)*num_region_edges+offset_edge_length+1;
	if (iend > axis_size[d_last]) { iend = axis_size[d_last]; };

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, istart, iend, j*num_region_edges,
	   axis_increment[d_last]);

	region_min[iregion] = minmax_grid.Min(0);
	region_max[iregion] = minmax_grid.Max(0);
	iregion++;
      }
    }

    if (iregion != num_regions) {
      error.AddMessage("Programming error.  Computed min value for ", iregion,
		       " regions.");
      error.AddMessage("Number of regions is ", num_regions, ".");
      throw error;
    }
  }


  // **************************************************
  // TEMPLATE CLASS GRID MEMBER FUNCTIONS
  // **************************************************

  template <class DTYPE, class ATYPE, class VTYPE> 
  GRID<DTYPE,ATYPE,VTYPE>::GRID
  (const DTYPE dimension, const ATYPE * axis_size)
  // constructor
  {
    Init(dimension, axis_size);
  }

  template <class DTYPE, class ATYPE, class VTYPE> 
  GRID<DTYPE,ATYPE,VTYPE>::~GRID<DTYPE,ATYPE,VTYPE>()
  // destructor
  {
    FreeAll();
  }

  template <class DTYPE, class ATYPE, class VTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE>::Init
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    this->axis_size = NULL;
    this->dimension = 0;
    this->num_vertices = 1;
    if (dimension > 0) 
      { Set(dimension, axis_size); };
  }

  template <class DTYPE, class ATYPE, class VTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE>::Set
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    FreeAll();
    this->dimension = dimension;
    this->axis_size = new ATYPE[dimension];
    for (DTYPE d = 0; d < dimension; d++)
      { this->axis_size[d] = axis_size[d]; }
    this->num_vertices = 
      compute_num_grid_vertices<VTYPE>(dimension, axis_size);

    if (dimension < 0) {
      IJK::PROCEDURE_ERROR error("Grid::Set");
      error.AddMessage("Programming error.  Illegal dimension ",
		       dimension, ".");
      error.AddMessage("Dimension should be non-negative.");
      throw error;
    }
  }

  template <class DTYPE, class ATYPE, class VTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE>::FreeAll()
  {
    if (axis_size == NULL) { delete [] axis_size; };
    axis_size = NULL;
    dimension = 0;
    num_vertices = 0;
  }

  template <class DTYPE, class ATYPE, class VTYPE> 
  GRID<DTYPE,ATYPE,VTYPE>::
  GRID(const GRID<DTYPE,ATYPE,VTYPE> & grid)
  // copy constructor
  {
    Init(grid.Dimension(), grid.AxisSize());
  }

  template <class DTYPE, class ATYPE, class VTYPE> 
  const GRID<DTYPE,ATYPE,VTYPE> & 
  GRID<DTYPE,ATYPE,VTYPE>::operator = (const GRID<DTYPE,ATYPE,VTYPE> & right)
  // copy assignment	
  {
    if (&right != this) {         // avoid self-assignment
      Set(right.Dimension(), right.AxisSize());
    }
  }

  template <class DTYPE, class ATYPE, class VTYPE>
  GRID_SIZE_TYPE GRID<DTYPE,ATYPE,VTYPE>::
  ComputeNumVertices(const VTYPE iv0, const VTYPE iv1) const
  {
    return(compute_num_vertices(Dimension(), AxisSize(), iv0, iv1));
  }

  template <class DTYPE, class ATYPE, class VTYPE>
  GRID_SIZE_TYPE GRID<DTYPE,ATYPE,VTYPE>::ComputeNumCubes() const
  {
    return(compute_num_grid_cubes(Dimension(), AxisSize()));
  }

  template <class DTYPE, class ATYPE, class VTYPE>
  template <class GTYPE>
  VTYPE GRID<DTYPE,ATYPE,VTYPE>::ComputeVertexIndex(const GTYPE * coord) const
  {
    VTYPE iv = compute_vertex_index(coord, Dimension(), AxisSize());
  }

  template <class DTYPE, class ATYPE, class VTYPE>
  template <class GTYPE>
  void GRID<DTYPE,ATYPE,VTYPE>::
  ComputeCoord(const VTYPE iv, GTYPE * coord) const
  // Precondition: coord[] is preallocated to length at least dimension
  {
    compute_coord(iv, Dimension(), AxisSize(), coord);
  }

  template <class DTYPE, class ATYPE, class VTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE>::
  CompareSize(const DTYPE dimension, const ATYPE * axis_size)
  // return true if grid dimension and axis size match parameters
  {
    if (dimension != this->Dimension()) { return(false); };
    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->AxisSize(d)) { return(false); };
    }

    return(true);
  }

  template <class DTYPE, class ATYPE, class VTYPE> 
  bool GRID<DTYPE,ATYPE,VTYPE>::Check
  (const DTYPE dimension, const ATYPE * axis_size, IJK::ERROR & error)
  // return true if grid dimension and axis size match parameter dimension
  //      and parameter axis size, respectively
  {
    if (dimension != this->dimension) {
      error.AddMessage("Incorrect grid dimension ", this->dimension, ".");
      return(false);
    }

    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->axis_size[d]) {
	error.AddMessage("Illegal axis size[", d, "] = ", 
			 this->axis_size[d], ".");
	return(false);
      }
    }

    VTYPE num_vertices = 
      compute_num_grid_vertices<VTYPE>(dimension, axis_size);

    if (num_vertices != this->num_vertices) {
      error.AddMessage("Incorrect number of grid vertices ", 
		       this->num_vertices, ".");
      return(false);
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID_BASE MEMBER FUNCTIONS
  // **************************************************

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE>::Replace
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid2, 
   const VTYPE iv0, const VTYPE iv1)
  // replace region from iv0 to iv1 with values from scalar_grid2
  // Precondition: scalar_grid2 has exactly the same dimension
  //   and axis_size as the current scalar grid.
  {
    const DTYPE dimension = this->dimension;
    const ATYPE * axis_size = this->axis_size;
    const STYPE * scalar2 = scalar_grid2.ScalarPtrConst();
    IJK::PROCEDURE_ERROR error("SCALAR_GRID_BASE::Replace");

    if (!Check(scalar_grid2.Dimension(), scalar_grid2.AxisSize(), error))
      { throw error; };

    if (dimension <= 0) { return; }

    ATYPE c0 = iv0%(this->AxisSize(0));
    ATYPE c1 = iv1%(this->AxisSize(0));

    // iv2 is the projection of iv1 onto the hyperplane orthogonal 
    // to axis 0 through iv0
    VTYPE iv2 = iv1 - c1 + c0;
  
    GRID_SIZE_TYPE num_slice_vertices =
      compute_num_grid_vertices(dimension, axis_size, iv0, iv2);

    // vlist contains a region slice orthogonal to the 0 axis
    VTYPE * vlist = new VTYPE[num_slice_vertices];
    get_grid_vertices(dimension, axis_size, iv0, iv2, vlist);

    for (ATYPE c = c0; c <= c1; c++) {
      for (VTYPE j = 0; j < num_slice_vertices; j++) {
	VTYPE jv = vlist[j];
	scalar[jv] = scalar2[jv];
      }

      // go to next slice
      for (VTYPE j = 0; j < num_slice_vertices; j++) 
	{ vlist[j]++; }
    }

    delete [] vlist;
  }

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID MEMBER FUNCTIONS
  // **************************************************

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  SCALAR_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  SCALAR_GRID():
    SCALAR_GRID_BASE <DTYPE,ATYPE,VTYPE,STYPE> (0, NULL)
  // constructor
  {
    this->scalar = NULL;
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  SCALAR_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  SCALAR_GRID
  (const DTYPE dimension, const ATYPE * axis_size) :
    SCALAR_GRID_BASE <DTYPE,ATYPE,VTYPE,STYPE> (dimension, axis_size)
  // constructor
  {
    this->scalar = new STYPE[this->NumVertices()];
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void SCALAR_GRID<DTYPE,ATYPE,VTYPE,STYPE>::FreeAll()
  {
    delete [] this->scalar;
    this->scalar = NULL;
    GRID<DTYPE,ATYPE,VTYPE>::FreeAll();
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void SCALAR_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  SetSize(const DTYPE dimension, const ATYPE * axis_size)
  {
    GRID_SIZE_TYPE numv = this->NumVertices();
    if (this->scalar == NULL || !CompareSize(dimension, axis_size)) {
      GRID<DTYPE,ATYPE,VTYPE>::Set(dimension, axis_size);
      if (this->scalar == NULL || numv != this->NumVertices()) {
	if (this->scalar != NULL) { delete [] this->scalar; };
	this->scalar = new STYPE[this->NumVertices()];
      }
    }
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void SCALAR_GRID<DTYPE,ATYPE,VTYPE,STYPE>::Subsample
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
   const ATYPE resolution)
  {
    const DTYPE dimension = scalar_grid.Dimension();
    VTYPE axis_increment[dimension];
    VTYPE vprimary[dimension];
    ATYPE subsampled_axis_size[dimension];
    IJK::PROCEDURE_ERROR error("SCALAR_GRID::Subsample");

    for (DTYPE d = 0; d < dimension; d++) {
      subsampled_axis_size[d] = 
	compute_subsample_size(scalar_grid.AxisSize(d), resolution);
    }

    SetSize(dimension, subsampled_axis_size);

    if (this->NumVertices() < 1) { return; };

    compute_increment(dimension, scalar_grid.AxisSize(), axis_increment);

    for (DTYPE d = 0; d < dimension; d++) { vprimary[d] = 0; };

    VTYPE iv = 0;
    while (vprimary[0] < scalar_grid.NumVertices()) {
      VTYPE jv = vprimary[0];
      this->scalar[iv] = scalar_grid.Scalar(jv);
      iv++;

      DTYPE d = 0;
      while (d+1 < dimension && 
	     (vprimary[d] + resolution*axis_increment[d] >=
	      vprimary[d+1] + axis_increment[d+1])) {
	d++;
      }

      // go to next subspace of dimension d
      vprimary[d] = vprimary[d] + resolution*axis_increment[d];

      for (DTYPE d2 = 0; d2 < d; d2++) 
	{ vprimary[d2] = vprimary[d]; };
    }

    if (iv != this->NumVertices()) {
      error.AddMessage("Programming error.  Subsampled ", iv, 
		       " grid vertices.");
      error.AddMessage("Subsampled grid should have ", this->NumVertices(),
		       " vertices.");
      throw error;
    }

  }

  // ******************************************************
  // TEMPLATE CLASS SCALAR_GRID_WRAPPER MEMBER FUNCTIONS
  // ******************************************************

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  SCALAR_GRID_WRAPPER<DTYPE,ATYPE,VTYPE,STYPE>::
  SCALAR_GRID_WRAPPER
  (const DTYPE dimension, const ATYPE * axis_size, STYPE * scalar):
    SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE>(dimension, axis_size)
  { this->scalar = scalar; };

  // **************************************************
  // TEMPLATE CLASS MINMAX_BASE MEMBER FUNCTIONS
  // **************************************************

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>::MINMAX_BASE()
  // constructor
  {
    scalar_min = NULL;
    scalar_max = NULL;
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>::MINMAX_BASE
  (const DTYPE dimension, const ATYPE * axis_size) :
    GRID<DTYPE,ATYPE,VTYPE>(dimension, axis_size)
  // constructor
  {
    scalar_min = new STYPE[this->NumVertices()];
    scalar_max = new STYPE[this->NumVertices()];
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>::FreeAll()
  {
    delete [] scalar_min;
    delete [] scalar_max;
    scalar_min = NULL;
    scalar_max = NULL;
    GRID<DTYPE,ATYPE,VTYPE>::FreeAll();
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>::
  SetSize(const DTYPE dimension, const ATYPE * axis_size)
  {
    GRID_SIZE_TYPE numv = this->NumVertices();
    if (this->scalar_min == NULL || !CompareSize(dimension, axis_size)) {
      GRID<DTYPE,ATYPE,VTYPE>::Set(dimension, axis_size);
      if (this->scalar_min == NULL || numv != this->NumVertices()) {
	if (this->scalar_min != NULL) { delete [] this->scalar_min; };
	if (this->scalar_max != NULL) { delete [] this->scalar_max; };
	this->scalar_min = new STYPE[this->NumVertices()];
	this->scalar_max = new STYPE[this->NumVertices()];
      }
    }
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>::
  Copy(const STYPE * scalar_min, const STYPE * scalar_max)
  {
    const GRID_SIZE_TYPE numv = this->NumVertices();
    IJK::PROCEDURE_ERROR error("MINMAX_BASE::Copy");

    if (numv > 0 && (scalar_min == NULL || scalar_max == NULL)) {
      error.AddMessage("Programming error. Memory for scalar_min or scalar_max has not been allocated.");
      throw error;
    }

    std::copy(scalar_min, scalar_min+numv, this->scalar_min);
    std::copy(scalar_max, scalar_max+numv, this->scalar_max);
  }

  // **************************************************
  // TEMPLATE CLASS MINMAX_GRID MEMBER FUNCTIONS
  // **************************************************

  /// Compute minimum and maximum of each region and store in primary vertices
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  OverwriteWithMinMax(const ATYPE region_edge_length)
  {
    ATYPE elength[this->dimension];
    ATYPE axis_increment[this->dimension];
    IJK::PROCEDURE_ERROR error("OverwriteWithMinMax");

    if (region_edge_length <= 0) {
      error("Programming error.  Region edge length must be positive.");
      throw error;
    }

    compute_increment(*this, axis_increment);

    for (DTYPE codim = 0; codim < this->Dimension(); codim++) {
      DTYPE d = this->Dimension()-codim-1;

      for (int d2 = 0; d2 < this->dimension; d2++) {
	if (d2 < d) { elength[d2] = 1; }
	else if (d2 > d) { elength[d2] = region_edge_length; }
	else { elength[d2] = this->axis_size[d2]; }
      }

      GRID_SIZE_TYPE subsample_size = 
	compute_subsample_size(this->dimension, this->axis_size, elength);

      IJK::ARRAY<VTYPE> vlist(subsample_size);
      subsample_grid_vertices(this->dimension, this->axis_size, elength,
			      vlist.Ptr());

      ATYPE num_subsample_vertices_along_axis =
	compute_subsample_size(AxisSize(d), region_edge_length);

      for (ATYPE j = 0; j < num_subsample_vertices_along_axis; j++) {
	ATYPE num_edges = region_edge_length;
	if ((j+1)*region_edge_length +1 > AxisSize(d)) {
	  num_edges = AxisSize(d) - j*region_edge_length - 1;
	}

	VTYPE facet_increment = (j*region_edge_length)*axis_increment[d];
	for (VTYPE k = 0; k < subsample_size; k++) {
	  VTYPE iv0 = vlist[k]+ facet_increment;
	  VTYPE iv1 = iv0;
	  for (ATYPE j2 = 1; j2 <= num_edges; j2++) {
	    iv1 += axis_increment[d];

	    STYPE min1 = this->scalar_min[iv1];
	    STYPE max1 = this->scalar_max[iv1];

	    if (this->scalar_min[iv0] > min1)
	      { this->scalar_min[iv0] = min1; };

	    if (this->scalar_max[iv0] < max1)
	      { this->scalar_max[iv0] = max1; };
	  }
	}
      }
    }
  }

  /// Compute minimum and maximum of each region and store in primary vertices
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  OverwriteWithMinMax(const ATYPE region_edge_length,
		      const ATYPE offset_edge_length)
  // Precondition: offset_edge_length < region_edge_length
  {
    ATYPE elength[this->dimension];
    ATYPE axis_increment[this->dimension];
    IJK::PROCEDURE_ERROR error("OverwriteWithMinMax");

    if (region_edge_length <= 0) {
      error("Programming error.  Region edge length must be positive.");
      throw error;
    }

    if (offset_edge_length >= region_edge_length) {
      error.AddMessage("Programming error.  Offset edge length must be less than region edge length.");
      throw error;
    }

    compute_increment(*this, axis_increment);

    for (DTYPE codim = 0; codim < this->Dimension(); codim++) {
      DTYPE d = this->Dimension()-codim-1;

      for (int d2 = 0; d2 < this->dimension; d2++) {
	if (d2 < d) { elength[d2] = 1; }
	else if (d2 > d) { elength[d2] = region_edge_length; }
	else { elength[d2] = this->axis_size[d2]; }
      }

      GRID_SIZE_TYPE subsample_size = 
	compute_subsample_size(this->dimension, this->axis_size, elength);

      IJK::ARRAY<VTYPE> vlist(subsample_size);
      subsample_grid_vertices(this->dimension, this->axis_size, elength,
			      vlist.Ptr());

      ATYPE num_subsample_vertices_along_axis =
	compute_subsample_size(AxisSize(d), region_edge_length);

      for (ATYPE j = 0; j < num_subsample_vertices_along_axis; j++) {
	ATYPE istart = 0;
	if (j > 0) { istart = j*region_edge_length - offset_edge_length; };
	ATYPE iend = (j+1)*region_edge_length+offset_edge_length+1;
	if (iend > AxisSize(d)) { iend = AxisSize(d); };
	ATYPE ionto = j*region_edge_length;

	for (VTYPE k = 0; k < subsample_size; k++) {
	  VTYPE iv0 = vlist[k]+ ionto*axis_increment[d];
	  VTYPE iv1 = vlist[k]+ istart*axis_increment[d];
	  for (ATYPE j2 = istart; j2 < iend; j2++) {

	    STYPE min1 = this->scalar_min[iv1];
	    STYPE max1 = this->scalar_max[iv1];

	    if (this->scalar_min[iv0] > min1)
	      { this->scalar_min[iv0] = min1; };

	    if (this->scalar_max[iv0] < max1)
	      { this->scalar_max[iv0] = max1; };

	    iv1 += axis_increment[d];
	  }
	}
      }
    }
  }

  /// Project scalar grid onto minmax grid
  /// Compute min and max of projection
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  ProjectMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE num_region_edges, const ATYPE facet_index, 
   const ATYPE axis_increment)
  // num_region_edges = number of grid edges per region edge
  // facet_index = index of facet.  Range [0..axis_size-1]
  // axis_size = axis size of axis (dimension-1)
  // axis_increment = axis increment of axis (dimension-1)
  // Precondition: minmax_grid.Dimension()+1 == dimension
  //               minmax_grid.AxisSize(d) == axis_size[d]
  //                   for all d < minmax_grid.Dimension()
  {
    // d_last = last dimension of scalar_grid
    const DTYPE d_last = this->Dimension();

    ATYPE num_edges = num_region_edges;
    ATYPE axis_size_d_last = axis_size[d_last];
    if (facet_index + num_region_edges + 1 > axis_size_d_last) {
      num_edges = axis_size_d_last - facet_index - 1;
    }

    VTYPE facet_increment = facet_index*axis_increment;
    for (VTYPE k = 0; k < this->NumVertices(); k++) {
      VTYPE iv0 = k + facet_increment;
      STYPE s0 = scalar[iv0];
      this->scalar_min[k] = s0;
      this->scalar_max[k] = s0;
      VTYPE iv1 = iv0;
      for (ATYPE j2 = 1; j2 <= num_edges; j2++) {
	iv1 += axis_increment;

	STYPE s1 = scalar[iv1];

	if (this->scalar_min[k] > s1) { this->scalar_min[k] = s1; };
	if (this->scalar_max[k] < s1) { this->scalar_max[k] = s1; };
      }
    }
  }

  /// Project scalar grid onto minmax grid
  /// Compute min and max of projection
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  ProjectMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE istart, const ATYPE iend, const ATYPE ionto,
   const ATYPE axis_increment)
  // istart = starting hyperplane
  // iend = ending hyperplane
  // ionto = project onto hyperplane ionto
  // facet_index = index of facet.  Range [0..axis_size-1]
  // axis_size = axis size of axis (dimension-1)
  // axis_increment = axis increment of axis (dimension-1)
  // Precondition: minmax_grid.Dimension()+1 == dimension
  //               minmax_grid.AxisSize(d) == axis_size[d]
  //                   for all d < minmax_grid.Dimension()
  {
    for (VTYPE k = 0; k < this->NumVertices(); k++) {
      VTYPE iv0 = k + ionto*axis_increment;
      STYPE s0 = scalar[iv0];
      this->scalar_min[k] = s0;
      this->scalar_max[k] = s0;
      VTYPE iv1 = k + istart*axis_increment;
      for (ATYPE j2 = istart; j2 < iend; j2++) {
	STYPE s1 = scalar[iv1];

	if (this->scalar_min[k] > s1) { this->scalar_min[k] = s1; };
	if (this->scalar_max[k] < s1) { this->scalar_max[k] = s1; };

	iv1 += axis_increment;
      }
    }
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_GRID<DTYPE,ATYPE,VTYPE,STYPE>::
  ProjectMinMax
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
   const ATYPE num_region_edges, const ATYPE facet_index, 
   const ATYPE axis_increment)
  {
    ProjectMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
		  scalar_grid.ScalarPtrConst(), num_region_edges,
		  facet_index, axis_increment);
  }



  // **************************************************
  // TEMPLATE CLASS MINMAX_REGIONS MEMBER FUNCTIONS
  // **************************************************

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>::
  MINMAX_REGIONS()
  // constructor
  {
    region_edge_length = 0;
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>::
  MINMAX_REGIONS
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length) :
    MINMAX_BASE<DTYPE,ATYPE,VTYPE,STYPE>(dimension, axis_size)
  // constructor
  {
    this->region_edge_length = region_edge_length;
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE> 
  MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>::
  ~MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>()
  // destructor
  {
    region_edge_length = 0;
  }

  /// Compute min and max of each region
  /// Creates a grid of regions and stores min and max for each region
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>::ComputeMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE region_edge_length)
  {
    ATYPE num_regions_along_axis[dimension];

    this->region_edge_length = region_edge_length;

    for (DTYPE d = 0; d < dimension; d++) {
      num_regions_along_axis[d] = 
	compute_num_regions_along_axis
	(axis_size[d], region_edge_length);
    }

    SetSize(dimension, num_regions_along_axis);

    compute_region_minmax
      (dimension, axis_size, scalar, region_edge_length, 
       this->scalar_min, this->scalar_max);
  }

  /// Compute min and max of each region
  /// Creates a grid of regions and stores min and max for each region
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>::ComputeMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE region_edge_length, const ATYPE offset_edge_length)
  {
    ATYPE num_regions_along_axis[dimension];

    this->region_edge_length = region_edge_length;

    for (DTYPE d = 0; d < dimension; d++) {
      num_regions_along_axis[d] = 
	compute_num_regions_along_axis
	(axis_size[d], region_edge_length);
    }

    SetSize(dimension, num_regions_along_axis);

    compute_region_minmax
      (dimension, axis_size, scalar, region_edge_length, offset_edge_length,
       this->scalar_min, this->scalar_max);
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>::ComputeMinMax
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid, 
   const ATYPE region_edge_length)
  {
    ComputeMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
			scalar_grid.ScalarPtrConst(), region_edge_length);
  }

  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void MINMAX_REGIONS<DTYPE,ATYPE,VTYPE,STYPE>::ComputeMinMax
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid, 
   const ATYPE region_edge_length, const ATYPE offset_edge_length)
  {
    ComputeMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
			scalar_grid.ScalarPtrConst(), region_edge_length,
			offset_edge_length);
  }

}

#endif
