void OCTREE::SetMinMax(const SCALAR_GRID_BASE & scalar_grid)
// set min and max for each octree node
{
  is_minmax_set = true;

  if (NumLevels() == 0) return;

  MINMAX_GRID minmax_grid(scalar_grid.Dimension(), scalar_grid.AxisSize());
  minmax_grid.Copy(scalar_grid.ScalarPtrConst(), 
		   scalar_grid.ScalarPtrConst());
  minmax_grid.OverwriteWithMinMax(2);

  int ileaf_level = NumLevels()-1;
  NODE_PTR leaf_node0 = level[ileaf_level].node;
  AXIS_SIZE_TYPE region_length = level[ileaf_level].child_length;
  for (int i = 0; i < NumLevelNodes(ileaf_level); i++) {
    NODE_PTR leaf_node = leaf_node0 + i;

    VERTEX_INDEX vertex0 = leaf_node->Vertex0();
    leaf_node->SetMin(minmax_grid.Min(vertex0));
    leaf_node->SetMax(minmax_grid.Max(vertex0));
  }

  for (int j = 0; j+1 < NumLevels(); j++) {
    int ilevel = NumLevels()-j-2;
    NODE_PTR node0 = level[ilevel].node;
    for (int i = 0; i < NumLevelNodes(ilevel); i++) {
      NODE_PTR node = node0+i;
      node->SetMin(node->ChildMin());
      node->SetMax(node->ChildMax());
    }
  }
}

************************************************

  /* DEBUG
  /// Compute min and max scalar values of vertices 
  ///   in region with primary vertex iv0
  /// Note: Some region edges may be shorter than region_edge_length
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void compute_region_minmax
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
   const VTYPE iv0, const ATYPE region_edge_length,
   const ATYPE * axis_increment, STYPE & region_min, STYPE & region_max)
  {
    IJK::PROCEDURE_ERROR error("compute_region_minmax");
    const DTYPE dimension = scalar_grid.Dimension();
    const ATYPE * axis_size = scalar_grid.AxisSize();
    const DTYPE d_last = dimension-1;
    VTYPE vprimary[dimension];

    region_max = scalar_grid.Scalar(iv0);
    region_min = scalar_grid.Scalar(iv0);
    if (dimension < 1) { return; };

    for (DTYPE d = 0; d < dimension; d++)
      { vprimary[d] = iv0; };

    VTYPE v_end = scalar_grid.NumVertices();
    VTYPE v_end2 = iv0 + (region_edge_length+1)*axis_increment[d_last];
    if (v_end > v_end2) { v_end = v_end2; };

    while (vprimary[0] < v_end) {
      VTYPE iv = vprimary[0];
      STYPE s = scalar_grid.Scalar(iv);

      if (region_min > s) { region_min = s; };
      if (region_max < s) { region_max = s; };

      DTYPE d = 0;
      while (d+1 < dimension && 
	     ((vprimary[d] >= vprimary[d+1] + axis_increment[d+1]) ||
	      (vprimary[d] >=
	       vprimary[d+1] + region_edge_length*axis_increment[d])))
	{ d++; }

      // go to next subspace of dimension d
      vprimary[d] = vprimary[d] + axis_increment[d];

      for (DTYPE d2 = 0; d2 < d; d2++) 
	{ vprimary[d2] = vprimary[d]; };
    }

  }

  /// Compute minimum and maximum scalar value of each region
  template <class DTYPE, class ATYPE, class VTYPE, class STYPE>
  void compute_region_minmax
  (const SCALAR_GRID_BASE<DTYPE,ATYPE,VTYPE,STYPE> & scalar_grid,
   const ATYPE num_region_edges, STYPE * region_min, STYPE * region_max)
  // num_region_edges = number of grid edges per region edge
  // region_min[j] = min scalar value of region j
  // region_max[j] = max scalar value of region j
  //   Precondition: region_min[] and region_max[] are preallocated to size 
  //                 at least number of regions
  {
    const DTYPE dimension = scalar_grid.Dimension();
    const ATYPE * axis_size = scalar_grid.AxisSize();
    const GRID_SIZE_TYPE num_vertices = scalar_grid.NumVertices();
    const DTYPE d_last = dimension-1;
    VTYPE axis_increment[dimension];
    VTYPE vprimary[dimension];
    bool full[dimension];
    IJK::PROCEDURE_ERROR error("compute_region_minmax");

    if (num_region_edges <= 0) {
      error("Programming error.  Number of region edges must be positive.");
      throw error;
    }

    const GRID_SIZE_TYPE num_regions = 
      compute_num_regions(dimension, axis_size, num_region_edges);

    if (num_regions < 1 || num_vertices < 1) { return; };

    region_min[0] = scalar_grid.Scalar(0);
    region_max[0] = scalar_grid.Scalar(0);
    
    if (dimension < 1) { return;  };

    compute_increment(scalar_grid, axis_increment);

    const GRID_SIZE_TYPE num_region_vertices =
      compute_num_region_vertices(dimension, num_region_edges);

    VTYPE * region_vertex_increment = new VTYPE[num_region_vertices];
    compute_region_vertex_increment
      (dimension, axis_size, num_region_edges, region_vertex_increment);

    bool full0 = true;  // true if region (0,0,0,...) is full
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] < num_region_edges+1) { full0 = false; };
    }

    for (DTYPE d = 0; d < dimension; d++) {
      vprimary[d] = 0;
      full[d] = full0;
    };

    VTYPE iregion = 0;
    VTYPE v_end = 0;
    if (num_vertices > axis_increment[d_last]) {
      // Note: This condition should always be true
      v_end = num_vertices - axis_increment[d_last];
    }

    STYPE smin, smax, s;
    while (vprimary[d_last] < v_end) {

      if (!full[0]) {
	compute_region_minmax
	  (scalar_grid, vprimary[0], num_region_edges, axis_increment, 
	   smin, smax);
      }
      else {
	VTYPE iv0 = vprimary[0];
	smin = scalar_grid.Scalar(iv0+region_vertex_increment[0]);
	smax = smin;
	for (int i = 1; i < num_region_vertices; i++) {
	  s = scalar_grid.Scalar(iv0+region_vertex_increment[i]);
	  if (smin > s) { smin = s; };
	  if (smax < s) { smax = s; };
	}
      }

      region_min[iregion] = smin;
      region_max[iregion] = smax;
      iregion++;
				      
      DTYPE d = 0;
      while (d+1 < dimension && 
	     (vprimary[d] + (num_region_edges+1)*axis_increment[d] >=
	      vprimary[d+1] + axis_increment[d+1])) {
	d++;
      }

      // go to next subspace of dimension d
      vprimary[d] = vprimary[d] + num_region_edges*axis_increment[d];
      if (d+1 < dimension) {
	if (full[d+1] == false) { full[d] = false; }
	else if (vprimary[d] + num_region_edges*axis_increment[d] >=
		 vprimary[d+1] + axis_increment[d+1])
	  { full[d] = false; };
      }
      else {
	if (vprimary[d] + num_region_edges*axis_increment[d] >=
	    num_vertices)
	  { full[d] = false; };
      }

      for (DTYPE d2 = 0; d2 < d; d2++) {
	vprimary[d2] = vprimary[d]; 
	full[d2] = full[d];
      };
    }

    delete [] region_vertex_increment;
    region_vertex_increment = NULL;

    if (iregion != num_regions) {
      error.AddMessage("Programming error.  Computed min value for ", iregion,
		       " regions.");
      error.AddMessage("Number of regions is ", num_regions, ".");
      throw error;
    }
  }
  */


**************************************

/* OBSOLETE
// set min and max for each octree node
{
  is_minmax_set = true;

  if (NumLevels() == 0) return;

  int ileaf_level = NumLevels()-1;
  NODE_PTR leaf_node0 = level[ileaf_level].node;
  AXIS_SIZE_TYPE region_length = level[ileaf_level].child_length;
  for (int i = 0; i < NumLevelNodes(ileaf_level); i++) {
    NODE_PTR leaf_node = leaf_node0 + i;
    SCALAR_TYPE min_value;
    SCALAR_TYPE max_value;
    if (leaf_node->NumChildren() == NumFullRegionCubes()) {
      VERTEX_INDEX vertex0 = leaf_node->Vertex0();
      min_value = scalar[vertex0];
      max_value = scalar[vertex0];
      for (int i = 1; i < NumFullRegionVertices(); i++) {
	VERTEX_INDEX iv = vertex0+full_region_vertex_increment[i];
	SCALAR_TYPE s = scalar[iv];
	if (min_value > s) { min_value = s; };
	if (max_value < s) { max_value = s; };
      }
    }
    else {
      KEY_TYPE key = level[ileaf_level].Key(i);
      VERTEX_INDEX leaf_vertex0 = leaf_node->Vertex0();
      VERTEX_INDEX child_vertex0 =
	bon_octree_table.ChildVertex0(key, 0, leaf_vertex0, region_length);
      min_value = ComputeMinCube(scalar, child_vertex0);
      max_value = ComputeMaxCube(scalar, child_vertex0);
      
      for (int m = 1; m < leaf_node->NumChildren(); m++) {
	VERTEX_INDEX child_vertex0 =
	  bon_octree_table.ChildVertex0(key, m, leaf_vertex0, region_length);
	SCALAR_TYPE smin = ComputeMinCube(scalar, child_vertex0);
	SCALAR_TYPE smax = ComputeMaxCube(scalar, child_vertex0);
	if (smin < min_value) { min_value = smin; };
	if (smax > max_value) { max_value = smax; };
      };
    }

    leaf_node->SetMin(min_value);
    leaf_node->SetMax(max_value);
  }

  for (int j = 0; j+1 < NumLevels(); j++) {
    int ilevel = NumLevels()-j-2;
    NODE_PTR node0 = level[ilevel].node;
    for (int i = 0; i < NumLevelNodes(ilevel); i++) {
      NODE_PTR node = node0+i;
      node->SetMin(node->ChildMin());
      node->SetMax(node->ChildMax());
    }
  }
}
*/



