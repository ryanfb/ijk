// \file ijkoctree.cxx
/// Scalar grid octree

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

#include <assert.h>
#include <climits>
#include <sstream>
#include <string>

#include "ijkoctree.h"
#include "ijkgrid.txx"

using namespace IJK;
using namespace IJKGRID;
using namespace IJKOCTREE;


// **************************************************
// LOCAL FUNCTIONS
// **************************************************

/// Return coordinate of specified vertex
namespace {

  VERTEX_INDEX compute_num_cubes
  (const int dimension, const AXIS_SIZE_TYPE length)
  // return number of cubes in a region
  {
    VERTEX_INDEX num_cubes = 1;
    for (int d = 0; d < dimension; d++) {
      num_cubes = num_cubes * length;
    }
    return(num_cubes);
  }

  void compute_vertex_increment
  (const int dimension, const AXIS_SIZE_TYPE * axis_size,
   const AXIS_SIZE_TYPE length, VERTEX_INDEX * vertex_increment)
  {
    GRID_COORD_TYPE coord[dimension];
    VERTEX_INDEX axis_increment[dimension];

    const VERTEX_INDEX numv = compute_num_region_vertices(dimension, length);

    for (int d = 0; d < dimension; d++)
      coord[d] = 0;

    compute_increment(dimension, axis_size, axis_increment);

    VERTEX_INDEX increment = 0;
    for (VERTEX_INDEX i = 0; i < numv; i++) {
      vertex_increment[i] = increment;
      increment++;
      coord[0]++;
      int d = 0;
      while (d < dimension && coord[d] > length) {
	increment = increment - coord[d]*axis_increment[d];
	coord[d] = 0;
	d++;
	if (d < dimension) {
	  coord[d]++;
	  increment = increment + axis_increment[d];
	}
      }
    }
  }

  void compute_cube_increment
  (const int dimension, const AXIS_SIZE_TYPE * axis_size,
   const AXIS_SIZE_TYPE length, VERTEX_INDEX * cube_increment)
  {
    GRID_COORD_TYPE coord[dimension];
    VERTEX_INDEX axis_increment[dimension];

    const VERTEX_INDEX num_cubes = compute_num_cubes(dimension, length);

    if (num_cubes == 0) { return; };

    for (int d = 0; d < dimension; d++)
      coord[d] = 0;

    compute_increment(dimension, axis_size, axis_increment);

    VERTEX_INDEX increment = 0;
    for (VERTEX_INDEX i = 0; i < num_cubes; i++) {
      cube_increment[i] = increment;
      increment++;
      coord[0]++;
      int d = 0;
      while (d < dimension && coord[d] >= length) {
	increment = increment - coord[d]*axis_increment[d];
	coord[d] = 0;
	d++;
	if (d < dimension) {
	  coord[d]++;
	  increment = increment + axis_increment[d];
	}
      }
    }
  }

  int compute_num_levels(const int dimension, const AXIS_SIZE_TYPE * axis_size)
  // Note: Individual cubes are not in octree
  // Leaf nodes in octree represent up to 2^dimension cubes
  {
    int num_levels = 0;

    if (dimension == 0) return(0);

    AXIS_SIZE_TYPE max_axis_size = axis_size[0];
    AXIS_SIZE_TYPE min_axis_size = axis_size[0];

    for (int d = 1; d < dimension; d++) {
      if (max_axis_size < axis_size[d])
	max_axis_size = axis_size[d];
      if (min_axis_size > axis_size[d])
	min_axis_size = axis_size[d];
    }

    if (min_axis_size < 2) return(0);

    num_levels = 1;
    AXIS_SIZE_TYPE l = 2;
    while (l+1 < max_axis_size) {
      l = l << 1;
      num_levels++;
    }

    return(num_levels);
  }

}


// **************************************************
// OCTREE DATA STRUCTURE
// **************************************************

OCTREE::OCTREE(const int dimension, const AXIS_SIZE_TYPE * axis_size):
  bon_octree_table(dimension, axis_size)
{
  PROCEDURE_ERROR error("OCTREE::OCTREE");

  // initialize
  this->dimension = 0;
  num_levels = 0;
  level = NULL;
  full_region_vertex_increment = NULL;
  num_full_region_vertices = 0;
  num_full_region_cubes = 0;
  is_minmax_set = false;

  if (dimension+1 > sizeof(int)*CHAR_BIT) {
    error.AddMessage("Octree dimension is too large.");
    throw error;
  }

  this->dimension = dimension;

  this->axis_size = new AXIS_SIZE_TYPE[dimension];
  for (int d = 0; d < dimension; d++) {
    this->axis_size[d] = axis_size[d];
  }

  num_cube_vertices = (1L << dimension);
  cube_vertex_increment = new VERTEX_INDEX[num_cube_vertices];
  compute_vertex_increment(dimension, axis_size, 1, cube_vertex_increment);
 
  num_full_region_cubes = (1L << dimension);
  num_full_region_vertices = compute_num_region_vertices(dimension, 2);
  full_region_vertex_increment = 
    new VERTEX_INDEX[num_full_region_vertices];
  compute_vertex_increment
    (dimension, axis_size, 2, full_region_vertex_increment);

  ComputeNumLevels();

  level = new LEVEL[num_levels];
  for (int ilevel = 0; ilevel < NumLevels(); ilevel++) {
    int height = NumLevels()-ilevel-1;
    AXIS_SIZE_TYPE region_length = 2;
    region_length = region_length << height;
    VERTEX_INDEX num_regions = 
      compute_num_regions(dimension, axis_size, region_length);
    level[ilevel].Create(dimension, height, num_regions);
  }

  CreateNodes();
}


void OCTREE::CreateNodes()
{
  PROCEDURE_ERROR error("OCTREE::OCTREE");

  if (num_levels == 0) { return; };

  // check that all levels have been created
  if (!CheckLevel(error)) { throw(error); };

  // store root node at level 0
  NODE_PTR root = level[0].node;
  root->SetRange(ComputeRootRange());
  KEY_TYPE root_key = level[0].Key(0);
  root->SetNumChildren(bon_octree_table.NumChildren(root_key));
  root->SetVertex0(0);

  for (int ilevel = 0; ilevel+1 < NumLevels(); ilevel++) {
    AXIS_SIZE_TYPE child_length = level[ilevel].child_length;

    int k = 0;
    NODE_PTR cur_node0 = level[ilevel].node;
    NODE_PTR next_node0 = level[ilevel+1].node;
    for (int j = 0; j < NumLevelNodes(ilevel); j++) {
      NODE_PTR cur_node = cur_node0 + j;
      KEY_TYPE key = level[ilevel].Key(j);
      cur_node->SetChild0(next_node0+k);
      VERTEX_INDEX cur_vertex0 = cur_node->Vertex0();

      for (int m = 0; m < cur_node->NumChildren(); m++) {
	NODE_PTR next_node = next_node0 + k + m;
	VERTEX_INDEX range = ComputeChildRange(ilevel, cur_node, m);
	VERTEX_INDEX child_vertex0 =
	  bon_octree_table.ChildVertex0(key, m, cur_vertex0, child_length);

	next_node->SetRange(range);
	KEY_TYPE child_key = level[ilevel+1].Key(k+m);
	next_node->SetNumChildren(bon_octree_table.NumChildren(child_key));
	next_node->SetVertex0(child_vertex0);
      }
      k = k + cur_node->NumChildren();
    }

    if (k != NumLevelNodes(ilevel+1)) {
      error.AddMessage("Programming error in creating nodes.");
      error.AddMessage("Error in creating nodes at level ", ilevel+1, ".");
      throw error;
    }
  }

  // set child pointer to NULL for each leaf node
  int ilast_level = NumLevels()-1;
  level[ilast_level].ClearChild0();
}


void OCTREE::ComputeNumLevels() 
// Note: Individual cubes are not in octree
// Leaf nodes in octree represent up to 2^dimension cubes
{
  num_levels = 0;

  if (dimension == 0) return;

  AXIS_SIZE_TYPE max_axis_size = axis_size[0];
  AXIS_SIZE_TYPE min_axis_size = axis_size[0];

  for (int d = 1; d < dimension; d++) {
    if (max_axis_size < axis_size[d])
      max_axis_size = axis_size[d];
    if (min_axis_size > axis_size[d])
      min_axis_size = axis_size[d];
  }

  if (min_axis_size < 2) return;

  num_levels = 1;
  AXIS_SIZE_TYPE l = 2;
  while (l+1 < max_axis_size) {
    l = l << 1;
    num_levels++;
  }
  
}


VERTEX_INDEX OCTREE::ComputeRootRange() const
{
  AXIS_SIZE_TYPE axis_range[dimension];
  const AXIS_SIZE_TYPE ONE = 1;

  for (int d = 0; d < dimension; d++) {
    if (axis_size[d] < 2) { return(0); };

    axis_range[d] = axis_size[d]-2;
  }

  AXIS_SIZE_TYPE max_axis_range = 0;
  for (int d = 0; d < dimension; d++) {
    if (max_axis_range < axis_range[d])
      max_axis_range = axis_range[d];
  }

  VERTEX_INDEX range = 0;
  VERTEX_INDEX bit = 1;
  AXIS_SIZE_TYPE l = max_axis_range;
  while (l != 0) {
    for (int d = 0; d < dimension; d++) {
      if ((axis_range[d] & ONE) == 1) { range = range | bit; };
      axis_range[d] = axis_range[d] >> 1;
      bit = bit << 1;
    }
    l = l >> 1;
  }

  return(range);
}

VERTEX_INDEX OCTREE::ComputeChildRange
(const int ilevel, const NODE_PTR node, const int m) const
{
  VERTEX_INDEX key = node->Range() >> level[ilevel].first_head_bit;
  VERTEX_INDEX range = node->Range() & level[ilevel].tailmask;
  range = range | bon_octree_table.ChildMask(key,m);
  range = range & level[ilevel].tailmask;
  return(range);
}

VERTEX_INDEX OCTREE::ComputeChildVertex0
(int ilevel, const NODE * node, const int m) const
{
  KEY_TYPE key = node->Range() >> level[ilevel].first_head_bit;
  VERTEX_INDEX vertex0 = 
    bon_octree_table.ChildVertex0(key, m, node->Vertex0(),
				  level[ilevel].child_length);
  return(vertex0);
}

SCALAR_TYPE OCTREE::ComputeMinCube
(const SCALAR_TYPE * scalar, const VERTEX_INDEX vertex0) const
// Precondition: vertex0 is the lowest coordinate vertex of a grid cube
{
  SCALAR_TYPE s = scalar[vertex0];
  for (int i = 1; i < num_cube_vertices; i++) {
    VERTEX_INDEX iv = vertex0 + cube_vertex_increment[i];
    if (s > scalar[iv])
      { s = scalar[iv]; };
  }

  return(s);
}


SCALAR_TYPE OCTREE::ComputeMaxCube
(const SCALAR_TYPE * scalar, const VERTEX_INDEX vertex0) const
// Precondition: vertex0 is the lowest coordinate vertex of a grid cube
{
  SCALAR_TYPE s = scalar[vertex0];
  for (int i = 1; i < num_cube_vertices; i++) {
    VERTEX_INDEX iv = vertex0 + cube_vertex_increment[i];
    if (s < scalar[iv])
      { s = scalar[iv]; };
  }

  return(s);
}

void OCTREE::SetMinMax(const MINMAX_REGIONS & minmax_regions)
// set min and max for each octree node
{
  const MASK_TYPE ZERO = 0;
  const MASK_TYPE ONE = 1;
  const KEY_TYPE num_keys = bon_octree_table.NumKeys();
  VERTEX_INDEX axis_increment[dimension];
  VERTEX_INDEX region_increment[num_keys*num_keys];

  is_minmax_set = true;

  VERTEX_INDEX ** region;
  typedef VERTEX_INDEX * VERTEX_INDEX_PTR;

  region = new VERTEX_INDEX_PTR[NumLevels()];
  for (int j = 0; j < NumLevels(); j++) {
    region[j] = new VERTEX_INDEX[NumLevelNodes(j)];
  }

  compute_increment(minmax_regions.Dimension(), minmax_regions.AxisSize(), 
		    axis_increment);

  region[0][0] = 0;
  for (int ilevel = 0; ilevel+1 < NumLevels(); ilevel++) {
    NODE_PTR node0 = level[ilevel].node;

    AXIS_SIZE_TYPE child_length = level[ilevel].child_length;

    for (KEY_TYPE i = 0; i < num_keys*num_keys; i++)
      { region_increment[i] = 0; };

    // compute region increment
    for (KEY_TYPE k = 0; k < num_keys; k++) {
      KEY_TYPE branch = 0;
      KEY_TYPE m = 0;
      while (branch < num_keys) {
	KEY_TYPE i = k*num_keys+m;
	for (int d = 0; d < dimension; d++) {
	  // get d'th bit
	  KEY_TYPE b = (branch >> d) & ONE;
	  if (b == 1) {
	    region_increment[i] += axis_increment[d]*(child_length/2);
	  }
	}
	m = m+1;

	// get next branch
	branch = branch+1;
	while (branch < num_keys && ((branch & ~k) != ZERO)) {
	  branch = branch+1;
	}
      }
    }

    for (int i = 0; i < NumLevelNodes(ilevel); i++) {
      const NODE * node = node0+i;
      for (int m = 0; m < node->NumChildren(); m++) {
	CONST_NODE_PTR child = node->Child(m);
	int k = child - level[ilevel+1].node;
	KEY_TYPE key = level[ilevel].Key(i);
	region[ilevel+1][k] = region[ilevel][i]+
	  region_increment[num_keys*key+m];
      }
    }
  }

  int leaf_level = NumLevels()-1;
  NODE_PTR leaf_node0 = level[leaf_level].node;
  for (int i = 0; i < NumLevelNodes(leaf_level); i++) {
    NODE_PTR leaf_node = leaf_node0+i;
    VERTEX_INDEX iregion = region[leaf_level][i];
    leaf_node->SetMin(minmax_regions.Min(iregion));
    leaf_node->SetMax(minmax_regions.Max(iregion));
  }

  // free array region[]
  for (int j = 0; j < NumLevels(); j++) {
    delete [] region[j];
  }
  delete [] region;


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

void OCTREE::SetMinMax(const SCALAR_TYPE * scalar)
// set min and max for each octree node
{
  is_minmax_set = true;

  if (NumLevels() == 0) return;

  MINMAX_REGIONS minmax_regions;

  minmax_regions.ComputeMinMax(Dimension(), AxisSize(), scalar, 2);

  SetMinMax(minmax_regions);
}

bool OCTREE::Check(ERROR & error) const
{
  const VERTEX_INDEX num_vertices = 
    compute_num_grid_vertices(dimension, axis_size);

  if (num_vertices == 0 && NumLevels() > 0) {
    error.AddMessage("Illegal value ", NumLevels(), " for number of levels.");
    error.AddMessage("Number of levels should be zero since number of vertices equals zero.");
    return(false);
  }

  if (!CheckLevel(error)) { return(false); };

  for (int ilevel = 0; ilevel < NumLevels(); ilevel++) {
    NODE_PTR node0 = level[ilevel].node;
    for (int j = 0; j < NumLevelNodes(ilevel); j++) {
      NODE_PTR node = node0 + j;
      if (node->Vertex0() < 0 || node->Vertex0() >= num_vertices) {
	error.AddMessage("Illegal vertex0 value for node ", j, " in level ", ilevel, ".");
	error.AddMessage("Vertex0 value is ", node->Vertex0(), 
			 ".  Should be between 0 and ", num_vertices-1, ".");
	return(false);
      }

      if (node->NumChildren() <= 0) {
	error.AddMessage("Illegal number of children for node ", j, " in level ", ilevel, ".");
	error.AddMessage("Node has ", node->NumChildren(), 
			 " children.  Value should be positive.");
	return(false);
      }
    }
  }

  for (int ilevel = 0; ilevel+1 < NumLevels(); ilevel++) {
    NODE_PTR node0 = level[ilevel].node;
    NODE_PTR next_node0 = level[ilevel+1].node;
    for (int j = 0; j < NumLevelNodes(ilevel); j++) {
      NODE_PTR node = node0 + j;
      CONST_NODE_PTR child_node = node->Child(0);
      DEGREE_TYPE num_children = node->NumChildren();

      VERTEX_INDEX num_next_level_nodes = NumLevelNodes(ilevel+1);

      if (child_node <  next_node0 ||
	  child_node+num_children > next_node0 + num_next_level_nodes) {
	error.AddMessage("Children of node ", j, " in level ", ilevel,
			 " are not in level ", ilevel+1, ".");
	return(false);
      }
    }
  }

  if (NumLevels() > 0) {
    int ileaf_level = NumLevels()-1;
    NODE_PTR leaf_node0 = level[ileaf_level].node;
    for (int j = 0; j < NumLevelNodes(ileaf_level); j++) {
      NODE_PTR leaf_node = leaf_node0 + j;
      if (leaf_node->Child(0) != NULL) {
	error.AddMessage("Leaf node ", j, " has non-Null pointer to child nodes.");
	error.AddMessage("Child node should be NULL.");
	return(false);
      }
    }
    
  }

  if (IsMinMaxSet()) {
    if (!CheckMinMax(error)) { return(false);  };
  }

  return(true);
}

bool OCTREE::CheckMinMax(ERROR & error) const
// check min/max of each internal node is min/max of children
{
  if (!IsMinMaxSet()) {
    error.AddMessage("Min/Max values have not been set.");
    return(false);
  }

  for (int ilevel = 0; ilevel+1 < NumLevels(); ilevel++) {
    NODE_PTR node0 = level[ilevel].node;
    for (int i = 0; i < NumLevelNodes(ilevel); i++) {
      NODE_PTR node = node0+i;
      if (node->MinValue() != node->ChildMin()) {
	error.AddMessage("Min value at ", i, "'th node of level ", ilevel, " does not equal min value of its children.");
	return(false);
      }

      if (node->MaxValue() != node->ChildMax()) {
	error.AddMessage("Max value at ", i, "'th node of level ", ilevel, " does not equal max value of its children.");
	return(false);
      }
    }
  }

  return(true);
}

bool OCTREE::CheckLevel(ERROR & error) const
{
  for (int ilevel = 0; ilevel < NumLevels(); ilevel++) {
    if (level[ilevel].node == NULL) {
      error.AddMessage("Level ", ilevel, " was not created.");
      return(false);
    }
  }
  return(true);
}

OCTREE::~OCTREE()
// destructor
{
  delete [] level;
  level = NULL;
  num_levels = 0;
  delete [] axis_size;
  axis_size = NULL;
  delete [] cube_vertex_increment;
  cube_vertex_increment = NULL;
  delete [] full_region_vertex_increment;
  full_region_vertex_increment = NULL;
  num_full_region_vertices = 0;
}

// **************************************************
// NODE STRUCTURE
// **************************************************

SCALAR_TYPE NODE::ChildMin() const
// return min of all children of node
// Precondition: node is an internal tree node
{
  SCALAR_TYPE child_min = (Child(0))->MinValue();

  for (int i = 1; i < NumChildren(); i++) {
    SCALAR_TYPE v = Child(i)->MinValue();
    if (child_min > v) { child_min = v; };
  }

  return(child_min);
}


SCALAR_TYPE NODE::ChildMax() const
// return max of all children of node
// Precondition: node is an internal tree node
{
  SCALAR_TYPE child_max = (Child(0))->MaxValue();

  for (int i = 1; i < NumChildren(); i++) {
    SCALAR_TYPE v = Child(i)->MaxValue();
    if (child_max < v) { child_max = v; };
  }

  return(child_max);
}

// **************************************************
// OCTTREE LEVEL STRUCTURE
// **************************************************

LEVEL::LEVEL()
// constructor
{
  height = 0;
  node = NULL;
  num_nodes = 0;
  first_head_bit = 0;
  tailmask = 0;
  child_length = 1;
}

void LEVEL::Create(const int dimension, 
		   const int height, const VERTEX_INDEX num_nodes)
// create level nodes and level information
// dimension = volume dimension 
// height = height of level in tree
// num_nodes = number of nodes in level
{
  const MASK_TYPE ONE = 1;

  if (node != NULL) {
    PROCEDURE_ERROR error("LEVEL::Create");
    error.AddMessage("Member function LEVEL::Create() can only be called once.");
    throw error;
  }

  this->height = height;
  node = new NODE[num_nodes];
  if (node != NULL) { this->num_nodes = num_nodes; }
  else { this->num_nodes = 0; };

  child_length = 1;
  child_length = (child_length << height);

  first_head_bit = dimension*height;

  tailmask = 0;
  for (KEY_TYPE i = 0; i < first_head_bit; i++) {
    tailmask = (tailmask << 1) | ONE;
  }
}

void LEVEL::ClearChild0()
{
  for (VERTEX_INDEX i = 0; i < num_nodes; i++) {
    node[i].SetChild0(NULL);
  }
}

LEVEL::~LEVEL()
// destructor
{
  delete [] node;
  node = NULL;
  num_nodes = 0;
}

// **************************************************
// BRANCH ON NEED OCTREE LOOKUP TABLE
// **************************************************

// Branch On Need (BON) Octree Lookup Table
// Lookup table for fast processing of BON information
// Branch On Need Octree based on paper by Wilhelms and Gelder,
//   ACM Transactions on Graphics, Vol. 11, 1992

BON_OCTREE_TABLE::BON_OCTREE_TABLE
(const int dimension, const AXIS_SIZE_TYPE * axis_size)
// constructor
{
  const MASK_TYPE ZERO = 0;
  const MASK_TYPE ONE = 1;
  VERTEX_INDEX axis_increment[dimension];
  MASK_TYPE direction_mask[dimension];

  num_keys = (1L << dimension);
  childmask = NULL;
  child_vertex_increment = NULL;
  num_children = NULL;

  childmask = new MASK_TYPE[num_keys*num_keys];
  child_vertex_increment = new VERTEX_INDEX[num_keys*num_keys];
  num_children = new int[num_keys];

  // compute num_children[]
  for (KEY_TYPE k = 0; k < num_keys; k++) {
    num_children[k] = 1;
    KEY_TYPE k2 = k;
    while (k2 != 0) {
      if ((k2 & ONE) == 1) 
	{ num_children[k] = num_children[k] << 1; }
      k2 = k2 >> 1;
    }
  }

  // compute num_levels in octree
  int num_levels = compute_num_levels(dimension, axis_size);

  // set direction_mask[0]
  direction_mask[0] = ZERO;
  for (int ilevel = 0; ilevel < num_levels; ilevel++) {
    direction_mask[0] = (direction_mask[0] << dimension) | ONE;
  }

  // set direction_mask[i] for i > 0
  for (int d = 1; d < dimension; d++) {
    direction_mask[d] = direction_mask[0] << d;
  }

  // compute axis_increment
  compute_increment(dimension, axis_size, axis_increment);

  // intialize all child_vertex_increments and childmasks to zero
  for (KEY_TYPE i = 0; i < num_keys*num_keys; i++) {
    child_vertex_increment[i] = 0;
    childmask[i] = ZERO;
  };

  // compute child_vertex_increment[] and childmask[]
  for (KEY_TYPE k = 0; k < num_keys; k++) {
    KEY_TYPE branch = 0;
    KEY_TYPE m = 0;
    while (branch < num_keys) {
      KEY_TYPE i = k*num_keys+m;
      for (int d = 0; d < dimension; d++) {
	// get d'th bit
	KEY_TYPE b = (branch >> d) & ONE;
	if (b == 1) {
	  child_vertex_increment[i] = 
	    child_vertex_increment[i] + axis_increment[d];
	}

	// get d'th bit
	b = ((~branch & k) >> d) & ONE;
	if (b == 1) {
	  childmask[i] = childmask[i] | direction_mask[d];
	}

      }
      m = m+1;

      // get next branch
      branch = branch+1;
      while (branch < num_keys && ((branch & ~k) != ZERO)) {
	branch = branch+1;
      }
    }
  }
}

BON_OCTREE_TABLE::~BON_OCTREE_TABLE()
// destructor
{
  num_keys = 0;
  delete [] childmask;
  childmask = NULL;
  delete [] child_vertex_increment;
  child_vertex_increment = NULL;
  delete [] num_children;
  num_children = NULL;
}

// **************************************************
// OCTREE STACK
// **************************************************

IJKOCTREE::OCTREE_STACK::OCTREE_STACK(const int num_levels) 
{
  size = 0;
  stack = new OCTREE_STACK_ELEMENT[num_levels];
  this->num_levels = num_levels;
};

IJKOCTREE::OCTREE_STACK::~OCTREE_STACK()
{
  size = 0;
  delete [] stack;
  stack = NULL;
  num_levels = 0;
}

