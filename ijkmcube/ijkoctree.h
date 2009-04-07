/// \file ijkoctree.h
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

#ifndef _IJKOCTREE_
#define _IJKOCTREE_

#include <string>
#include <vector>

#include "ijk.txx"
#include "ijkgrid.txx"

namespace IJKOCTREE {

  // SHOULD TURN THIS INTO TEMPLATES
  typedef float SCALAR_TYPE;
  typedef int GRID_COORD_TYPE;
  typedef int AXIS_SIZE_TYPE;
  typedef int VERTEX_INDEX;
  typedef VERTEX_INDEX MASK_TYPE;
  typedef int KEY_TYPE;
  typedef int DEGREE_TYPE;

  typedef IJKGRID::MINMAX_REGIONS<int,AXIS_SIZE_TYPE,VERTEX_INDEX,SCALAR_TYPE>
    MINMAX_REGIONS;

// **************************************************
// OCTTREE NODES
// **************************************************

  class MINMAX {
  protected:
    SCALAR_TYPE min_value;
    SCALAR_TYPE max_value;
    MINMAX(){};

  public:

    // get functions
    SCALAR_TYPE MinValue() const { return(min_value); };
    SCALAR_TYPE MaxValue() const { return(max_value); };

    // set functions
    void SetMin(const SCALAR_TYPE s) { min_value = s; };
    void SetMax(const SCALAR_TYPE s) { max_value = s; };
  };

  class NODE:public MINMAX {
  protected:
    const NODE * child0;
    DEGREE_TYPE num_children;
    VERTEX_INDEX range;
    VERTEX_INDEX vertex0;

  public:
    NODE(){};

    // get functions
    const NODE * Child(const int i) const { return(child0+i); };
    DEGREE_TYPE NumChildren() const { return(num_children); };
    VERTEX_INDEX Range() const { return(range); };
    VERTEX_INDEX Vertex0() const { return(vertex0); };

    // Precondition for ChildMin() and ChildMax()
    // Node must be an internal tree node (i.e. Child(0) != NULL)
    SCALAR_TYPE ChildMin() const;
    SCALAR_TYPE ChildMax() const;

    // set functions
    void SetRange(const VERTEX_INDEX range) { this->range = range; };
    void SetNumChildren(const int num_children)
    { this->num_children = num_children; };
    void SetChild0(const NODE * child) { this->child0 = child; };
    void SetVertex0(const VERTEX_INDEX v0) { vertex0 = v0; };
  };
  typedef NODE * NODE_PTR;
  typedef const NODE * CONST_NODE_PTR;

// **************************************************
// OCTTREE LEVEL STRUCTURE
// **************************************************

  class LEVEL {
  public:
    int height;
    NODE * node;
    VERTEX_INDEX num_nodes;
    KEY_TYPE first_head_bit;
    MASK_TYPE tailmask;
    AXIS_SIZE_TYPE child_length;

  public:
    LEVEL();
    ~LEVEL();
    void Create(const int dimension, 
		const int height, const VERTEX_INDEX num_nodes);
    void ClearChild0();  // set Child(0) to NULL for all nodes
                         // does not affect num_children

    // get functions
    KEY_TYPE Key(const int i) const
    { return(node[i].Range() >> first_head_bit); };
  };

// **************************************************
// BRANCH ON NEED OCTREE LOOKUP TABLE
// **************************************************

// Branch On Need (BON) Octree Lookup Table
// Lookup table for fast processing of BON information
// Branch On Need Octree based on paper by Wilhelms and Gelder,
//   ACM Transactions on Graphics, Vol. 11, 1992

  class BON_OCTREE_TABLE {

  protected:
    int dimension;
    KEY_TYPE num_keys;
    MASK_TYPE * childmask;
    VERTEX_INDEX * child_vertex_increment;
    DEGREE_TYPE * num_children;

  public:
    BON_OCTREE_TABLE(const int dimension, const AXIS_SIZE_TYPE * axis_size);
    ~BON_OCTREE_TABLE();

    // get functions
    KEY_TYPE NumKeys() const { return(num_keys); };
    MASK_TYPE ChildMask(const KEY_TYPE key, const DEGREE_TYPE m) const
    { return(childmask[key*num_keys + m]); };
    VERTEX_INDEX ChildVertex0
      (const KEY_TYPE key, const DEGREE_TYPE m, const VERTEX_INDEX vertex0, 
       const AXIS_SIZE_TYPE child_length) const
    { return(vertex0 + child_vertex_increment[key*num_keys+m]*child_length); };
    DEGREE_TYPE NumChildren(const KEY_TYPE key) const
    { return(num_children[key]); };
  };
  
// **************************************************
// OCTREE DATA STRUCTURE
// **************************************************

  class OCTREE {

  protected:
    int dimension;
    AXIS_SIZE_TYPE * axis_size;
    int num_cube_vertices;
    VERTEX_INDEX * cube_vertex_increment;
    int num_levels;
    LEVEL * level;
    bool is_minmax_set;

    // data for constructing nodes
    BON_OCTREE_TABLE bon_octree_table;
    VERTEX_INDEX * full_region_vertex_increment;
    VERTEX_INDEX num_full_region_vertices;
    VERTEX_INDEX num_full_region_cubes;

    void ComputeNumLevels();
    void CreateNodes();

  public:
    OCTREE(const int dimension, const AXIS_SIZE_TYPE * axis_size);
    ~OCTREE();

    // get functions
    const NODE * Root() const { return(level[0].node); };
    int Dimension() const { return(dimension); };
    const AXIS_SIZE_TYPE * AxisSize() const { return(axis_size); }
    VERTEX_INDEX NumCubeVertices() const
    { return(num_cube_vertices); };
    VERTEX_INDEX NumFullRegionCubes() const
    { return(num_full_region_cubes); };
    VERTEX_INDEX NumFullRegionVertices() const
    { return(num_full_region_vertices); };
    VERTEX_INDEX NumLevelNodes(const int i) const 
    { return(level[i].num_nodes); };
    int NumLevels() const { return(num_levels); };
    VERTEX_INDEX NumLeaves() const {
      if (NumLevels() == 0) { return(0); }
      else { return(NumLevelNodes(NumLevels()-1)); };
    };
    bool IsMinMaxSet() const { return(is_minmax_set); };

    // set functions
    void SetMinMax(const MINMAX_REGIONS & minmax);
    void SetMinMax(const SCALAR_TYPE * scalar);

    // compute functions
    VERTEX_INDEX ComputeChildRange
      (const int ilevel, const NODE_PTR node, const int i) const;
    VERTEX_INDEX ComputeChildVertex0
      (const int ilevel, const NODE * node, const int i) const;
    VERTEX_INDEX ComputeRootRange() const;
    SCALAR_TYPE ComputeMinCube
      (const SCALAR_TYPE * scalar, const VERTEX_INDEX vertex0) const;
    SCALAR_TYPE ComputeMaxCube
      (const SCALAR_TYPE * scalar, const VERTEX_INDEX vertex0) const;

    // check function
    bool CheckMinMax(IJK::ERROR & error) const;
    bool CheckLevel(IJK::ERROR & error) const;
    bool Check(IJK::ERROR & error) const;
  };

// **************************************************
// OCTREE STACK
// **************************************************

  struct OCTREE_STACK_ELEMENT {
    IJKOCTREE::CONST_NODE_PTR node;
    IJKOCTREE::CONST_NODE_PTR child_node;
    IJKOCTREE::DEGREE_TYPE num_unprocessed_children;
  };

  class OCTREE_STACK {

  protected:
    int size;
    int num_levels;
    OCTREE_STACK_ELEMENT * stack;

  public:
    OCTREE_STACK(const int num_levels);
    ~OCTREE_STACK();

    // get functions
    int Size() const { return(size); };
    bool IsEmpty() const { return(size == 0); };
    bool TopIsLeaf() const { return(size == num_levels); };
    bool NoNextChild() const {
      int itop = Size() - 1;
      return(stack[itop].num_unprocessed_children == 0);
    };
    CONST_NODE_PTR TopNode()
    { return(stack[Size()-1].node); };

    // push/pop functions
    void PushRoot(CONST_NODE_PTR root) {
      // Precondition: IsEmpty() == true
      stack[0].node = root;
      stack[0].child_node = root->Child(0);
      stack[0].num_unprocessed_children = root->NumChildren();
      size = 1;
    };

    void PushChild() {
      // Preconditions: IsEmpty() == false && IsLeaf() == false
      int ilevel = Size()-1;
      IJKOCTREE::CONST_NODE_PTR child_node = stack[ilevel].child_node;
      stack[size].node = child_node;
      stack[size].num_unprocessed_children = 
	child_node->NumChildren();
      stack[size].child_node = child_node->Child(0);
      size++;

      // change stack[ilevel].child_node to next child
      stack[ilevel].child_node++;
      stack[ilevel].num_unprocessed_children--;
    };
    void Pop() { size--; };
  };

};

#endif
