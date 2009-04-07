/// \file ijktable.h
/// Class containing a table of isosurface patches in a given polyhedron.
/// All 2^numv +/- patterns are stored in the table 
///   where numv = # polyhedron vertices.
/// Version 0.3.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2007, 2006, 2003, 2001 Rephael Wenger

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

#ifndef _IJKTABLE_
#define _IJKTABLE_

#include <iostream>
#include <vector>

#include "ijk.txx"

namespace IJKTABLE {

  typedef unsigned char u_char;
  typedef unsigned long u_long;
  typedef unsigned int u_int;
  typedef u_char ISOSURFACE_VERTEX_INDEX;
  typedef u_char EDGE_INDEX;
  typedef u_char FACET_INDEX;
  typedef int TABLE_INDEX;
  typedef int FACET;

  const int NO_VERTEX = -1;

//**************************************************
// ISOSURFACE TABLE POLYHEDRON
//**************************************************

class ISOSURFACE_TABLE_POLYHEDRON {

 protected:
  int dimension;
  int num_vertices;
  int num_edges;
  int num_facets;
  int * vertex_coord;
  int * edge_endpoint;
  FACET * facet;
  void Init();           // initialize

 public:
  ISOSURFACE_TABLE_POLYHEDRON(const int d);  // constructor
  ~ISOSURFACE_TABLE_POLYHEDRON();            // destructor
  ISOSURFACE_TABLE_POLYHEDRON(const ISOSURFACE_TABLE_POLYHEDRON & init);
  // copy
  const ISOSURFACE_TABLE_POLYHEDRON & operator = 
    (const ISOSURFACE_TABLE_POLYHEDRON &);  // assign 


  // get functions
  int Dimension() const { return(dimension); };
  int NumVertices() const { return(num_vertices); };
  int NumEdges() const { return(num_edges); };
  int NumFacets() const { return(num_facets); };
  int NumFacetVertices(const FACET_INDEX jf) const;
    // jf = facet index
  int VertexCoord(const int iv, const int ic) const
    // iv = vertex index. ic = coordinate index.
    { return(vertex_coord[iv*dimension + ic]); };
  int EdgeEndpoint(const EDGE_INDEX ie, const int j) const
    // ie = edge index. j = 0 or 1.
    { return(edge_endpoint[int(ie)*2 + j]); };
  int MidpointCoord(const EDGE_INDEX ie, const int ic) const;
  FACET Facet(const FACET_INDEX jf) const
    // jf = facet index
    { return(facet[jf]); };
  bool FacetVertexFlag(const FACET_INDEX jf, const int iv) const
    { return(facet[jf] & ((1L) << iv)); };

  // set functions
  void SetDimension(const int d);
  void SetNumVertices(const int numv);
  void SetNumEdges(const int nume);
  void SetNumFacets(const int numf);
  void SetSize(const int numv, const int nume, const int numf)
    { SetNumVertices(numv); SetNumEdges(nume); SetNumFacets(numf); };
  void SetVertexCoord(const int iv, const int ic, const int coord);
  // Note: SetNumVertices or SetSize must be called before SetVertexCoord
  void SetEdge(const EDGE_INDEX ie, const int iv0, const int iv1);
  // Note: SetNumEdges or SetSize must be called before SetEdge
  void SetFacetVertex(const FACET_INDEX jf, const int iv, const bool in_facet);
  // Note: SetNumFacets must be called before SetFacetVertex

  // free memory
  void FreeAll();                            // free all memory

  // check functions
  bool CheckDimension() const;
  bool Check(IJK::ERROR & error_msg) const;

  // generate polyhedron
  void GenCube(const int cube_dimension);
  void GenSimplex(const int simplex_dimension);
  void GenPyramid(const int pyramid_dimension);
};

typedef ISOSURFACE_TABLE_POLYHEDRON * ISOSURFACE_TABLE_POLYHEDRON_PTR;


//**************************************************
// ISOSURFACE VERTEX
//**************************************************

class ISOSURFACE_VERTEX {

 public:
  typedef enum {VERTEX, EDGE, FACET, POINT} ISOSURFACE_VERTEX_TYPE;
  typedef float COORD_TYPE;

 protected:
  ISOSURFACE_VERTEX_TYPE vtype;
  int face;
  int num_coord;
  COORD_TYPE * coord;
  std::string label;
  bool is_label_set;

 public:
  ISOSURFACE_VERTEX();        // constructor
  ~ISOSURFACE_VERTEX();       // destructor

  // Get Functions
  ISOSURFACE_VERTEX_TYPE Type() const { return(vtype); };
  int Face() const { return(face); };
  COORD_TYPE Coord(const int d) const
  { return(coord[d]); };
  int NumCoord() const { return(num_coord); };
  std::string Label() const { return(label); };
  bool IsLabelSet() const { return(is_label_set); };

  // Set Functions
  void SetType(const ISOSURFACE_VERTEX_TYPE t) { vtype = t; };
  void SetFace(const int index) { face = index; };
  void SetNumCoord(const int numc);
  void SetCoord(const int ic, const COORD_TYPE c) 
  { coord[ic] = c; };
  void SetLabel(const std::string & s) 
  { label = s; is_label_set = true; };
};


//**************************************************
// ISOSURFACE TABLE
//**************************************************

class ISOSURFACE_TABLE {

 protected:
  class ISOSURFACE_TABLE_ENTRY {

  public:
    int num_simplices;
    ISOSURFACE_VERTEX_INDEX * simplex_vertex_list;
    ISOSURFACE_TABLE_ENTRY();                  // constructor
    ~ISOSURFACE_TABLE_ENTRY();                 // destructor

    bool Check(IJK::ERROR & error_msg) const;
    void FreeAll();                            // free all memory
  };


 public:
  typedef enum {BINARY, BASE3, NONSTANDARD} ENCODING;
  // Programmer's note: Update standard_encoding_name[] in ijktable.cxx
  //   whenever ENCODING is changed.

 protected:
  ENCODING encoding;           // type of encoding
  std::string encoding_name;   // string storing encoding name

  ISOSURFACE_TABLE_POLYHEDRON polyhedron;
  int simplex_dimension;
  ISOSURFACE_VERTEX * isosurface_vertex;
  int num_isosurface_vertices;
  ISOSURFACE_TABLE_ENTRY * entry;
  long num_table_entries;

  // max # vertices allowed for table polyhedron
  int max_num_vertices;

  bool is_table_allocated;

  // throw error if not enough isosurface vertices are allocated
  void CheckIsoVerticesAlloc
    (const char * procname, const int vstart, const int numv);

  // initialization routine
  void Init(const int dimension, const int simplex_dimension);

 public:
  // constructors
  ISOSURFACE_TABLE(const int d);
  ISOSURFACE_TABLE(const int dimension, const int simplex_dimension);

  ~ISOSURFACE_TABLE();                // destructor

  // get functions
  ENCODING Encoding() const { return(encoding); };
  std::string EncodingName() const { return(encoding_name); };
  static std::string StandardEncodingName(const ENCODING encoding);
  const ISOSURFACE_TABLE_POLYHEDRON & Polyhedron() const
    { return(polyhedron); };
  const ISOSURFACE_VERTEX & IsosurfaceVertex(const int i) const
    { return(isosurface_vertex[i]); };
  int Dimension() const { return(polyhedron.Dimension()); };
  int SimplexDimension() const { return(simplex_dimension); };
  int NumVerticesPerSimplex() const { return(SimplexDimension()+1); };
  int NumIsosurfaceVertices() const { return(num_isosurface_vertices); };
  int NumTableEntries() const { return(num_table_entries); };
  int NumSimplices(const TABLE_INDEX it) const
    { return(entry[it].num_simplices); };
  ISOSURFACE_VERTEX_INDEX SimplexVertex
    (const TABLE_INDEX it, const int is, const int k) const
    // it = table entry index. is = simplex index. k = k'th vertex
    // Note: Returns index of edge containing k'th simplex vertex
    { return(entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+k]); };
  int MaxNumVertices() const { return(max_num_vertices); };
  // Note: Even tables for polyhedra of this size are probably impossible 
  //   to compute/store
  bool IsTableAllocated() const
    { return(is_table_allocated); };

  // set polyhedron functions
  void SetDimension(const int d) { polyhedron.SetDimension(d); };
  void SetNumPolyVertices(const int numv) 
    { polyhedron.SetNumVertices(numv); };
  void SetNumPolyEdges(const int nume) { polyhedron.SetNumEdges(nume); };
  void SetNumPolyFacets(const int numf) { polyhedron.SetNumFacets(numf); };
  void SetPolySize(const int numv, const int nume, const int numf)
    { SetNumPolyVertices(numv); SetNumPolyEdges(nume); 
      SetNumPolyFacets(numf); };
  void SetPolyVertexCoord(const int iv, const int ic, const int coord)
    { polyhedron.SetVertexCoord(iv, ic, coord); };
  // Note: SetNumPolyVertices or SetPolySize must be called before 
  //   SetPolyVertexCoord
  void SetPolyEdge(const int ie, const int iv0, const int iv1)
    { polyhedron.SetEdge(ie, iv0, iv1); };
  // Note: SetNumPolyEdges or SetPolySize must be called before SetPolyEdge
  void SetPolyFacetVertex(const int jf, const int iv, const bool in_facet)
    { polyhedron.SetFacetVertex(jf, iv, in_facet); };
  // Note: SetPolyNumFacetVertices must be called before SetPolyFacetVertex
  void Set(const ISOSURFACE_TABLE_POLYHEDRON & polyhedron)
    { this->polyhedron = polyhedron; };

  // set isosurface vertices functions
  void SetNumIsosurfaceVertices(const int num_vertices);
  void SetIsoVertexType(const int i, 
			const ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE t) 
  { isosurface_vertex[i].SetType(t); };
  void SetIsoVertexFace(const int i, const int index) 
  { isosurface_vertex[i].SetFace(index); };
  void SetIsoVertexNumCoord(const int i, const int numc)
  { isosurface_vertex[i].SetNumCoord(numc); };
  void SetIsoVertexCoord(const int i, 
			 const int ic, const ISOSURFACE_VERTEX::COORD_TYPE c) 
  { isosurface_vertex[i].SetCoord(ic, c); };
  void SetIsoVertexLabel(const int i, const std::string & s) 
  { isosurface_vertex[i].SetLabel(s); };

  // store polyhedron vertices, edges or faces as isosurface vertices
  void StorePolyVerticesAsIsoVertices(const int vstart);
  void StorePolyEdgesAsIsoVertices(const int vstart);
  void StorePolyFacetsAsIsoVertices(const int vstart);

  // set isosurface table functions
  void SetSimplexDimension(const int d) { this->simplex_dimension = d; };
  void SetEncoding(const ENCODING encoding);
  void SetBinaryEncoding() { SetEncoding(BINARY); };
  void SetBase3Encoding() { SetEncoding(BASE3); };
  void SetNonstandardEncoding(const std::string & name);
  virtual void SetNumTableEntries(const int num_table_entries);
  void SetNumSimplices(const TABLE_INDEX it, const int nums);
  void SetSimplexVertex(const TABLE_INDEX it, const int is, 
			const int iv, const ISOSURFACE_VERTEX_INDEX isov);
  // generate polyhedron
  void GenCube(const int cube_dimension) 
    { polyhedron.GenCube(cube_dimension); };
  // Note: Cubes of dimension > 4 will have too many vertices
  void GenSimplex(const int simplex_dimension) 
    { polyhedron.GenSimplex(simplex_dimension); };
  void GenPyramid(const int pyramid_dimension) 
    { polyhedron.GenPyramid(pyramid_dimension); };

  // check functions
  bool CheckDimension(const int d) const;
  bool CheckDimension() const
    { return(CheckDimension(Dimension())); };
  bool CheckTable(IJK::ERROR & error_msg) const;
  bool Check(IJK::ERROR & error_msg) const;

  // free memory
  virtual void FreeAll();                     // free all memory
};

typedef ISOSURFACE_TABLE * ISOSURFACE_TABLE_PTR;

//**************************************************
// ISOSURFACE EDGE TABLE
//**************************************************

class ISOSURFACE_EDGE_TABLE:public ISOSURFACE_TABLE {

 protected:
  class ISOSURFACE_EDGE_TABLE_ENTRY {

  public:
    int num_edges;
    EDGE_INDEX * edge_endpoint_list;

    ISOSURFACE_EDGE_TABLE_ENTRY();                  // constructor
    ~ISOSURFACE_EDGE_TABLE_ENTRY();                 // destructor

    bool Check(IJK::ERROR & error_msg) const;
    void FreeAll();                            // free all memory
  };

 protected:
  ISOSURFACE_EDGE_TABLE_ENTRY * edge_entry;

  // initialization routine
  void Init(const int d);

 public:
  ISOSURFACE_EDGE_TABLE(const int d);      // constructor
  ~ISOSURFACE_EDGE_TABLE();                // destructor

  // get functions
  int NumEdges(const TABLE_INDEX it) const 
  { return(edge_entry[it].num_edges); };
  EDGE_INDEX EdgeEndpoint
    (const TABLE_INDEX it, const int ie, const int iend) const
    // it = table entry index. ie = edge index. iend = endpoint index (0 or 1)
  { return(edge_entry[it].edge_endpoint_list[2*ie+iend]); };

  // set isosurface table functions
  virtual void SetNumTableEntries(const int num_table_entries);

  // generate edge lists
  void GenEdgeLists();

  // check functions
  bool CheckTable(IJK::ERROR & error_msg) const;
  bool Check(IJK::ERROR & error_msg) const;

  // free memory
  virtual void FreeAll();                     // free all memory
};

typedef ISOSURFACE_EDGE_TABLE * ISOSURFACE_EDGE_TABLE_PTR;

//**************************************************
// UTILITY FUNCTIONS
//**************************************************

// calculate number of entries required in ISOSURFACE_TABLE
u_long calculate_num_entries(const int num_vert, const int num_colors);

// convert integer to base "base"
void convert2base(const u_long ival, const u_int base, int * digit,
		  const u_int max_num_digits);


//**************************************************
// ROUTINES FOR GENERATING POLYHEDRA
//**************************************************

void generate_prism(const ISOSURFACE_TABLE_POLYHEDRON & base_polyhedron,
		    ISOSURFACE_TABLE_POLYHEDRON & prism);

};

#endif
