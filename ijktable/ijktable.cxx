/// \file ijktable.cxx
/// Class containing a table of isosurface patches in a given polyhedron.
/// Version 0.2.0

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

#include <ctype.h>
#include <limits>
#include <limits.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>

#include <set>
#include <vector>
#include <algorithm>

#include "ijktable.h"

using namespace IJK;
using namespace IJKTABLE;
using namespace std;

#ifndef LONG_BIT

#define LONG_BIT (CHAR_BIT * sizeof(long))

#endif

// local namespace
namespace {
  static char * standard_encoding_name[] =    
    { "BINARY", "BASE3", "NONSTANDARD" };
}

//**************************************************
// ISOSURFACE_TABLE_POLYHEDRON
//**************************************************

ISOSURFACE_TABLE_POLYHEDRON::ISOSURFACE_TABLE_POLYHEDRON(const int d)
// constructor
{
  Init();
  dimension = d;
}

ISOSURFACE_TABLE_POLYHEDRON::~ISOSURFACE_TABLE_POLYHEDRON()
// destructor
{
  FreeAll();
}

ISOSURFACE_TABLE_POLYHEDRON::ISOSURFACE_TABLE_POLYHEDRON
  (const ISOSURFACE_TABLE_POLYHEDRON & init)
// copy
{
  Init();

  *this = init;
}

const ISOSURFACE_TABLE_POLYHEDRON & ISOSURFACE_TABLE_POLYHEDRON::operator = 
  (const ISOSURFACE_TABLE_POLYHEDRON & right)  
// assign 
{
  if (&right != this) {         // avoid self-assignment
    FreeAll();
    dimension = right.Dimension();
    SetSize(right.NumVertices(), right.NumEdges(), right.NumFacets());

    // copy vertices
    for (int iv = 0; iv < NumVertices(); iv++)
      for (int ic = 0; ic < Dimension(); ic++)
	SetVertexCoord(iv, ic, right.VertexCoord(iv, ic));

    // copy edges
    for (int ie = 0; ie < NumEdges(); ie++)
      SetEdge(ie, right.EdgeEndpoint(ie, 0), right.EdgeEndpoint(ie, 1));

    // copy facets
    for (int jf = 0; jf < NumFacets(); jf++)
      facet[jf] = right.Facet(jf);
  };

  return *this;
}

void ISOSURFACE_TABLE_POLYHEDRON::FreeAll()
// free all memory
{
  num_vertices = 0;
  num_edges = 0;
  num_facets = 0;
  delete [] vertex_coord;
  vertex_coord = NULL;
  delete [] edge_endpoint;
  edge_endpoint = NULL;
  delete [] facet;
  facet = NULL;
}

void ISOSURFACE_TABLE_POLYHEDRON::Init()
// initialize
{
  dimension = 0;
  num_vertices = 0;
  num_edges = 0;
  num_facets = 0;
  vertex_coord = NULL;
  edge_endpoint = NULL;
  facet = NULL;
}

int ISOSURFACE_TABLE_POLYHEDRON::NumFacetVertices(const FACET_INDEX jf) const
// return number of facet vertices
// jf = facet index
{
  int numf = 0;
  for (int iv = 0; iv < NumVertices(); iv++) {
    if (FacetVertexFlag(jf, iv) != 0)
      numf++;
  };

  return(numf);
}

void ISOSURFACE_TABLE_POLYHEDRON::SetDimension(const int d)
// set polyhedron dimension
{
  FreeAll();
  num_vertices = num_edges = 0;

  dimension = d;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetNumVertices(const int numv)
// set number of vertices
// Must be called before setting polyhedron vertices
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetNumVertices";

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal polyhedron dimension.");

  FreeAll();
  num_vertices = num_edges = 0;

  if (numv == 0)
    throw PROCEDURE_ERROR(procname, "Number of vertices must be non-zero.");

  // Note that even if numv <= LONG_BIT, there may not be enough 
  //   memory to store the isosurface table.
  if (numv > LONG_BIT)
    throw PROCEDURE_ERROR
      (procname, "Number of polyhedron vertices is too large.");

  num_vertices = numv;
  vertex_coord = new int[num_vertices*Dimension()];
}

void ISOSURFACE_TABLE_POLYHEDRON::SetNumEdges(const int nume)
// set number of edges
// Must be called before setting polyhedron edges
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetNumEdges";

  delete [] edge_endpoint;
  edge_endpoint = NULL;
  num_edges = 0;

  if (!CheckDimension()) 
    throw PROCEDURE_ERROR(procname, "Illegal dimension.");

  if (NumVertices() == 0)
    throw PROCEDURE_ERROR
      (procname, "Number of vertices must be set before number of edges.");

  if (nume < 1)
    throw PROCEDURE_ERROR(procname, "Number of edges must be non-zero.");

  if (nume > std::numeric_limits<EDGE_INDEX>::max())
    throw PROCEDURE_ERROR
      (procname, "Number of polyhedron edges is too large.");

  num_edges = nume;
  edge_endpoint = new int[num_edges*2];
}


void ISOSURFACE_TABLE_POLYHEDRON::SetNumFacets(const int numf)
// set number of facets
// Must be called before setting polyhedron facets
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetNumFacets";

  delete [] facet;
  facet = NULL;
  num_facets = 0;

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal dimension.");

  if (NumVertices() == 0)
    throw PROCEDURE_ERROR
      (procname, "Number of vertices must be set before number of facets.");

  if (numf < 1)
    throw PROCEDURE_ERROR(procname, "Number of facets must be non-zero.");

  if (numf > std::numeric_limits<FACET_INDEX>::max())
    throw PROCEDURE_ERROR
      (procname, "Number of polyhedron facets is too large.");

  num_facets = numf;
  facet = new FACET[numf];

  if (facet == NULL)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for list of facets.");

  // initialize each facet to 0
  for (int jf = 0; jf < numf; jf++) {
    facet[jf] = 0;
  };

}

void ISOSURFACE_TABLE_POLYHEDRON::SetVertexCoord
  (const int iv, const int ic, const int coord)
// set polyhedron vertex coordinate
// iv = vertex index.  In range [0..NumVertices()-1].
// ic = coordinate index. In range [0..Dimension()-1].
// coord = coordinate.  Must be even.
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetVertexCoord";

  if (iv < 0 || iv >= NumVertices())
    throw PROCEDURE_ERROR(procname, "Illegal vertex index.");

  if (ic < 0 || ic >= Dimension())
    throw PROCEDURE_ERROR(procname, "Illegal vertex coordinate index.");

  if (coord%2 != 0)
    throw PROCEDURE_ERROR
      (procname, "Illegal vertex coordinate.  Vertex coordinate must be even.");

  if (vertex_coord == NULL)
    throw PROCEDURE_ERROR
      (procname, "Vertex coordinate memory not allocated.");

  vertex_coord[iv*Dimension() + ic] = coord;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetEdge
  (const EDGE_INDEX ie, const int iv0, const int iv1)
// set polyhedron edge coordinate
// ie = edge index.  In range [0..NumEdges()-1].
// iv0 = vertex 0 index.  In range [0..NumVertices()-1].
// iv1 = vertex 1 index.  In range [0..NumVertices()-1].
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetEdge";

  int ie_int = int(ie);
  if (ie_int < 0 || ie_int >= NumEdges())
    throw PROCEDURE_ERROR(procname, "Illegal edge index.");

  if (iv0 < 0 || iv0 > NumVertices() ||
      iv1 < 0 || iv1 > NumVertices())
    throw PROCEDURE_ERROR(procname, "Illegal vertex index.");

  if (edge_endpoint == NULL)
    throw PROCEDURE_ERROR(procname, "Edge endpoint memory not allocated.");

  edge_endpoint[2*int(ie)] = iv0;
  edge_endpoint[2*int(ie)+1] = iv1;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetFacetVertex
  (const FACET_INDEX jf, const int iv, const bool in_facet)
// set facet vertex to in or not in facet
// jf = facet index
// iv = vertex index
// in_facet = true if vertex iv is in facet
// Note: SetNumFacets must be called before SetFacetVertex
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetFacetVertex";

  int jf_int = int(jf);
  if (jf_int < 0 || jf_int >= NumFacets())
    throw PROCEDURE_ERROR(procname, "Illegal facet index.");

  if (iv < 0 || iv >= NumVertices()) {
    throw PROCEDURE_ERROR(procname, "Illegal vertex index.");
  };

  if (facet == NULL)
    throw PROCEDURE_ERROR
      (procname, "Number of facets must be set before facet vertices.");

  long mask = 1L << iv;
  if (in_facet) {
    facet[jf_int] = facet[jf_int] | mask;
  }
  else {
    mask = ~mask;
    facet[jf_int] = facet[jf_int] & mask;
  };
}

int ISOSURFACE_TABLE_POLYHEDRON::MidpointCoord
  (const EDGE_INDEX ie, const int ic) const
// return ic'th coordinate of midpoint of edge ie
// Note: vertice coordinates must all be even so midpoint coordinate
//   is an integer
{
  int iv0 = EdgeEndpoint(ie, 0);
  int iv1 = EdgeEndpoint(ie, 1);
  int coord0 = VertexCoord(iv0, ic);
  int coord1 = VertexCoord(iv1, ic);
  
  return((coord0+coord1)/2);
}

void ISOSURFACE_TABLE_POLYHEDRON::GenCube(const int cube_dimension)
// generate a polyhedron
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::GenCube";

  dimension = cube_dimension;

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal cube dimension.");

  int numv = 1L << Dimension();
  int nume = (numv * Dimension()) / 2;
  int numf = 2*Dimension();

  SetSize(numv, nume, numf);

  // set vertex coordinates
  for (int iv = 0; iv < NumVertices(); iv++) {
    long mask = 1L;
    for (int ic = 0; ic < Dimension(); ic++) {
      int bit = mask & iv;
      int coord = 0;
      if (bit != 0)
	coord = 2;
      SetVertexCoord(iv, ic, coord);
      mask = mask << 1;
    };
  };

  // generate edges in lexicographic order
  int ie = 0;
  long control = 0;
  while (ie < NumEdges()) {

    // find first 0 bit in control
    int ic = 0;
    long mask = 1L;
    while ((mask & control) != 0) {
      mask = mask << 1;
      ic++;
    };

    // find start vertex by stripping ic bits from control
    int start = control;
    start = start >> ic;
    start = start << ic;

    for (int i = 0; i < (1L << ic); i++) {
      int iv0 = start + i;
      int iv1 = iv0 + (1L << ic);
      SetEdge(ie, iv0, iv1);
      ie++;
    };

    control++;
  };

  if (control+1 != (1L << Dimension()))
    throw PROCEDURE_ERROR(procname, "Programming error in edge generation.");

  // generate facets
  int num_vertices_per_facet = NumVertices()/2;
  long mask = 1L;
  for (int ic = 0; ic < Dimension(); ic++) {
    int jf0 = 2*ic;
    int jf1 = 2*ic+1;

    for (int iv = 0; iv < NumVertices(); iv++) {
      long bit = mask & iv;
      if (bit == 0) {
	SetFacetVertex(jf0, iv, true);
      }
      else {
	SetFacetVertex(jf1, iv, true);
      };
    };

    if (NumFacetVertices(jf0) != num_vertices_per_facet ||
	NumFacetVertices(jf1) != num_vertices_per_facet) {
      throw PROCEDURE_ERROR
	(procname, "Programming error in facet generation.");
    };

    mask = mask << 1;
  };

}


void ISOSURFACE_TABLE_POLYHEDRON::GenSimplex(const int simplex_dimension)
// generate a simplex
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::GenSimplex";

  dimension = simplex_dimension;

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal simplex dimension.");

  int numv = Dimension() + 1;
  int nume = (numv * Dimension()) / 2;
  int numf = Dimension() + 1;

  SetSize(numv, nume, numf);

  // initialize all vertex coordinates to 0
  for (int iv = 0; iv < NumVertices(); iv++)
    for (int ic = 0; ic < Dimension(); ic++)
      SetVertexCoord(iv, ic, 0);

  // set vertex coordinates
  int ic = 0;
  for (int jv = 1; jv < NumVertices(); jv++) {
    SetVertexCoord(jv, ic, 2);
    ic++;
  };

  // generate edges in lexicographic order
  int ie = 0;
  for (int iv0 = 0; iv0 < NumVertices()-1; iv0++)
    for (int iv1 = iv0+1; iv1 < NumVertices(); iv1++) {
      SetEdge(ie, iv0, iv1);
      ie++;
    };

  if (ie != nume)
    throw PROCEDURE_ERROR(procname, "Programming error in edge generation.");

  // generate facets
  int num_vertices_per_facet = Dimension();
  for (int jf = 0; jf < numf; jf++) {

    for (int jv = 0; jv < num_vertices_per_facet; jv++) {
      int iv = (jf + jv) % numv;
      SetFacetVertex(jf, iv, true);
    };

    if (NumFacetVertices(jf) != num_vertices_per_facet) 
      throw PROCEDURE_ERROR
	(procname, "Programming error in facet generation.");
  };
}


void ISOSURFACE_TABLE_POLYHEDRON::GenPyramid(const int pyramid_dimension)
// generate a pyramid
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::GenPyramid";

  dimension = pyramid_dimension;

  if (!CheckDimension() || dimension < 2)
    throw PROCEDURE_ERROR(procname, "Illegal pyramid dimension.");

  int numv = (1L << (Dimension()-1)) + 1;
  int num_base_edges = (numv-1) * (Dimension()-1)/2;
  int nume = num_base_edges + (numv-1);
  int numf = 2*Dimension()-1;

  int apex = numv-1;

  SetSize(numv, nume, numf);

  // set vertex coordinates
  for (int iv = 0; iv < NumVertices()-1; iv++) {
    long mask = 1L;
    for (int ic = 0; ic < Dimension(); ic++) {
      int bit = mask & iv;
      int coord = 0;
      if (bit != 0)
	coord = 4;
      SetVertexCoord(iv, ic, coord);
      mask = mask << 1;
    };
  };

  // Set coordinate of pyramid
  for (int ic = 0; ic < Dimension(); ic++) {
    SetVertexCoord(apex, ic, 2);
  }

  // generate edges in lexicographic order
  int ie = 0;
  long control = 0;
  while (ie < num_base_edges) {

    // find first 0 bit in control
    int ic = 0;
    long mask = 1L;
    while ((mask & control) != 0) {
      mask = mask << 1;
      ic++;
    };

    // find start vertex by stripping ic bits from control
    int start = control;
    start = start >> ic;
    start = start << ic;

    for (int i = 0; i < (1L << ic); i++) {
      int iv0 = start + i;
      int iv1 = iv0 + (1L << ic);
      SetEdge(ie, iv0, iv1);
      ie++;
    };

    control++;
  };

  if (control+1 != (1L << (Dimension()-1)))
    throw PROCEDURE_ERROR(procname, "Programming error in edge generation.");

  for (int iv = 0; iv < NumVertices()-1; iv++) {
    SetEdge(num_base_edges+iv, iv, apex);
  }

  // generate facets containing apex
  int num_vertices_per_facet = 1+(numv-1)/2;
  long mask = 1L;
  for (int ic = 0; ic < Dimension()-1; ic++) {
    int jf0 = 2*ic;
    int jf1 = 2*ic+1;

    SetFacetVertex(jf0, apex, true);
    SetFacetVertex(jf1, apex, true);
    for (int iv = 0; iv < NumVertices()-1; iv++) {
      long bit = mask & iv;
      if (bit == 0) {
	SetFacetVertex(jf0, iv, true);
      }
      else {
	SetFacetVertex(jf1, iv, true);
      };
    };

    if (NumFacetVertices(jf0) != num_vertices_per_facet ||
	NumFacetVertices(jf1) != num_vertices_per_facet) {
      throw PROCEDURE_ERROR
	(procname, "Programming error in facet generation.");
    };

    mask = (mask << 1L);
  };

  // generate base facet
  for (int iv = 0; iv < NumVertices()-1; iv++) {
    SetFacetVertex(numf-1, iv, true);
  };

}

bool ISOSURFACE_TABLE_POLYHEDRON::CheckDimension() const
// check dimension
{
  if (dimension < 1)
    return(false);
  else
    return(true);
}

bool ISOSURFACE_TABLE_POLYHEDRON::Check(ERROR & error_msg) const
// check polyhedron
{
  if (!CheckDimension()) {
    error_msg.ClearAll();
    error_msg.AddMessage("Illegal polyhedron dimension ",
			 Dimension(), ".");
    return(false);
  }

  if (NumVertices() < 1) {
    error_msg.ClearAll();
    error_msg.AddMessage("Illegal number of vertices.");
    return(false);
  };

  if (NumEdges() < 1) {
    error_msg.ClearAll();
    error_msg.AddMessage("Illegal number of edges.");
    return(false);
  };

  if (vertex_coord == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for vertex coordinate list not allocated.");
    return(false);
  };

  if (edge_endpoint == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for edge endpoint list not allocated.");
    return(false);
  }

  for (int iv = 0; iv < NumVertices(); iv++) {
    for (int ic = 0; ic < Dimension(); ic++) {
      if ((VertexCoord(iv, ic) % 2) != 0) {
	error_msg.ClearAll();
	error_msg.AddMessage("Vertex coordinates must be even integers.");
	return(false);
      };
    };
  };

  for (int ie = 0; ie < NumEdges(); ie++) {
    for (int ip = 0; ip < 2; ip++) {
      int iv = EdgeEndpoint(ie, ip);
      if (iv < 0 || iv >= NumVertices()) {
	error_msg.ClearAll();
	error_msg.AddMessage("Illegal edge endpoint ", iv,
			     " for edge ", ie, ".");
	return(false);
      };
    };
  };

  if (NumFacets() > 0) {
    if (facet == NULL) {
      error_msg.ClearAll();
      error_msg.AddMessage("Memory for facet list not allocated.");
      return(false);
    };
  };

  return(true);
}

//**************************************************
// ISOSURFACE_VERTEX
//**************************************************

ISOSURFACE_VERTEX::ISOSURFACE_VERTEX()
  // constructor
{ 
  coord = NULL; 
  num_coord = 0; 
  is_label_set = false; 
}

ISOSURFACE_VERTEX::~ISOSURFACE_VERTEX()
  // destructor
{
  if (coord != NULL) { 
    delete [] coord;
    coord = NULL;
  };

  num_coord = 0;
  is_label_set = false;
}

void ISOSURFACE_VERTEX::SetNumCoord(const int numc)
{
  if (coord != NULL)
    delete [] coord;
  coord = NULL;

  num_coord = numc;
  coord = new COORD_TYPE[numc];
}

//**************************************************
// ISOSURFACE_TABLE
//**************************************************

ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::ISOSURFACE_TABLE_ENTRY()
// constructor
{
  num_simplices = 0;
  simplex_vertex_list = NULL;
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::~ISOSURFACE_TABLE_ENTRY()
// destructor
{
  FreeAll();
}

bool ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::Check
(ERROR & error_msg) const
{
  if (num_simplices < 0) {
    error_msg.ClearAll();
    error_msg.AddMessage
      ("Number of simplices in isosurface table entry must be non-negative.");
    return(false);
  }

  if (num_simplices > 0 && simplex_vertex_list == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for simplex vertex list not allocated.");
    return(false);
  }

  return(true);
}

void ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::FreeAll()
// free all memory
{
  delete [] simplex_vertex_list;
  simplex_vertex_list = NULL;
  num_simplices = 0;
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE(const int d):polyhedron(d)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(d, d-1);
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE
(const int dimension, const int simplex_dimension):polyhedron(dimension)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(dimension, simplex_dimension);
}

void ISOSURFACE_TABLE::Init(const int dimension, const int simplex_dimension)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  const char * procname = "ISOSURFACE_TABLE::Init";

  this->simplex_dimension = simplex_dimension;

  max_num_vertices = 20;
  // Note: Even tables for polyhedron of this size are probably impossible 
  //   to compute/store

  num_isosurface_vertices = 0;
  isosurface_vertex = NULL;

  encoding = NONSTANDARD;
  encoding_name = StandardEncodingName(NONSTANDARD);

  num_table_entries = 0;
  entry = NULL;
  is_table_allocated = false;
  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal polyhedron dimension.");
}

ISOSURFACE_TABLE::~ISOSURFACE_TABLE()
// destructor
{
  FreeAll();
}

std::string ISOSURFACE_TABLE::StandardEncodingName
(const ENCODING encoding)
{
  return(standard_encoding_name[encoding]);
}

void ISOSURFACE_TABLE::SetEncoding(const ENCODING encoding)
{
  this->encoding = encoding;
  encoding_name = StandardEncodingName(encoding);
}

void ISOSURFACE_TABLE::SetNonstandardEncoding(const std::string & name)
{
  this->encoding = NONSTANDARD;
  encoding_name = name;
}

void ISOSURFACE_TABLE::SetNumIsosurfaceVertices(const int num_vertices)
{
  if (isosurface_vertex != NULL)
    delete [] isosurface_vertex;
  isosurface_vertex = NULL;

  num_isosurface_vertices = num_vertices;
  isosurface_vertex = new ISOSURFACE_VERTEX[num_vertices];
}

void ISOSURFACE_TABLE::CheckIsoVerticesAlloc
(const char * procname, const int vstart, const int numv)
// check allocation of array isosurface_vertices
// procname = calling procedure name, for error messages
// vstart = first vertex
// numv = number of vertices required
{
  if (numv == 0) return;

  if (isosurface_vertex == NULL) {
    throw PROCEDURE_ERROR
      (procname, "Set number of isosurface vertices before storing vertices.");
  }

  if (numv+vstart > NumIsosurfaceVertices()) {
    throw PROCEDURE_ERROR
      (procname, "Illegal isosurface vertex index.");
  }
}

void ISOSURFACE_TABLE::StorePolyVerticesAsIsoVertices(const int vstart)
// store polyhedron vertices as isosurface vertices
// store polyhedron vertices starting at isosurface vertex vstart
{
  const int num_polyv = Polyhedron().NumVertices();
  const char * procname = "ISOSURFACE_TABLE::StorePolyVerticesAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polyv);

  for (int iv = 0; iv < num_polyv; iv++) {
    SetIsoVertexType(iv+vstart, ISOSURFACE_VERTEX::VERTEX);
    SetIsoVertexFace(iv+vstart, iv);
  }
}

void ISOSURFACE_TABLE::StorePolyEdgesAsIsoVertices(const int vstart)
// store polyhedron edges as isosurface vertices
// store polyhedron edges starting at isosurface vertex vstart
{
  const int num_polye = Polyhedron().NumEdges();
  const char * procname = "ISOSURFACE_TABLE::StorePolyEdgesAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polye);

  for (int ie = 0; ie < num_polye; ie++) {
    SetIsoVertexType(ie+vstart, ISOSURFACE_VERTEX::EDGE);
    SetIsoVertexFace(ie+vstart, ie);
  }
}

void ISOSURFACE_TABLE::StorePolyFacetsAsIsoVertices(const int vstart)
// store polyhedron facets as isosurface vertices
// store polyhedron facets starting at isosurface vertex vstart
{
  const int num_polyf = Polyhedron().NumFacets();
  const char * procname = "ISOSURFACE_TABLE::StorePolyFacetsAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polyf);

  for (int jf = 0; jf < num_polyf; jf++) {
    SetIsoVertexType(jf+vstart, ISOSURFACE_VERTEX::FACET);
    SetIsoVertexFace(jf+vstart, jf);
  }
}

void ISOSURFACE_TABLE::SetNumTableEntries(const int num_table_entries)
  // allocate table
{
  const char * procname = "ISOSURFACE_TABLE::SetNumTableEntries";

  if (entry != NULL) delete [] entry;
  entry = NULL;
  this->num_table_entries = 0;
  is_table_allocated = false;

  entry = new ISOSURFACE_TABLE_ENTRY[num_table_entries];
  if (entry == NULL && num_table_entries > 0)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for isosurface table.");

  this->num_table_entries = num_table_entries;
  is_table_allocated = true;
}


void ISOSURFACE_TABLE::SetNumSimplices(const TABLE_INDEX it, const int nums)
// set number of simplices in table entry it
// it = table entry
// nums = number of simplices
{
  const char * procname = "ISOSURFACE_TABLE::SetNumSimplices";

  if (!IsTableAllocated() || entry == NULL) {
    throw PROCEDURE_ERROR
      (procname, "Table must be allocated before entering table entries.");
  };

  if (it < 0 || it >= NumTableEntries())
    throw PROCEDURE_ERROR(procname, "Illegal table index.");
  if (nums < 0)
    throw PROCEDURE_ERROR
      (procname, "Number of simplices must be non-negative.");

  entry[it].num_simplices = 0;
  delete entry[it].simplex_vertex_list;
  entry[it].simplex_vertex_list = NULL;

  if (nums > 0)
    entry[it].simplex_vertex_list = 
      new ISOSURFACE_VERTEX_INDEX[nums*NumVerticesPerSimplex()];

  entry[it].num_simplices = nums;
}


void ISOSURFACE_TABLE::SetSimplexVertex
  (const TABLE_INDEX it, const int is, const int k, 
   const ISOSURFACE_VERTEX_INDEX isov)
// set simplex vertex
// it = index table entry.  In range [0..NumTableEntries()-1].
// is = index simplex.  
// k = k'th simplex vertex.  In range [0..NumVerticesPerSimplex()-1].
// isov = index of isosurface vertex
{
  entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+k] = isov;
}

bool ISOSURFACE_TABLE::CheckDimension(const int d) const
// check dimension
{
  if (d < 1)
    return(false);
  else
    return(true);
}

bool ISOSURFACE_TABLE::CheckTable(ERROR & error_msg) const
// check table
{
  if (polyhedron.NumVertices() >= LONG_BIT) {
    error_msg.ClearAll();
    error_msg.AddMessage("Too many polyhedron vertices");
    return(false);
  }

  if (polyhedron.NumVertices() > MaxNumVertices()) {
    error_msg.ClearAll();
    error_msg.AddMessage("Too many polyhedron vertices");
    return(false);
  }

  if (polyhedron.NumVertices() < 1) {
    error_msg.ClearAll();
    error_msg.AddMessage("Polyhedron must have at least one vertex.");
    return(false);
  }

  if (entry == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for isosurface table not allocated.");
    return(false);
  }

  for (int it = 0; it < NumTableEntries(); it++)
    if (!entry[it].Check(error_msg)) {
      error_msg.AddMessage("Error detected at isosurface table entry ", 
			   it, ".");
      return(false);
    }

  for (int jt = 0; jt < NumTableEntries(); jt++)
    for (int is = 0; is < NumSimplices(jt); is++)
      for (int iv = 0; iv < NumVerticesPerSimplex(); iv++) {
	int iso_v = int(SimplexVertex(jt, is, iv));
	if (iso_v < 0 || iso_v >= NumIsosurfaceVertices()) {
	  error_msg.ClearAll();
	  error_msg.AddMessage
	    ("Illegal isosurface vertex ", iso_v, " in isosurface table entry ", jt, "."); 
	  return(false);
	}
      };

  return(true);
}

bool ISOSURFACE_TABLE::Check(ERROR & error_msg) const
{
  if (!Polyhedron().Check(error_msg)) return(false);
  if (!CheckTable(error_msg)) return(false);
  return(true);
}

void ISOSURFACE_TABLE::FreeAll()
// free all memory
{
  if (entry != NULL) {
    for (int i = 0; i < num_table_entries; i++)
      entry[i].FreeAll();
    delete [] entry;
    entry = NULL;
  };
  num_table_entries = 0;
  is_table_allocated = false;

  polyhedron.FreeAll();

  delete [] isosurface_vertex;
  isosurface_vertex = NULL;
  num_isosurface_vertices = 0;
}


//**************************************************
// ISOSURFACE_EDGE_TABLE
//**************************************************

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::ISOSURFACE_EDGE_TABLE_ENTRY()
// constructor
{
  num_edges = 0;
  edge_endpoint_list = NULL;
}

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::~ISOSURFACE_EDGE_TABLE_ENTRY()
// destructor
{
  FreeAll();
}

bool ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::Check
(ERROR & error_msg) const
{
  if (num_edges < 0) {
    error_msg.ClearAll();
    error_msg.AddMessage
      ("Number of edges in isosurface table entry must be non-negative.");
    return(false);
  }

  if (num_edges > 0 && edge_endpoint_list == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for edge endpoint list not allocated.");
    return(false);
  }

  return(true);
}

void ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::FreeAll()
// free all memory
{
  delete [] edge_endpoint_list;
  edge_endpoint_list = NULL;
  num_edges = 0;
}

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE(const int d):ISOSURFACE_TABLE(d)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(d);
}

void ISOSURFACE_EDGE_TABLE::Init(const int d)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  edge_entry = NULL;
}

ISOSURFACE_EDGE_TABLE::~ISOSURFACE_EDGE_TABLE()
// destructor
{
  FreeAll();
}

void ISOSURFACE_EDGE_TABLE::SetNumTableEntries(const int num_table_entries)
  // allocate table
{
  const char * procname = "ISOSURFACE_EDGE_TABLE::SetNumTableEntries";

  ISOSURFACE_TABLE::SetNumTableEntries(num_table_entries);

  edge_entry = new ISOSURFACE_EDGE_TABLE_ENTRY[num_table_entries];
  if (edge_entry == NULL)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for isosurface edge table.");
}

bool ISOSURFACE_EDGE_TABLE::CheckTable(ERROR & error_msg) const
// check table
{
  if (!ISOSURFACE_TABLE::CheckTable(error_msg)) return(false);

  if (edge_entry == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for isosurface table edge entries not allocated.");
    return(false);
  }

  for (int it = 0; it < NumTableEntries(); it++)
    if (!edge_entry[it].Check(error_msg)) {
      error_msg.AddMessage("Error detected at isosurface table entry ", 
			   it, ".");
      return(false);
    }

  for (int jt = 0; jt < NumTableEntries(); jt++)
    for (int ie = 0; ie < NumEdges(jt); ie++)
      for (int iend = 0; iend < 2; iend++) {
	int iso_v = int(EdgeEndpoint(jt, ie, iend));
	if (iso_v < 0 || iso_v >= NumIsosurfaceVertices()) {
	  error_msg.ClearAll();
	  error_msg.AddMessage
	    ("Illegal isosurface vertex ", iso_v, " in isosurface table edge entry ", jt, "."); 
	  return(false);
	}
      };

  return(true);
}

bool ISOSURFACE_EDGE_TABLE::Check(ERROR & error_msg) const
{
  if (!Polyhedron().Check(error_msg)) return(false);
  if (!CheckTable(error_msg)) return(false);
  return(true);
}

void ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE::FreeAll()
// free all memory
{
  if (edge_entry != NULL) {
    for (int i = 0; i < num_table_entries; i++) {
      edge_entry[i].FreeAll();
    }
    delete [] edge_entry;
    edge_entry = NULL;  
  };

  ISOSURFACE_TABLE::FreeAll();
}

namespace {
  class EDGE_CONTAINER { 
  public:
    EDGE_INDEX v[2]; 

    EDGE_CONTAINER(){};
    bool operator < (const EDGE_CONTAINER & e) const {
      if (v[0] < e.v[0]) { return(true); }
      else if (v[0] == e.v[0] && v[1] < e.v[1]) { return(true); }
      else {return(false); }
    }
  };

}

void ISOSURFACE_EDGE_TABLE::GenEdgeLists()
// for each table entry, generate edge lists from simplex lists
{
  PROCEDURE_ERROR error("ISOSURFACE_EDGE_TABLE::GenEdgeLists");
  EDGE_INDEX vlist[Dimension()];
  typedef set<EDGE_CONTAINER> EDGE_SET;
  EDGE_SET eset;

  if (!IsTableAllocated()) {
    error.AddMessage("Programming error: Isosurface table not allocated.");
    throw error;
  };

  for (int it = 0; it < NumTableEntries(); it++) {
    eset.clear();

    for (int is = 0; is < NumSimplices(it); is++) {

      // get simplex vertices
      for (int iv = 0; iv < Dimension(); iv++) {
	vlist[iv] = SimplexVertex(it, is, iv);
      }
      sort(vlist, vlist+Dimension());

      // store simplex edges
      for (int i0 = 0; i0 < Dimension(); i0++)
	for (int i1 = i0+1; i1 < Dimension(); i1++) {
	  EDGE_CONTAINER e;
	  e.v[0] = vlist[i0];
	  e.v[1] = vlist[i1];

	  eset.insert(e);
	}
    }

    edge_entry[it].FreeAll();
    if (eset.size() > 0) {
      edge_entry[it].edge_endpoint_list = new EDGE_INDEX[2*eset.size()];
      if (edge_entry[it].edge_endpoint_list == NULL) {
	error.AddMessage("Unable to allocate memory for edge list in table entry ", it, ".");
	throw error;
      };
      edge_entry[it].num_edges = eset.size();

      int ie = 0;
      for (EDGE_SET::iterator edge_iter = eset.begin(); 
	   edge_iter != eset.end(); 
	   edge_iter++) {
	edge_entry[it].edge_endpoint_list[2*ie] = edge_iter->v[0];
	edge_entry[it].edge_endpoint_list[2*ie+1] = edge_iter->v[1];
	ie++;
      }
    }
  }

}



//**************************************************
// Routines for generating polyhedra
//**************************************************

void IJKTABLE::generate_prism
  (const ISOSURFACE_TABLE_POLYHEDRON & base_polyhedron,
   ISOSURFACE_TABLE_POLYHEDRON & prism)
// generate a prism with base base_polyhedron
// first numv vertices have last coordinate = 0
// last numv vertices have last coordinate = 2
// first nume edges connect first numv vertices
// second nume edges connect second numv vertices
// third numv edges connect first to second set of vertices
// (numv = # vertices in base_polyhedron; nume = # edges in base_polyhedron)
{
  int dim = base_polyhedron.Dimension();
  int numc = dim;
  int numv = base_polyhedron.NumVertices();
  int nume = base_polyhedron.NumEdges();
  int numf = base_polyhedron.NumFacets();
  int prism_dim = dim + 1;
  int prism_numc = prism_dim;
  int prism_lastc = prism_numc - 1;
  int prism_numv = numv * 2;
  int prism_nume = nume * 2 + numv;
  int prism_numf = 2 + numf;
  prism.SetDimension(prism_dim);
  prism.SetSize(prism_numv, prism_nume, prism_numf);

  // set prism vertex coord
  for (int iv = 0; iv < numv; iv++) {
    for (int ic = 0; ic < prism_lastc; ic++) {
      int coord = base_polyhedron.VertexCoord(iv, ic);
      prism.SetVertexCoord(iv, ic, coord);
      prism.SetVertexCoord(iv+numv, ic, coord);
    };
    prism.SetVertexCoord(iv, prism_lastc, 0);
    prism.SetVertexCoord(iv+numv, prism_lastc, 2);
  };

  // set edges
  for (int ie = 0; ie < base_polyhedron.NumEdges(); ie++) {
    int iv0 = base_polyhedron.EdgeEndpoint(ie, 0);
    int iv1 = base_polyhedron.EdgeEndpoint(ie, 1);
    prism.SetEdge(ie, iv0, iv1);
    prism.SetEdge(ie+nume, iv0+numv, iv1+numv);
  };

  for (int iv = 0; iv < base_polyhedron.NumVertices(); iv++) {
    int ie = 2*nume + iv;
    prism.SetEdge(ie, iv, iv+numv);
  };

  // set facets
  for (int iv = 0; iv < numv; iv++) {
    // set two base facets
    prism.SetFacetVertex(0, iv, true);
    prism.SetFacetVertex(1, iv+numv, true);
  };

  for (int jf = 0; jf < numf; jf++) {
    // set prism facet corresponding to original facet jf
    int prism_jf = 2 + jf;
    for (int iv = 0; iv < numv; iv++) {
      if (base_polyhedron.FacetVertexFlag(jf, iv)) {
	prism.SetFacetVertex(prism_jf, iv, true);
	prism.SetFacetVertex(prism_jf, iv+numv, true);
      };
    };
  };


}

//**************************************************
// UTILITY FUNCTIONS
//**************************************************

u_long IJKTABLE::calculate_num_entries
(const int num_vert, const int num_colors)
  // calculate num table entries = (num_colors)^(num_vert)
{
  const char * procname = "calculate_num_entries";

  u_long num_table_entries = 0;

  if (num_colors < 1)
    throw PROCEDURE_ERROR(procname, "Number of colors must be positive.");

  const u_long max2 = ULONG_MAX/num_colors;

  num_table_entries = 1;
  for (int iv = 0; iv < num_vert; iv++) {
    if (num_table_entries > max2)
      throw PROCEDURE_ERROR(procname, "Number of entries is too large.");

    num_table_entries = num_table_entries * num_colors;
  };

  return(num_table_entries);
}

void IJKTABLE::convert2base
(const u_long ival, const u_int base, int * digit, 
 const u_int max_num_digits)
{
  u_long jval = ival;
  for (int i = 0; i < max_num_digits; i++) {
    digit[i] = jval % base;
    jval = jval/base;
  };

  if (jval != 0) {
    PROCEDURE_ERROR error("convert2base");
    error.AddMessage("Error converting ", ival, " to base ", base, ".");
    error.AddMessage("Output has more than ", max_num_digits, " digits.");

    throw error;
  };
}

