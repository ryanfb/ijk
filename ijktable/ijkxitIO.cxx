/// \file ijkxitIO.cxx
/// I/O routines for .xit file
/// .xit is an xml format 
/// Version 0.1.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2007, 2006 Rephael Wenger

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

#include <iostream>
#include <sstream>
#include <ctime>

#include "expat.h"

#include <string.h>

#include "ijkxitIO.h"

using namespace std;
using namespace IJK;
using namespace IJKXIO;

//**************************************************
// read_xit
//**************************************************

namespace {
  // global error class
  ERROR read_error;
}

// types
typedef enum { START_TAG, END_TAG } XML_TAG_TYPE;
typedef enum { NO_ELEMENT,
	       ISOTABLE, VERSION, CREATION_DATE, DIMENSION, 
	       POLY, VERTICES, NUM_VERTICES, C,
	       EDGES, NUM_EDGES, V, W,
	       FACETS, NUM_FACETS, F,
	       ISO_VERTICES, IN_E, IN_V, IN_F, L,
	       TABLE, ENCODING, NUM_ENTRIES, S,
	       UNKNOWN_ELEMENT } XML_ELEMENT;

typedef struct {
  XML_TAG_TYPE tag_type;
  XML_ELEMENT element;
  bool flag_eof;
  /* NOT YET USED
     string char_data; */
} USER_DATA;

const char * xml_strings[] = 
  { "No element", "isotable", "version",
    "creationDate", "dimension",
    "poly", "vertices", "numVertices", "c",
    "edges", "numEdges", "v", "w",
    "facets", "numFacets", "f",
    "isoVertices", "inE", "inV", "inF", "L",
    "table", "encoding", "numEntries", "s",
    "Unknown element" };

// local functions
namespace {
  void read_tag(istream & in, XML_Parser parser, USER_DATA & user_data);
  void read_start_tag(istream & in, XML_Parser parser, 
		      const XML_ELEMENT element, USER_DATA & user_data);
  void read_end_tag(istream & in, XML_Parser parser, 
		    const XML_ELEMENT element, USER_DATA & user_data);
  std::string read_string(istream & in, XML_Parser parser, 
			  USER_DATA & user_data);
  void read_isosurface_vertex
  (istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKTABLE::ISOSURFACE_TABLE & table, const int isov);
  void XMLCALL startElement(void *data, const char *name, 
			    const char **attr);
  void XMLCALL endElement(void *data, const char *ename);
  void XMLCALL startPolyElement(void *data, const char *name, 
				const char **attr);
  void XMLCALL startIsoVerticesElement(void *data, const char *name,
				       const char **attr);
  void XMLCALL endIsoVerticesElement(void *data, const char *ename);
  void XMLCALL startTableElement(void *data, const char *name, 
				 const char **attr);
  void initUserData(USER_DATA & user_data);
}


void IJKXIO::read_xit(istream & in, IJKTABLE::ISOSURFACE_TABLE & table)
{
  XML_Parser parser = NULL;
  USER_DATA user_data;
  int polyDimension = 0;
  int surfaceDimension = 0;
  int numVertices = 0;
  int numEdges = 0;
  int numFacets = 0;
  int numIsosurfaceVertices = 0;
  int numEntries = 0;
  int d = 0;
  ostringstream msg;

  read_error.AddProcMessage("read_xit");

  try {

  initUserData(user_data);

  parser = XML_ParserCreate(NULL);

  if (!parser) 
    throw read_error("Couldn't allocate memory for XML parser.");

  XML_SetElementHandler(parser, startElement, endElement);
  XML_SetUserData(parser, (void *) &user_data);

  // read in version, creationDate
  read_start_tag(in, parser, VERSION, user_data);
  read_start_tag(in, parser, CREATION_DATE, user_data);

  // read dimension
  read_start_tag(in, parser, DIMENSION, user_data);
  in >> polyDimension;
  in >> surfaceDimension;
  if (in.fail()) throw read_error("Error reading dimension.");
  table.SetDimension(polyDimension);
  table.SetSimplexDimension(surfaceDimension);

  // read in polytope description
  XML_SetElementHandler(parser, startPolyElement, endElement);
  read_start_tag(in, parser, POLY, user_data);

  // read vertices
  read_start_tag(in, parser, VERTICES, user_data);
  read_start_tag(in, parser, NUM_VERTICES, user_data);
  in >> numVertices;
  if (in.fail()) throw read_error("Error reading number of vertices.");
  table.SetNumPolyVertices(numVertices);

  for (int i = 0; i < numVertices; i++) {
    read_start_tag(in, parser, C, user_data);
    if (in.fail()) {
      read_error.AddMessage
	("Error reading vertex coordinates of polyhedron vertex ", i, ".");
      throw read_error;
    };
    
    int coord;
    for (int d = 0; d < polyDimension; d++) {
      in >> coord;
      if (in.fail()) {
	read_error.AddMessage
	  ("Error reading polyhedron vertex ", i, ", coordinate ", d, ".");
	throw read_error;
      }
      table.SetPolyVertexCoord(i, d, coord);
    }
  }

  // read edges
  read_start_tag(in, parser, EDGES, user_data);
  read_start_tag(in, parser, NUM_EDGES, user_data);
  in >> numEdges;
  if (in.fail()) throw read_error("Error reading number of edges.");
  table.SetNumPolyEdges(numEdges);

  for (int i = 0; i < numEdges; i++) {
    read_start_tag(in, parser, V, user_data);
    int end0, end1;
    in >> end0;
    in >> end1;
    if (in.fail()) {
      read_error.AddMessage("Error reading polyhedron edge ", i, ".");
      throw read_error;
    };
    if (end0 < 0 || end1 < 0 || end0 >= numVertices || end1 >= numVertices) {
      msg.str("");
      msg << "Error reading polyhedron edge " << i 
	  << ".  Illegal vertex indices: " << end0 << " " << end1 << "."
	  << ends;
      throw read_error(msg.str());
    };
    table.SetPolyEdge(i, end0, end1);
  }

  // read facets
  read_start_tag(in, parser, FACETS, user_data);
  read_start_tag(in, parser, NUM_FACETS, user_data);
  in >> numFacets;
  if (in.fail()) throw read_error("Error reading number of facets.");
  table.SetNumPolyFacets(numFacets);

  for (int i = 0; i < numFacets; i++) {
    read_start_tag(in, parser, F, user_data);
    int numv;
    in >> numv;
    if (in.fail()) {
      read_error.AddMessage("Error reading polyhedron facet ", i, ".");
      throw read_error;
    }

    for (int j = 0; j < numv; j++) {
      int v;
      in >> v;
      if (in.fail()) {
	msg.str("");
	msg << "Error reading polyhedron facet " << i << ", vertex " 
	    << j << ends;
	throw read_error(msg.str());
      }
      if (v < 0 || v >= numVertices) {
	msg.str("");
	msg << "Error reading polyhedron facet " << i << ", vertex " 
	    << j << ".  Illegal vertex index: " << v << "." << ends;
	throw read_error(msg.str());
      }
      table.SetPolyFacetVertex(i, v, true);
    }
  }

  // read isosurface vertices
  XML_SetElementHandler(parser, 
			startIsoVerticesElement, endIsoVerticesElement);
  read_start_tag(in, parser, ISO_VERTICES, user_data);
  read_start_tag(in, parser, NUM_VERTICES, user_data);
  in >> numIsosurfaceVertices;
  table.SetNumIsosurfaceVertices(numIsosurfaceVertices);
  for (int iv = 0; iv < numIsosurfaceVertices; iv++) {
    read_isosurface_vertex(in, parser, user_data, table, iv);
  }

  // read isosurface lookup table
  XML_SetElementHandler(parser, startTableElement, endElement);
  read_start_tag(in, parser, TABLE, user_data);

  // read encoding
  read_start_tag(in, parser, ENCODING, user_data);
  std::string encoding_name = read_string(in, parser, user_data);
  if (encoding_name == table.StandardEncodingName(IJKTABLE::ISOSURFACE_TABLE::BINARY)) {
    table.SetEncoding(IJKTABLE::ISOSURFACE_TABLE::BINARY);
  }
  else if (encoding_name == table.StandardEncodingName(IJKTABLE::ISOSURFACE_TABLE::BASE3)) {
    table.SetEncoding(IJKTABLE::ISOSURFACE_TABLE::BASE3);
  }
  else {
    table.SetNonstandardEncoding(encoding_name);
  }
  
  read_start_tag(in, parser, NUM_ENTRIES, user_data);
  in >> numEntries;
  if (in.fail())
    throw read_error("Error reading number of lookup table entries.");
  table.SetNumTableEntries(numEntries);

  if (table.NumTableEntries() != numEntries)
    throw read_error("Error in number of table entries.");

  for (int it = 0; it < numEntries; it++) {
    read_start_tag(in, parser, S, user_data);
    int nums;
    in >> nums;
    if (in.fail()) {
      read_error.AddMessage("Error reading table entry ", it, ".");
      throw read_error;
    }

    table.SetNumSimplices(it, nums);

    for (int js = 0; js < nums; js++) {
      for (int k = 0; k <= surfaceDimension; k++) {
	int ie;
	in >> ie;
	if (in.fail()) {
	  msg.str("");
	  msg << "Error reading table entry " << it
	      << ", simplex " << js << ", "
	      << " vertex/edge " << k << ends;
	  throw read_error(msg.str());
	}

	// NOTE: There is no check on the value of ie

	table.SetSimplexVertex(it, js, k, ie);
      }
    };
  }
  } catch(...) {
    if (parser != NULL)
      XML_ParserFree(parser);

    // re-throw the exception
    throw;
  }

  XML_ParserFree(parser);
}


namespace {

void read_tag
(istream & in, XML_Parser parser, USER_DATA & user_data)
{
  const int BUFFSIZE = 1024;
  static char buffer[BUFFSIZE];

  int j;
  char s;
  int done = 0;
  int length = 0;
  while (length < BUFFSIZE && in.get(s) && s != '\n' && s != '>') {
    buffer[length] = s;
    length++;
  };
  if (length < BUFFSIZE && s == '>' && in.good()) {
    buffer[length] = s;
    length++;
  };

  if (in.eof()) {
    user_data.flag_eof = 1;
    done = 1;
  }

  if (XML_Parse(parser, buffer, length, done) == XML_STATUS_ERROR) {
    string msg = "Parse error " + 
      string(XML_ErrorString(XML_GetErrorCode(parser)));
    throw read_error(msg);
  }
}

void read_start_tag(istream & in, XML_Parser parser, 
		    const XML_ELEMENT element, USER_DATA & user_data)
{
  user_data.element = NO_ELEMENT;
  while(in.good() &&
	(user_data.element != element || user_data.tag_type != START_TAG)) {
    read_tag(in, parser, user_data);
  }

  if (in.eof() || in.fail()) {
    string msg = "Missing start tag: <" + string(xml_strings[element]) + ">";
    throw read_error(msg);
  };
}

void read_end_tag(istream & in, XML_Parser parser, 
		  const XML_ELEMENT element, USER_DATA & user_data)
{
  user_data.element = NO_ELEMENT;
  while(in.good() &&
	(user_data.element != element || user_data.tag_type != END_TAG)) {
    read_tag(in, parser, user_data);
  }

  if (in.eof() || in.fail()) {
    string msg = "Missing end tag: <" + string(xml_strings[element]) + ">";
    throw read_error(msg);
  };
}

std::string read_string(istream & in, XML_Parser parser, 
			USER_DATA & user_data)
// read string, starting with first non-blank character and ending
// with next blank character
{
  string s;
  char c;

  // read white space before string
  in.get(c);
  while (in.good() && isspace(c)) { in.get(c); }

  // read string
  while (in.good() && !isspace(c) && c != '<') {
    s = s + c;
    in.get(c);
  }
  if (in.good() && c == '<')
    in.putback(c);

  return(s);
}

void read_isosurface_vertex
(istream & in, XML_Parser parser, USER_DATA & user_data,
 IJKTABLE::ISOSURFACE_TABLE & table, const int isov)
{
  user_data.element = NO_ELEMENT;
  read_start_tag(in, parser, W, user_data);

  while (in.good() && 
	 ((user_data.element != IN_V && user_data.element != IN_E &&
	   user_data.element != IN_F && user_data.element != C) ||
	  user_data.tag_type != START_TAG)) {
    read_tag(in, parser, user_data);
  }

  switch(user_data.element) {
  case IN_V:
    int iv;
    in >> iv;
    table.SetIsoVertexType(isov, IJKTABLE::ISOSURFACE_VERTEX::VERTEX);
    table.SetIsoVertexFace(isov, iv);
    read_end_tag(in, parser, IN_V, user_data);
    break;

  case IN_E:
    int ie;
    in >> ie;
    table.SetIsoVertexType(isov, IJKTABLE::ISOSURFACE_VERTEX::EDGE);
    table.SetIsoVertexFace(isov, ie);
    read_end_tag(in, parser, IN_E, user_data);
    break;

  case IN_F:
    int jf;
    in >> jf;
    table.SetIsoVertexType(isov, IJKTABLE::ISOSURFACE_VERTEX::FACET);
    table.SetIsoVertexFace(isov, jf);
    read_end_tag(in, parser, IN_F, user_data);
    break;

  case C:
    int num_coord = table.Dimension();
    table.SetIsoVertexNumCoord(isov, num_coord);
    for (int d = 0; d < num_coord; d++) {
      IJKTABLE::ISOSURFACE_VERTEX::COORD_TYPE coord;
      in >> coord;
      table.SetIsoVertexCoord(isov, d, coord);
    }
    read_end_tag(in, parser, C, user_data);
    break;
  }

  read_tag(in, parser, user_data);
  while (in.good() && 
	 ((user_data.element != L || user_data.tag_type != START_TAG) &&
	  (user_data.element != W || user_data.tag_type != END_TAG))) {
    read_tag(in, parser, user_data);
  }

  if (in.eof() || in.fail()) {
    string msg = "Missing end tag for isosurface vertex tag <w>.";
    throw read_error(msg);
  };

  if (user_data.element == L && user_data.tag_type == START_TAG) {
    std::string label = read_string(in, parser, user_data);
    table.SetIsoVertexLabel(isov, label);
  }
}

void XMLCALL startElement(void *data, const char *name, const char **attr)
{
  USER_DATA * user_data = (USER_DATA *) data;

  if (user_data == NULL) return;

  user_data->tag_type = START_TAG;

  if (strcmp(name, xml_strings[VERSION]) == 0) {
    user_data->element = VERSION;
  }
  else if (strcmp(name, xml_strings[CREATION_DATE]) == 0) {
    user_data->element = CREATION_DATE;
  }
  else if (strcmp(name, xml_strings[DIMENSION]) == 0) {
    user_data->element = DIMENSION;
  }
  else {
    user_data->element = UNKNOWN_ELEMENT;
  };
}

void XMLCALL endElement(void *data, const char *ename)
{
  USER_DATA * user_data = (USER_DATA *) data;

  if (user_data == NULL) return;

  user_data->tag_type = END_TAG;
}

/* NOT YET USED
static void XMLCALL
charDataHandler(void *data, const XML_Char *s, int length)
{
  USER_DATA * user_data = (USER_DATA *) data;

  if (user_data == NULL) return;

  user_data->char_data += s;
}
*/

void XMLCALL startPolyElement(void *data, const char *name, 
			      const char **attr)
{
  USER_DATA * user_data = (USER_DATA *) data;

  if (user_data == NULL) return;

  user_data->tag_type = START_TAG;

  if (strcmp(name, xml_strings[VERTICES]) == 0) {
    user_data->element = VERTICES;
  }
  else if (strcmp(name, xml_strings[NUM_VERTICES]) == 0) {
    user_data->element = NUM_VERTICES;
  }
  else if (strcmp(name, xml_strings[C]) == 0) {
    user_data->element = C;
  }
  else if (strcmp(name, xml_strings[EDGES]) == 0) {
    user_data->element = EDGES;
  }
  else if (strcmp(name, xml_strings[NUM_EDGES]) == 0) {
    user_data->element = NUM_EDGES;
  }
  else if (strcmp(name, xml_strings[V]) == 0) {
    user_data->element = V;
  }
  else if (strcmp(name, xml_strings[FACETS]) == 0) {
    user_data->element = FACETS;
  }
  else if (strcmp(name, xml_strings[NUM_FACETS]) == 0) {
    user_data->element = NUM_FACETS;
  }
  else if (strcmp(name, xml_strings[F]) == 0) {
    user_data->element = F;
  }
  else if (strcmp(name, xml_strings[POLY]) == 0) {
    user_data->element = POLY;
  }
  else {
    user_data->element = UNKNOWN_ELEMENT;
  };

}

void setIsoVerticesElement(void *data, const char *name)
{
  USER_DATA * user_data = (USER_DATA *) data;

 if (strcmp(name, xml_strings[ISO_VERTICES]) == 0) {
    user_data->element = ISO_VERTICES;
  }
  else if (strcmp(name, xml_strings[NUM_VERTICES]) == 0) {
    user_data->element = NUM_VERTICES;
  }
  else if (strcmp(name, xml_strings[W]) == 0) {
    user_data->element = W;
  }
  else if (strcmp(name, xml_strings[IN_V]) == 0) {
    user_data->element = IN_V;
  }
  else if (strcmp(name, xml_strings[IN_E]) == 0) {
    user_data->element = IN_E;
  }
  else if (strcmp(name, xml_strings[IN_F]) == 0) {
    user_data->element = IN_F;
  }  
  else if (strcmp(name, xml_strings[C]) == 0) {
    user_data->element = C;
  }
  else if (strcmp(name, xml_strings[L]) == 0) {
    user_data->element = L;
  }
  else {
    user_data->element = UNKNOWN_ELEMENT;
  };
}


void XMLCALL startIsoVerticesElement(void *data, const char *name,
				     const char **attr)
{
  USER_DATA * user_data = (USER_DATA *) data;

  if (user_data == NULL) return;

  user_data->tag_type = START_TAG;

  setIsoVerticesElement(data, name);
}

void XMLCALL endIsoVerticesElement(void *data, const char *name)
{
  USER_DATA * user_data = (USER_DATA *) data;

  if (user_data == NULL) return;

  user_data->tag_type = END_TAG;

  setIsoVerticesElement(data, name);
}
						      

void XMLCALL startTableElement(void *data, const char *name, 
			       const char **attr)
{
  USER_DATA * user_data = (USER_DATA *) data;

  if (user_data == NULL) return;

  user_data->tag_type = START_TAG;
  if (strcmp(name, xml_strings[TABLE]) == 0) {
    user_data->element = TABLE;
  }
  else if (strcmp(name, xml_strings[ENCODING]) == 0) {
    user_data->element = ENCODING;
  }
  else if (strcmp(name, xml_strings[NUM_ENTRIES]) == 0) {
    user_data->element = NUM_ENTRIES;
  }
  else if (strcmp(name, xml_strings[S]) == 0) {
    user_data->element = S;
  }
  else {
    user_data->element = UNKNOWN_ELEMENT;
  };
}

void initUserData(USER_DATA & user_data)
{
  user_data.tag_type = START_TAG;
  user_data.element = NO_ELEMENT;
  user_data.flag_eof = false;
}

}

//**************************************************
// write_xit
//**************************************************

// local functions
namespace {
  void write_date(ostream & out);
}

void IJKXIO::write_xit
(ostream & out, const IJKTABLE::ISOSURFACE_TABLE & table)
{
  int dimension = table.Dimension();
  int numv = table.Polyhedron().NumVertices();
  int nume = table.Polyhedron().NumEdges();
  int numf = table.Polyhedron().NumFacets();
  int numisov = table.NumIsosurfaceVertices();
  int numTableEntries = table.NumTableEntries();

  out << "<?xml version=\"1.0\"?>" << endl;
  out << "<isotable>" << endl;
  out << "<!-- Isosurface lookup table -->" << endl;
  out << "<version> 1.0 </version>" << endl;
  out << "<creationDate> ";
  write_date(out);
  out << " </creationDate>" << endl;
  out << "<dimension> " << table.Dimension() << "  "
      << table.SimplexDimension() << " </dimension>" << endl;
  out << "<poly>" << endl;
  out << "<vertices>" << endl;
  out << "<numVertices> " << numv << " </numVertices>" << endl;
  for (int i = 0; i < numv; i++) {
    out << "<c> ";
    for (int d = 0; d < dimension; d++) {
      out << table.Polyhedron().VertexCoord(i, d) << " ";
    };
    out << "</c>" << endl;
  }
  out << "</vertices>" << endl;
  out << "<edges>" << endl;
  out << "<numEdges> " << nume << " </numEdges>" << endl;
  for (int i = 0; i < nume; i++) {
    out << "<v> ";
    for (int j = 0; j < 2; j++) {
      out << table.Polyhedron().EdgeEndpoint(i, j) << " ";
    };
    out << "</v>" << endl;
  };
  out << "</edges>" << endl;
  out << "<facets>" << endl;
  out << "<numFacets> " << numf << " </numFacets>" << endl;
  for (int i = 0; i < numf; i++) {
    out << "<f> ";
    out << table.Polyhedron().NumFacetVertices(i) << " ";
    for (int jv = 0; jv < numv; jv++) {
      if (table.Polyhedron().FacetVertexFlag(i, jv)) {
	out << jv << " ";
      };
    };
    out << "</f>" << endl;
  };
  out << "</facets>" << endl;
  out << "</poly>" << endl;

  out << "<isoVertices>" << endl;
  out << "<numVertices> " << numisov << " </numVertices>" << endl;
  for (int i = 0; i < numisov; i++) {
    out << "<w> ";
    IJKTABLE::ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE vtype =
      table.IsosurfaceVertex(i).Type();
    int iv, ie, jf, numc;
    switch(vtype) {
    case IJKTABLE::ISOSURFACE_VERTEX::VERTEX:
      iv = table.IsosurfaceVertex(i).Face();
      out << "<inV> " << iv << " </inV> ";
      break;

    case IJKTABLE::ISOSURFACE_VERTEX::EDGE:
      ie = table.IsosurfaceVertex(i).Face();
      out << "<inE> " << ie << " </inE> ";
      break;

    case IJKTABLE::ISOSURFACE_VERTEX::FACET:
      jf = table.IsosurfaceVertex(i).Face();
      out << "<inF> " << jf << " </inF> ";
      break;

    case IJKTABLE::ISOSURFACE_VERTEX::POINT:
      numc = table.IsosurfaceVertex(i).NumCoord();
      out << "<c> ";
      if (numc > 0) {
	out << "<c> ";
	for (int ic = 0; ic < numc; ic++) {
	  out << table.IsosurfaceVertex(i).Coord(ic) << " ";
	}
	out  << "</c> ";
      }
      break;
    }

    if (table.IsosurfaceVertex(i).IsLabelSet()) {
      out << "<L> " << table.IsosurfaceVertex(i).Label() << " </L> ";
    }
    out << "</w>" << endl;
  }
  out << "</isoVertices>" << endl;

  out << "<table>" << endl;
  out << "<encoding> " << table.EncodingName() << " </encoding>" << endl;
  out << "<numEntries> " << numTableEntries << " </numEntries>" 
      << endl;
  for (int it = 0; it < numTableEntries; it++) {
    int nums = table.NumSimplices(it);
    out << "<s> ";
    out << nums << " ";
    for (int js = 0; js < nums; js++) {
      for (int kv = 0; kv < table.SimplexDimension()+1; kv++) {
	out << int(table.SimplexVertex(it, js, kv)) << " ";
      }
    }
    out << "</s>" << endl;
  };
  out << "</table>" << endl;
  out << "</isotable>" << endl;
}

namespace {

void write_date(ostream & out)
{
  const size_t BUFFER_LENGTH = 16;
  char buffer[BUFFER_LENGTH];
  time_t rawtime;
  struct tm * timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, BUFFER_LENGTH, "%Y-%m-%d", timeinfo);
  out << buffer;
}

}

