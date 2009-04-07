/// \file ijkgentable.cxx
/// generate isosurface table
/// Version 0.3.0

/*
  IJK: Isosurface Jeneration Code
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

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "ijktable.h"
#include "ijkxitIO.h"

using namespace std;
using namespace IJK;
using namespace IJKTABLE;

// types
enum MESH_POLYHEDRON_TYPE { SIMPLEX, CUBE, PYRAMID, SPRISM };

// global variables
int dimension = 0;
MESH_POLYHEDRON_TYPE mesh_polyhedron_type = CUBE;
bool verbose = true;
int output_trigger = 5000;
const int MAX_HYPERCUBE_DIMENSION = 4;
const int MAX_PYRAMID_DIMENSION = 5;
const int MAX_SIMPLEX_DIMENSION = 20;
bool genivoltable = false;
bool nep_flag = false;
char * output_filename = NULL;
const char * table_string = NULL;
const char * isosurface_table_string = "isosurface table";
const char * interval_volume_table_string = "interval volume table";

// routines
void generate_isosurface_table(ISOSURFACE_TABLE & isotable);
void generate_nep_isosurface_table(ISOSURFACE_TABLE & isotable);
void generate_interval_volume_table(ISOSURFACE_TABLE & isotable);
void parse_command_line(int argc, char **argv);
void check_dimension(const int d, 
		     const MESH_POLYHEDRON_TYPE mesh_poly_type);
void check_num_vertices(const int d,
			const MESH_POLYHEDRON_TYPE mesh_poly_type,
			const int max_numv);

void memory_exhaustion();
void usage_error();
void help_msg();

// external routine
void ijkgenpatch
  (const ISOSURFACE_TABLE_POLYHEDRON & polyhedron,
   const int * const vertex_sign, 
   vector<ISOSURFACE_VERTEX_INDEX> & edge_list, int & num_simplices);

void ijkgenpatch_nep
  (const ISOSURFACE_TABLE_POLYHEDRON & polyhedron, 
   const int * const vertex_sign,
   vector<ISOSURFACE_VERTEX_INDEX> & isov_list, 
   vector<ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE> & isov_type,
   int & num_simplices);


int main(int argc, char **argv)
{
  std::set_new_handler(memory_exhaustion);
  ostringstream isotable_filename;

  try {

    parse_command_line(argc, argv);

    if (dimension == 0) {
      cout << "Enter dimension: ";
      cin >> dimension;
    };

    check_dimension(dimension, mesh_polyhedron_type);

    ISOSURFACE_TABLE isotable(dimension);

    if (!genivoltable) {
      table_string = isosurface_table_string;
      isotable.SetSimplexDimension(dimension-1);
      if (!nep_flag) { 
	isotable.SetBinaryEncoding(); 
      }
      else { 
	isotable.SetBase3Encoding(); 
      };
    }
    else {
      assert(!nep_flag);

      table_string = interval_volume_table_string;
      isotable.SetBase3Encoding();
      isotable.SetSimplexDimension(dimension);
    };

    check_num_vertices(dimension, mesh_polyhedron_type, 
		       isotable.MaxNumVertices());

    switch(mesh_polyhedron_type) {

    case CUBE:
    default:
      isotable.GenCube(dimension);
      break;

    case SIMPLEX:
      isotable.GenSimplex(dimension);
      break;

    case PYRAMID:
      isotable.GenPyramid(dimension);
      break;

    case SPRISM:
      ISOSURFACE_TABLE_POLYHEDRON base_polyhedron(dimension-1), 
	prism(dimension);
      base_polyhedron.GenSimplex(dimension-1);
      generate_prism(base_polyhedron, prism);
      isotable.Set(prism);
      break;
    };

    ERROR error;
    if (!isotable.Polyhedron().Check(error)) {
      cerr << "Problem constructing polyhedron." << endl;
      error.Print(cerr);
      cerr << "Exiting." << endl;
      exit(10);
    };

    u_long num_table_entries = 0;
    int nume = isotable.Polyhedron().NumEdges();
    int numv = isotable.Polyhedron().NumVertices();
    if (!genivoltable) {
      if (!nep_flag) {
	int num_isosurface_vertices = nume;
	isotable.SetNumIsosurfaceVertices(num_isosurface_vertices);
	isotable.StorePolyEdgesAsIsoVertices(0);
	num_table_entries = calculate_num_entries(numv, 2);
      }
      else {
	int num_isosurface_vertices = nume+numv;
	isotable.SetNumIsosurfaceVertices(num_isosurface_vertices);
	isotable.StorePolyVerticesAsIsoVertices(0);
	isotable.StorePolyEdgesAsIsoVertices(numv);
	num_table_entries = calculate_num_entries(numv, 3);
      }
    }
    else {
      int num_isosurface_vertices = 2*nume+numv;
      isotable.SetNumIsosurfaceVertices(num_isosurface_vertices);
      for (int ie = 0; ie < nume; ie++) {
	isotable.SetIsoVertexFace(2*ie, ie);
	isotable.SetIsoVertexType(2*ie, ISOSURFACE_VERTEX::EDGE);
	isotable.SetIsoVertexLabel(2*ie, "0");
	isotable.SetIsoVertexFace(2*ie+1, ie);
	isotable.SetIsoVertexType(2*ie+1, ISOSURFACE_VERTEX::EDGE);
	isotable.SetIsoVertexLabel(2*ie+1, "1");
      };

      for (int iv = 0; iv < numv; iv++) {
	isotable.SetIsoVertexFace(iv+2*nume, iv);
	isotable.SetIsoVertexType(iv+2*nume, ISOSURFACE_VERTEX::VERTEX);
      }
      num_table_entries = calculate_num_entries(numv, 3);
    };

    isotable.SetNumTableEntries(num_table_entries);

    if (verbose) {
      cout << isotable.Polyhedron().Dimension() << "D ";
      switch(mesh_polyhedron_type) {

      case SIMPLEX:
	cout << "simplex has ";
	break;

      case CUBE:
	cout << "hypercube has ";
	break;

      case PYRAMID:
	cout << "pyramid has ";
	break;

      case SPRISM:
	cout << "prism with simplex base has ";
	break;
      };

      cout << isotable.Polyhedron().NumVertices()
	   << " vertices and " << isotable.Polyhedron().NumEdges()
	   << " edges." << endl;
      string table_string2 = table_string;
      if (table_string2.length() > 0) {
	table_string2[0] = toupper(table_string2[0]);
      }
      cout << table_string2 << " has "
	   << isotable.NumTableEntries() << " entries." << endl;
    };

    if (verbose) {
      cout << "Generating " << table_string << "." << endl;
    }

    if (!genivoltable) {
      if (!nep_flag) {
	generate_isosurface_table(isotable);
      }
      else {
	generate_nep_isosurface_table(isotable);
      }
    }
    else {
      generate_interval_volume_table(isotable);
    };

    if (!isotable.CheckTable(error)) {
      cerr << "Error detected in isosurface table." << endl;
      error.Print(cerr);
      cerr << "Exiting." << endl;
      exit(10);
    };

    if (output_filename == NULL) {

      // generate isosurface table file name
      isotable_filename.str("");
      if (!genivoltable) {
	isotable_filename << "iso.";
      }
      else {
	isotable_filename << "ivol.";
      };
      if (nep_flag) {
	isotable_filename << "nep.";
      }
      switch (mesh_polyhedron_type) {

      case SIMPLEX:
	isotable_filename << "simplex";
	break;

      case PYRAMID:
	isotable_filename << "pyramid";
	break;

      case SPRISM:
	isotable_filename << "sprism";
	break;

      case CUBE:
      default:
	isotable_filename << "cube";
	break;
      };

      isotable_filename << "." << isotable.Dimension() 
			<< "D.xit" << ends; 
    }
    else {
      isotable_filename << output_filename;
    }

    if (verbose)
      cout << "Writing " << table_string << " to file " 
	   << isotable_filename.str() << "." << endl;

    ofstream isotable_file(isotable_filename.str().c_str(), ios::out);
    if (!isotable_file) {
      cerr << "Unable to open file " << isotable_filename.str() 
	   << endl;
      cerr << "Exiting." << endl;
      exit(15);
    };

    try {
      IJKXIO::write_xit(isotable_file, isotable);
    }
    catch(...) {
      cerr << "Error writing file " << isotable_filename.str() 
	   << "." << endl;
      throw;
    }
  }
  catch (ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

  exit(0);
}

void generate_isosurface_table(ISOSURFACE_TABLE & isotable)
{
  int vertex_sign[isotable.Polyhedron().NumVertices()];
  vector<ISOSURFACE_VERTEX_INDEX> edge_list;
  const int numv = isotable.Polyhedron().NumVertices();

  for (int it = 0; it < isotable.NumTableEntries(); it++) {
    edge_list.clear();
    int num_simplices = 0;

    convert2base(it, 2, vertex_sign, numv);
    for (int j = 0; j < numv; j++) {
      // convert 0 to -1
      vertex_sign[j] = 2*vertex_sign[j]-1;
      // if vertex_sign[j] < 0, vertex j is negative
      // if vertex_sign[j] > 0, vertex j is positive
    }

    ijkgenpatch(isotable.Polyhedron(), vertex_sign, edge_list, num_simplices);

    isotable.SetNumSimplices(it, num_simplices);

    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.Dimension(); j++) {
	int ie = edge_list[is*isotable.Dimension() + j];
	isotable.SetSimplexVertex(it, is, j, ie);
      };
    };

    if (verbose && it > 0 && it%output_trigger == 0) {
      cout << "  " << it << " out of " << isotable.NumTableEntries()
	   << " isosurface table entries completed." << endl; 
    };
  };

  if (verbose && isotable.NumTableEntries() && 
      isotable.NumTableEntries() > output_trigger) {
      cout << "  All " << isotable.NumTableEntries()
	   << " isosurface table entries completed." << endl; 
  }

}

void generate_nep_isosurface_table(ISOSURFACE_TABLE & isotable)
{
  int vertex_sign[isotable.Polyhedron().NumVertices()];
  vector<ISOSURFACE_VERTEX_INDEX> isov_list;
  vector<ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE> isov_type;
  const int numv = isotable.Polyhedron().NumVertices();

  for (int it = 0; it < isotable.NumTableEntries(); it++) {
    isov_list.clear();
    int num_simplices = 0;

    convert2base(it, 3, vertex_sign, numv);
    for (int j = 0; j < numv; j++) {
      // convert 0 to -1, 1 to 0 and 2 to 1
      vertex_sign[j] = vertex_sign[j]-1;
      // if vertex_sign[j] < 0, vertex j is negative
      // if vertex_sign[j] == 0, vertex j is zero
      // if vertex_sign[j] > 0, vertex j is positive
    }

    ijkgenpatch_nep(isotable.Polyhedron(), vertex_sign, 
		    isov_list, isov_type, num_simplices);
    assert(isov_list.size() == isov_type.size());

    isotable.SetNumSimplices(it, num_simplices);

    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.Dimension(); j++) {
	int k = is*isotable.Dimension() + j;
	int kf = isov_list[k];
	if (isov_type[k] == ISOSURFACE_VERTEX::VERTEX) {
	  isotable.SetSimplexVertex(it, is, j, kf);
	}
	else {
	  assert(isov_type[k] == ISOSURFACE_VERTEX::EDGE);
	  isotable.SetSimplexVertex(it, is, j, kf+numv);
	}
      };
    };

    if (verbose && it > 0 && it%output_trigger == 0) {
      cout << "  " << it << " out of " << isotable.NumTableEntries()
	   << " isosurface table entries completed." << endl; 
    };
  };

  if (verbose && isotable.NumTableEntries() && 
      isotable.NumTableEntries() > output_trigger) {
      cout << "  All " << isotable.NumTableEntries()
	   << " isosurface table entries completed." << endl; 
  }

}

void generate_interval_volume_table(ISOSURFACE_TABLE & isotable)
{
  int index_base3[isotable.Polyhedron().NumVertices()];
  int vertex_prism_sign[2*isotable.Polyhedron().NumVertices()];
  vector<ISOSURFACE_VERTEX_INDEX> edge_list;
  const int numv = isotable.Polyhedron().NumVertices();
  const int nume = isotable.Polyhedron().NumEdges();
  PROCEDURE_ERROR error("generate_interval_volume_table");

  ISOSURFACE_TABLE_POLYHEDRON prism(isotable.Dimension()+1);

  generate_prism(isotable.Polyhedron(), prism);
  const int num_prism_vertices = prism.NumVertices();

  for (int it = 0; it < isotable.NumTableEntries(); it++) {
    edge_list.clear();
    int num_simplices = 0;

    convert2base(it, 3, index_base3, numv);
    for (int i = 0; i < numv; i++) {
      if (index_base3[i] == 0) {
	vertex_prism_sign[i] = -1;
	vertex_prism_sign[i+numv] = -1;
      }
      else if (index_base3[i] == 2) {
	vertex_prism_sign[i] = 1;
	vertex_prism_sign[i+numv] = 1;
      }
      else {
	vertex_prism_sign[i] = 1;
	vertex_prism_sign[i+numv] = -1;
      }
    }

    ijkgenpatch(prism, vertex_prism_sign, edge_list, num_simplices);

    isotable.SetNumSimplices(it, num_simplices);

    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
	int prism_ie = edge_list[is*isotable.NumVerticesPerSimplex() + j];
	int ie;

	if (prism_ie < nume) {
	  ie = 2*prism_ie;
	}
	else if (prism_ie < 2*nume) {
	  ie = 2*(prism_ie-nume)+1;
	}
	else if (prism_ie < 2*nume+numv) {
	  ie = prism_ie;
	}
	else {
	  error.AddMessage("Programming error. Illegal edge ", ie, 
			   " in isosurface vertex list");
	  throw error;
	}
	isotable.SetSimplexVertex(it, is, j, ie);

      };
    };

    if (verbose && it > 0 && it%output_trigger == 0 &&
	isotable.NumTableEntries() > output_trigger) {
      cout << "  " << it << " out of " << isotable.NumTableEntries()
	   << " interval volume table entries completed." << endl; 
    };
  };

  if (verbose && isotable.NumTableEntries() && 
      isotable.NumTableEntries() > output_trigger) {
      cout << "  All " << isotable.NumTableEntries()
	   << " interval volume table entries completed." << endl; 
  }
}

void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc) {
    if (argv[iarg][0] != '-')
      usage_error();
    if (strcmp(argv[iarg], "-d") == 0) {
      iarg++;
      if (iarg >= argc)
	usage_error();

      sscanf(argv[iarg], "%d", &dimension);
    }
    else if (strcmp(argv[iarg], "-V") == 0) {
      verbose = true;
      iarg++;
      if (iarg >= argc)
	usage_error();

      sscanf(argv[iarg], "%d", &output_trigger);
    }
    else if (strcmp(argv[iarg], "-poly") == 0) {
      iarg++;
      if (iarg >= argc)
	usage_error();

      if (strcmp(argv[iarg], "cube") == 0)
	  mesh_polyhedron_type = CUBE;
      else if (strcmp(argv[iarg], "simplex") == 0)
	  mesh_polyhedron_type = SIMPLEX;
      else if (strcmp(argv[iarg], "pyramid") == 0)
	  mesh_polyhedron_type = PYRAMID;
      else if (strcmp(argv[iarg], "sprism") == 0)
	  mesh_polyhedron_type = SPRISM;
      else
	usage_error();
    }
    else if (strcmp(argv[iarg], "-ivol") == 0) {
      genivoltable = true;
    }
    else if (strcmp(argv[iarg], "-nep") == 0) {
      nep_flag = true;
    }
    else if (strcmp(argv[iarg], "-o") == 0) {
      iarg++;
      if (iarg >= argc)
	usage_error();

      output_filename = argv[iarg];
    }
    else {
      for (int j = 1; argv[iarg][j] != '\0'; j++) {
	if (argv[iarg][j] == 'c')
	  mesh_polyhedron_type = CUBE;
	else if (argv[iarg][j] == 'h')
	  help_msg();
	else if (argv[iarg][j] == 'q')
	  verbose = false;
	else if (argv[iarg][j] == 's')
	  mesh_polyhedron_type = SIMPLEX;
	else if (argv[iarg][j] == 'v')
	  verbose = true;
	else
	  usage_error();
      };
    };
    iarg++;
  };


  if (genivoltable && nep_flag) {
    cerr << "Error. Interval volume not yet implemented with -nep flag."
	 << endl;
    exit(10);
  }

}

void check_dimension(const int dimension,
		     const MESH_POLYHEDRON_TYPE mesh_polyhedron_type)
// exit if illegal dimension found
{
  if (dimension < 2) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Dimension must be at least 2." << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };

  if (mesh_polyhedron_type == CUBE && dimension > MAX_HYPERCUBE_DIMENSION) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Hypercube dimension is at most " 
	 << MAX_HYPERCUBE_DIMENSION << "." 
	 << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };

  if (mesh_polyhedron_type == SIMPLEX && dimension > MAX_SIMPLEX_DIMENSION) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Simplex dimension is at most " << MAX_SIMPLEX_DIMENSION << "." 
	 << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };

  if (mesh_polyhedron_type == PYRAMID && dimension > MAX_PYRAMID_DIMENSION) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Pyramid dimension is at most " << MAX_PYRAMID_DIMENSION << "." 
	 << endl;
    cerr << "Exiting." << endl;
    exit(10);
  }
}

void check_num_vertices(const int d,
			const MESH_POLYHEDRON_TYPE mesh_poly_type,
			const int max_numv)
{
  int numv = 0;
  if (mesh_polyhedron_type == CUBE) {
    numv = (1L << dimension);
  }
  else {
    numv = dimension+1;
  };

  if (numv > max_numv) {
    cerr << "Polyhedron has too many vertices." << endl;
    cerr << d << "D ";
    if (mesh_polyhedron_type == CUBE) {
      cerr << "hypercube ";
    }
    else {
      cerr << "simplex ";
    };
    cerr << "has " << numv << " vertices.  Maximum number of vertices = " 
	 << max_numv << "."
	 << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };
}

void memory_exhaustion()
{
   cerr << "Error: Out of memory.  Terminating program." << endl;
   exit(10);
}

void usage_msg()
{
  cout << "Usage: ijkgentable [-chqsv] [-ivol] [-nep] [-d <dimension>] [-V <n>]" 
       << endl;
  cout << "                   [-poly {cube|simplex|pyramid|sprism}]" << endl;
  cout << "                   [-o {output_filename}]" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help_msg()
{
  cout << "ijkgentable - generate isosurface or interval volume table" << endl;
  cout << "  Generate table of isosurface patches for all +/- vertex sign patterns." << endl;
  usage_msg();
  cout << "  -ivol: Generate interval volume lookup table." << endl;
  cout << "  -nep: Negative, equals, positive." << endl;
  cout << "        Differentiate scalar values equal to the isovalue" << endl;
  cout << "        from scalar values less (-) or greater (+) than the isovalue." << endl;
  cout << "  -poly cube: Cube polyhedron (default.)" << endl;
  cout << "  -poly simplex: Simplex polyhedron." << endl;
  cout << "  -poly pyramid: Pyramid with cube (square) base." << endl;
  cout << "  -poly sprism: Prism with simplex base." << endl;
  cout << "  -d <dimension> : Dimension is <dimension>." << endl;
  cout << "  -h : Print this help message (and exit.)" << endl;
  cout << "  -o {fname} : Output table to file {fname}." << endl;
  cout << "               Without this option, output file name is built " << endl;
  cout << "               from poly and dimension." << endl;
  cout << "  -q : Quiet mode." << endl;
  cout << "  -v : Verbose mode (default.)" << endl;
  cout << "  -V <n>: Verbose mode. Output after every <n> table entries." 
       << endl;
  exit(0);
}
