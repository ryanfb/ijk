/// \file ijkmcube.cxx
/// generate isosurface from scalar field
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

#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ijkmcube.h"
#include "ijkxitIO.h"
#include "ijkNrrd.h"
#include "ijkIO.txx"

using namespace IJK;
using namespace IJKGRID;
using namespace IJKMCUBE;
using namespace IJKTABLE;/// \file ijkmcube.h

using namespace std;

typedef enum OUTPUT_FORMAT { OFF, IV };
struct IO_TIME;

// global variables
vector<SCALAR_TYPE> isovalue;
vector<string> isovalue_string;
char * input_filename = NULL;
char * output_filename = NULL;
const char * isotable_directory = NULL;
bool genivol = false;
OUTPUT_FORMAT output_format = OFF;
bool report_time_flag = false;
bool use_stdout = false;
bool use_minmax_regions = false;
bool use_octree = false;
bool use_nep_isotable = false;
bool use_list = false;
int nep_num_dup = 2;
bool use_snapmc = false;
SNAP_TYPE snap_value = 0.0;
int region_length;
bool nowrite_flag = false;
bool flag_silent = false;
bool flag_subsample = false;
int subsample_resolution = 2;

// routines
void parse_command_line(int argc, char **argv);
void read_isosurface_table
(const int dimension, ISOSURFACE_TABLE & isotable, 
 string & isotable_filename, IO_TIME & io_time);
void check_octree(const IJKOCTREE::OCTREE & octree);
void report_time
(const IO_TIME & io_time, const MCUBE_TIME & mcube_time, 
 const double total_elapsed_time,
 const char * input_filename, const char * isotable_filename);
void report_iso_info
(const ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> &slist, 
 const GRID_INFO & mcube_grid);
void memory_exhaustion();
void help(), usage_error();
void split_string(const string & s, const char c,
		  string & prefix, string & suffix);
void rescale_coord(const int scale, vector<COORD_TYPE> & vertex_coord);
void set_minmax
(const MC_GRID & scalar_grid, const AXIS_SIZE_TYPE region_length,
 MINMAX_DATA & minmax);
void write_mesh
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist);
void write_mesh
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 IO_TIME & io_time);
void output_isosurface
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const GRID_INFO & grid_info, IO_TIME & io_time);
void output_isosurface
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const NEP_INFO & nep_info, IO_TIME & io_time);


// *** DEBUG ***
void write_snap_nrrd
(const MC_GRID & scalar_grid, const SCALAR_TYPE isovalue,  
 const SNAP_TYPE snap_value);


typedef enum PARAMETER 
  {REGION_PARAM, OCTREE_PARAM, LIST_PARAM, NEP_PARAM, NEP_DUP_PARAM,
   IVOL_PARAM, SNAP_PARAM,
   SUBSAMPLE_PARAM, HELP_PARAM, OFF_PARAM, IV_PARAM, 
   DIR_PARAM, OUTPUT_FILENAME_PARAM, STDOUT_PARAM, 
   NOWRITE_PARAM, SILENT_PARAM, TIME_PARAM, UNKNOWN_PARAM};
char * parameter_string[] = 
  {"-region", "-octree", "-list", 
   "-nep", "-nep_dup", "-ivol", "-snap", "-subsample",
   "-help", "-off", "-iv", "-dir", "-o", "-stdout",
   "-nowrite", "-s", "-time", "-unknown"};

// timing functions
class ELAPSED_CPU_TIME {

protected:
  clock_t t;

public:
  ELAPSED_CPU_TIME() { t = clock(); };

  clock_t getElapsed() {
    clock_t old_t = t;
    t = clock();
    return(t - old_t);
  };
};

class ELAPSED_TIME {

protected:
  time_t t;

public:
  ELAPSED_TIME() { time(&t);  };

  double getElapsed() {
    time_t old_t = t;
    time(&t);
    return(difftime(t,old_t));
  };
};

struct IO_TIME {
  double read_table_time;   // wall time to read isosurface lookup table
  double read_nrrd_time;    // wall time to read nrrd file
  double write_time;        // wall time to write output
};


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);
  string isotable_filename;
  MERGE_DATA * merge_data = NULL;

  ELAPSED_TIME wall_time;

  MCUBE_TIME mcube_time;
  IO_TIME io_time = {0.0, 0.0, 0.0};

  AXIS_SIZE_TYPE * grid_length = NULL;
  SCALAR_TYPE * scalar = NULL;

  try {

    std::set_new_handler(memory_exhaustion);

    isotable_directory = getenv("IJK_ISOTABLE_DIR");

#ifdef IJK_ISOTABLE_DIR
    if (isotable_directory == NULL) {
      isotable_directory = IJK_ISOTABLE_DIR;
    }
#endif

    parse_command_line(argc, argv);

    if (isotable_directory == NULL) {
      cerr << "Error.  Unknown isotable directory." << endl;
      cerr << "  Use -dir {isotable_directory} argument or set environment variable IJK_ISOTABLE_DIR." << endl;
      exit(5);
    };

    if (!genivol) {
      if (isovalue.size() > 1 && use_stdout) {
	cerr << "Error.  Cannot use stdout for more than one isovalue." << endl;
	cerr << "Exiting.";
	exit(10);
      };
      if (isovalue.size() > 1 && output_filename != NULL) {
	cerr << "Error.  Cannot specify output file for more than one isovalue." << endl;
	cerr << "Exiting.";
	exit(10);
      };
    }
    else {
      if (isovalue.size() > 2 && use_stdout) {
	cerr << "Error.  Cannot use stdout for more than one interval volume." 
	     << endl;
	cerr << "Exiting.";
	exit(10);
      };
      if (isovalue.size() > 2 && output_filename != NULL) {
	cerr << "Error.  Cannot specify output file for more than one interval volume." 
	     << endl;
	cerr << "Exiting.";
	exit(10);
      };
    };

    int dimension = 0;
    AXIS_SIZE_TYPE * full_grid_length = NULL;
    SCALAR_TYPE * full_scalar = NULL;
    Nrrd *nin;

    /* get scalar field data from nrrd file */
    nin = nrrdNew();
    if (nrrdLoad(nin, input_filename, NULL)) {
      char *err = biffGetDone(NRRD);
      cerr << "Error reading: " << input_filename << endl;
      cerr << "  Error: " << err << endl;
      exit(35);
    };
    dimension = nin->dim;

    if (dimension < 1) {
      cerr << "Illegal dimension.  Dimension must be at least 1." << endl;
      exit(20);
    };

    full_grid_length = new AXIS_SIZE_TYPE[dimension];
    size_t size[NRRD_DIM_MAX];

    nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size); 

    for (int d = 0; d < dimension; d++)
      full_grid_length[d] = size[d];

    VERTEX_INDEX num_grid_vertices = 
      compute_num_grid_vertices(dimension, full_grid_length);

    full_scalar = new SCALAR_TYPE[num_grid_vertices];

    nrrd2scalar(nin, full_scalar);


    // *** TEST ComputeRegionMinMax() ***
    IJKGRID::MINMAX_GRID<int, AXIS_SIZE_TYPE, VERTEX_INDEX, SCALAR_TYPE>
      minmax_grid;
    const int offset_length = 1;
    minmax_grid.ComputeRegionMinMax
      (dimension, full_grid_length, full_scalar, region_length, offset_length);

    if (dimension == 1) {
      for (int i0 = 0; i0 < minmax_grid.AxisSize(0); i0++) {
	cout << " (" << minmax_grid.Min(i0)
	     << "," << minmax_grid.Max(i0) << ")";
      }
      cout << endl;
    }
    else if (dimension == 2) {
      for (int i1 = 0; i1 < minmax_grid.AxisSize(1); i1++) {
	for (int i0 = 0; i0 < minmax_grid.AxisSize(0); i0++) {
	  VERTEX_INDEX iv = i1*minmax_grid.AxisSize(0) + i0;
	  cout << " (" << minmax_grid.Min(iv)
	       << "," << minmax_grid.Max(iv) << ")";
	}
	cout << endl;
      }
    }
    else if (dimension == 3) {
      for (int i2 = 0; i2 < minmax_grid.AxisSize(2); i2++) {
	for (int i1 = 0; i1 < minmax_grid.AxisSize(1); i1++) {
	  for (int i0 = 0; i0 < minmax_grid.AxisSize(0); i0++) {
	    VERTEX_INDEX iv = 
	      (i2*minmax_grid.AxisSize(1)+i1)*minmax_grid.AxisSize(0) + i0;
	    cout << " (" << minmax_grid.Min(iv)
		 << "," << minmax_grid.Max(iv) << ")";
	  }
	  cout << endl;
	}
	cout << endl << endl;
      }
    }

    cout << endl;
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

}


// **************************************************
// OUTPUT
// **************************************************

void output_isosurface
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const GRID_INFO & grid_info, IO_TIME & io_time)
{
  if (!nowrite_flag) 
    { write_mesh(i, isotable, vertex_coord, slist, io_time); }

  if (!use_stdout && !flag_silent) {
    report_iso_info(isotable, isovalue[i], vertex_coord, slist, 
		    grid_info);
  }
}


void output_isosurface
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const NEP_INFO & nep_info, IO_TIME & io_time)
{
  output_isosurface(i, isotable, vertex_coord, slist, nep_info.grid, io_time);

  if (!use_stdout && !flag_silent) {
    cerr << "    " << nep_info.num_non_empty_boundary_facets
	 << " non-empty boundary cube facets." << endl;

    if (nep_num_dup != 2) {
      cerr << "    " << nep_info.num_in_facet_cubes
	   << " cubes with isosurface patches contained in a facet." << endl;
      cerr << "    " << nep_info.num_dup_iso_patches
	   << " duplicate isosurface patches eliminated." << endl;
    }
  }
}


// **************************************************
// WRITE_MESH
// **************************************************

void write_mesh
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist)
{
  ERROR error_mcube("ijkmcube");
  const int numv_per_simplex = isotable.NumVerticesPerSimplex();
  const int dimension = isotable.Dimension();
  ofstream output_file;

  // create output filename
  string fname = input_filename;

  // remove path from file name
  string prefix, suffix;
  split_string(fname, '/', prefix, suffix);
  if (suffix != "") { fname = suffix; }

  string ofilename;

  if (output_filename != NULL) {
    ofilename = output_filename;
  }
  else {
    // construct output filename
    split_string(fname, '.', prefix, suffix);
    if (suffix == "nrrd" || suffix == "nhdr") {
      ofilename = prefix;
    }
    else {
      ofilename = input_filename;
    }
    if (!genivol) {
      ofilename += string(".") + string("isov=") + isovalue_string[i];
    }
    else {
      ofilename += string(".") + string("ivol=") + isovalue_string[i]
	+ "," + isovalue_string[i+1];
    }

    switch (output_format) {
    case OFF: 
      ofilename += ".off";
      break;

    case IV:
      ofilename += ".iv";
      break;
    };
  }

  switch (output_format) {

  case OFF:
    if (!use_stdout) {
      output_file.open(ofilename.c_str(), ios::out);
      IJKIO::ijkoutOFF(output_file, dimension, numv_per_simplex,
		       vertex_coord, slist);
      output_file.close();
    }
    else {
      IJKIO::ijkoutOFF(dimension, numv_per_simplex, vertex_coord, slist);
    };
    break;

  case IV:
    if (dimension == 3) {
      if (!use_stdout) {
	output_file.open(ofilename.c_str(), ios::out);
	IJKIO::ijkoutIV(output_file, dimension, vertex_coord, slist);
	output_file.close();
      }
      else {
	IJKIO::ijkoutOFF(dimension, vertex_coord, slist);
      }
    }
    else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
    break;

  default:
    throw error_mcube("Illegal output format.");
    break;
  }

  if (!use_stdout && !flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}

void write_mesh
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_mesh(i, isotable, vertex_coord, slist);

  if (report_time_flag) {
    io_time.write_time += wall_time.getElapsed();
  };
}

// **************************************************
// MISCELLANEOUS ROUTINES
// **************************************************

void read_isosurface_table
(const int dimension, ISOSURFACE_TABLE & isotable, 
 string & isotable_filename, IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  if (report_time_flag) { 
    io_time.read_table_time = wall_time.getElapsed(); 
  };

  string isotable_pathname;

  if (use_nep_isotable || use_snapmc) {
    isotable_filename = get_isotable_filename(NEP, dimension);
  }
  else if (genivol) {
    isotable_filename = get_isotable_filename(IVOL, dimension);
  }
  else {
    isotable_filename = get_isotable_filename(BINARY, dimension);
  }

  isotable_pathname = string(isotable_directory) + "/" + isotable_filename;

  ifstream isotable_file(isotable_pathname.c_str(), ios::in);
  if (!isotable_file) {
    isotable_file.clear();
    isotable_file.open(isotable_filename.c_str(), ios::in);
  };

  if (!isotable_file) {
    cerr << "Unable to obtain isosurface table file "
	 << isotable_pathname << "." << endl;
    exit(30);
  };

  try {
    IJKXIO::read_xit(isotable_file, isotable);
  }
  catch(...) {
    cerr << "Error reading file: " << isotable_filename << "." << endl;
    throw;
  };

  ERROR error_msg;
  if (!isotable.Check(error_msg)) {
    cerr << "Warning: Data structure inconsistency in isosurface table "
	 << isotable_pathname << "." << endl;
    for (int i = 0; i < error_msg.NumMessages(); i++) 
      cerr << error_msg.Message(i) << endl;
    cerr << "  Attempting to continue..." << endl << endl;
  }

  if (isotable.Dimension() != dimension) {
    cerr << "Error.  Dimension " << isotable.Dimension() 
	 << " of isosurface table in file " << isotable_filename << endl;
    cerr << "  does not match dimension " << dimension
	 << " of nrrd data in file " << input_filename << "." << endl;
    cerr << "  Exiting." << endl;
    exit(33);
  }

  if (report_time_flag) { 
    io_time.read_table_time = wall_time.getElapsed(); 
  };

}

void rescale_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
{
  for (int i = 0; i < vertex_coord.size(); i++) {
    vertex_coord[i] = scale * vertex_coord[i];
  };
}


void set_minmax
(const MC_GRID & scalar_grid, const AXIS_SIZE_TYPE region_length,
 MINMAX_DATA & minmax)
// set minmax data structure
{
  minmax.Create(scalar_grid, region_length);
}

void report_iso_info
(const ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> &slist, 
 const GRID_INFO & grid_info)
{
  VERTEX_INDEX numv = (vertex_coord.size())/isotable.Dimension();
  VERTEX_INDEX nums = (slist.size())/isotable.NumVerticesPerSimplex();
  VERTEX_INDEX num_grid_cubes = grid_info.num_cubes;
  VERTEX_INDEX num_mixed_cubes = grid_info.num_mixed_cubes;
  float percent = float(num_mixed_cubes)/float(num_grid_cubes);
  int ipercent = int(100*percent);
  cerr << "  Isovalue " << isovalue << ".  " 
       << numv << " isosurface vertices.  "
       << nums << " isosurface simplices." << endl;
  cerr << "    " << num_mixed_cubes
       << " (" << ipercent << "%) non-empty cubes." << endl;

}

void report_time
(const IO_TIME & io_time, const MCUBE_TIME & mcube_time, 
 const double total_elapsed_time,
 const char * input_filename, const char * isotable_filename)
{
  const char * ISOSURFACE_STRING = "isosurface";
  const char * INTERVAL_VOLUME_STRING = "interval volume";
  const char * mesh_type_string = NULL;
  
  if (!genivol) {
    mesh_type_string = ISOSURFACE_STRING;
  }
  else {
    mesh_type_string = INTERVAL_VOLUME_STRING;
  };

  cerr << "Time to read file " << input_filename << ": "
       << io_time.read_nrrd_time << " seconds." << endl;
  cerr << "Time to read " << mesh_type_string << " lookup table "
       << isotable_filename << ": " 
       << io_time.read_table_time << " seconds." << endl;
  cerr << "CPU time to run Marching Cubes: " 
       << mcube_time.total << " seconds." << endl;
  if (use_octree) {
    cerr << "    Time to create octree: "
	 << mcube_time.preprocessing << " seconds." << endl;
  }
  if (use_minmax_regions) {
    cerr << "    Time to create min/max regions: "
	 << mcube_time.preprocessing << " seconds." << endl;
  }
  if (use_snapmc) {
    cerr << "    Time to snap grid scalar values: "
	 << mcube_time.snap << " seconds." << endl;
  }
  if (use_list) {
    cerr << "    Time to create list of non-empty cubes: "
	 << mcube_time.preprocessing << " seconds." << endl;
  }
  cerr << "    Time to extract " << mesh_type_string << " triangles: "
       << mcube_time.extract << " seconds." << endl;
  cerr << "    Time to merge identical "
       << mesh_type_string << " vertices: " 
       << mcube_time.merge << " seconds." << endl;
  cerr << "    Time to position "
       << mesh_type_string << " vertices: "
       << mcube_time.position << " seconds." << endl;
  if (!nowrite_flag) {
    cerr << "Time to write "
	 << mesh_type_string << ": " 
	 << io_time.write_time << " seconds." << endl;
  };
  cerr << "Total elapsed time: " << total_elapsed_time
       << " seconds." << endl;

}

void memory_exhaustion()
{
   cerr << "Error: Out of memory.  Terminating program." << endl;
   exit(10);
}

PARAMETER get_parameter_token(char * s)
// convert string s into parameter token
{
  for (int i = 0; i < int(UNKNOWN_PARAM); i++)
    if (strcmp(parameter_string[i], s) == 0)
      return(PARAMETER(i));
  return(UNKNOWN_PARAM);
}

void parse_command_line(int argc, char **argv)
// parse command line
// control parameters, followed by one or more isovalues, 
// followed by input file name
{
  if (argc != 3) { usage_error(); };

  // remaining parameters should be edge list followed by input file name 
  sscanf(argv[1], "%d", &region_length);
  input_filename = argv[2];
}

void usage_msg()
{
  cerr << "Usage: test_project {region_edge_length} {input filename}" 
       << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help()
{
  usage_msg();
  exit(20);
}

void split_string(const string & s, const char c,
		  string & prefix, string & suffix)
// split string at last occurrence of character c into prefix and suffix
{
  string::size_type i = s.rfind(c);
  if (i == string::npos) {
    prefix = s;
    suffix = "";
  }
  else {
    if (i > 0) { prefix = s.substr(0,i); }
    else { prefix = ""; };

    if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
    else { suffix = ""; };
  }
}

