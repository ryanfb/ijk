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

// isosurface extraction routines
void apply_mcube
(const MC_GRID & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);
void apply_nep
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);
void apply_snapMC
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);


void extract_isosurfaces_snap
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, SCALAR_TYPE isovalue,
 vector<COORD_TYPE> & vertex_coord, vector<VERTEX_INDEX> & slist,
 NEP_INFO & nep_info);

// interval volume extraction routines
void extract_interval_volumes
(const MC_GRID & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MCUBE_TIME & mcube_time, IO_TIME & io_time);

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

    // *** TEST ***
    MC_GRID full_scalar_grid(dimension, full_grid_length, full_scalar);
    SCALAR_TYPE * scalar_min = new SCALAR_TYPE[num_grid_vertices];
    SCALAR_TYPE * scalar_max = new SCALAR_TYPE[num_grid_vertices];

    std::copy (full_scalar, full_scalar+num_grid_vertices, scalar_min);
    std::copy (full_scalar, full_scalar+num_grid_vertices, scalar_max);

    MINMAX_GRID<int,AXIS_SIZE_TYPE,VERTEX_INDEX,SCALAR_TYPE>
      min_grid(dimension, full_grid_length, scalar_min);
    MINMAX_GRID<int,AXIS_SIZE_TYPE,VERTEX_INDEX,SCALAR_TYPE>
      max_grid(dimension, full_grid_length, scalar_max);

    min_grid.project_min(2, 2);

    if (dimension == 2) {
      for (int i1 = 0; i1 < full_grid_length[1]; i1++) {
	for (int i0 = 0; i0 < full_grid_length[0]; i0++) {
	  VERTEX_INDEX iv = i1*full_grid_length[0] + i0;
	  cout << " " << setw(3) << scalar_min[iv];
	}
	cout << endl;
      }
    }
    else if (dimension == 3) {
      for (int i2 = 0; i2 < full_grid_length[2]; i2++) {
	for (int i1 = 0; i1 < full_grid_length[1]; i1++) {
	  for (int i0 = 0; i0 < full_grid_length[0]; i0++) {
	    VERTEX_INDEX iv = (i2*full_grid_length[1]+i1)*full_grid_length[0] + i0;
	    cout << " " << setw(3) << scalar_min[iv];
	  }
	  cout << endl;
	}
	cout << endl << endl;
      }
    }


    // *** TEST compute_min and compute_max ***
    const GRID_SIZE_TYPE num_regions = 
      compute_num_regions(dimension, full_grid_length, region_length);
    SCALAR_TYPE * region_min = new SCALAR_TYPE[num_regions];
    SCALAR_TYPE * region_max = new SCALAR_TYPE[num_regions];

    compute_region_min(full_scalar_grid, region_length, region_min);
    compute_region_max(full_scalar_grid, region_length, region_max);

    /* OBSOLETE
    AXIS_SIZE_TYPE num_regions_along_axis[dimension];
    compute_num_regions_along_axes(dimension, full_grid_length, region_length,
				   num_regions_along_axis);

    cout << "**** COMPUTE REGION MIN ***" << endl << endl;

    if (dimension == 2) {
      for (int i1 = 0; i1 < num_regions_along_axis[1]; i1++) {
	for (int i0 = 0; i0 < num_regions_along_axis[0]; i0++) {
	  VERTEX_INDEX iv = i1*num_regions_along_axis[0] + i0;
	  cout << " " << setw(3) << region_min[iv];
	}
	cout << endl;
      }
    }
    else if (dimension == 3) {
      for (int i2 = 0; i2 < num_regions_along_axis[2]; i2++) {
	for (int i1 = 0; i1 < num_regions_along_axis[1]; i1++) {
	  for (int i0 = 0; i0 < num_regions_along_axis[0]; i0++) {
	    VERTEX_INDEX iv = (i2*num_regions_along_axis[1]+i1)*num_regions_along_axis[0] + i0;
	    cout << " " << setw(3) << region_min[iv];
	  }
	  cout << endl;
	}
	cout << endl << endl;
      }
    }

    cout << endl;
    cout << "**** COMPUTE REGION MAX ***" << endl << endl;

    if (dimension == 2) {
      for (int i1 = 0; i1 < num_regions_along_axis[1]; i1++) {
	for (int i0 = 0; i0 < num_regions_along_axis[0]; i0++) {
	  VERTEX_INDEX iv = i1*num_regions_along_axis[0] + i0;
	  cout << " " << setw(3) << region_max[iv];
	}
	cout << endl;
      }
    }
    else if (dimension == 3) {
      for (int i2 = 0; i2 < num_regions_along_axis[2]; i2++) {
	for (int i1 = 0; i1 < num_regions_along_axis[1]; i1++) {
	  for (int i0 = 0; i0 < num_regions_along_axis[0]; i0++) {
	    VERTEX_INDEX iv = (i2*num_regions_along_axis[1]+i1)*num_regions_along_axis[0] + i0;
	    cout << " " << setw(3) << region_max[iv];
	  }
	  cout << endl;
	}
	cout << endl << endl;
      }
    }
    */

    return 0;
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
// MARCHING CUBES
// **************************************************

template <class DATASTRUCT> void mcube_using_datastruct
(const MC_GRID & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time,
 const DATASTRUCT & datastruct)
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  MCUBE_INFO mcube_info;
  mcube_info.grid.num_cubes = scalar_grid.ComputeNumCubes();

  io_time.write_time = 0;
  for (int i = 0; i < isovalue.size(); i++) {

    marching_cubes(scalar_grid, isotable, isovalue[i], slist, vertex_coord, 
		   merge_data, datastruct, mcube_info);
    mcube_time.Add(mcube_info.time);

    if (flag_subsample) 
      { rescale_coord(subsample_resolution, vertex_coord); };
    output_isosurface(i, isotable, vertex_coord, slist, 
		      mcube_info.grid, io_time);
  }
}

void mcube_no_datastruct
(const MC_GRID & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  MCUBE_INFO mcube_info;
  mcube_info.grid.num_cubes = scalar_grid.ComputeNumCubes();

  assert (!use_snapmc);

  io_time.write_time = 0;
  for (int i = 0; i < isovalue.size(); i++) {
    marching_cubes(scalar_grid, isotable, isovalue[i], slist, vertex_coord, 
		   merge_data, mcube_info);
    mcube_time.Add(mcube_info.time);

    if (flag_subsample) 
      { rescale_coord(subsample_resolution, vertex_coord); };
    output_isosurface(i, isotable, vertex_coord, slist, 
		      mcube_info.grid, io_time);
  }

}

void apply_mcube
(const MC_GRID & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  ELAPSED_CPU_TIME cpu_time;

  cpu_time.getElapsed();  // initialize cpu time

  if (use_octree) {
    // create octree
    IJKOCTREE::OCTREE octree(dimension, axis_size);
    octree.SetMinMax(scalar_grid.Scalar());
    check_octree(octree);

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;

    mcube_using_datastruct(scalar_grid, isotable, merge_data, 
			   mcube_time, io_time, octree);
  } 
  else if (use_minmax_regions) {
    MINMAX_DATA minmax;
    set_minmax(scalar_grid, region_length, minmax);

    if (use_snapmc) {
      cerr << "Error.  Minmax regions not implemented with snapMC." << endl;
      exit(31);
    }

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;

    mcube_using_datastruct(scalar_grid, isotable, merge_data, 
			   mcube_time, io_time, minmax);
  }
  else {
    // Marching Cubes algorithm with no preprocessing
    mcube_no_datastruct
      (scalar_grid, isotable, merge_data, mcube_time, io_time);
  }

}


// **************************************************
// NEP (NEGATIVE-EQUALS-POSITIVE) MARCHING CUBES 
// **************************************************

template <class DATASTRUCT> void nep_using_datastruct
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time,
 const DATASTRUCT & datastruct)
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  NEP_INFO nep_info;
  nep_info.grid.num_cubes = scalar_grid.ComputeNumCubes();

  io_time.write_time = 0;
  for (int i = 0; i < isovalue.size(); i++) {

    marching_cubes_nep
      (scalar_grid, isotable, isovalue[i], nep_num_dup,
       slist, vertex_coord, merge_data, datastruct, nep_info);
    mcube_time.Add(nep_info.time);

    if (flag_subsample) 
      { rescale_coord(subsample_resolution, vertex_coord); };
    output_isosurface(i, isotable, vertex_coord, slist, nep_info, io_time);
  }
}

void nep_no_datastruct
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  NEP_INFO nep_info;
  nep_info.grid.num_cubes = scalar_grid.ComputeNumCubes();

  io_time.write_time = 0;
  for (int i = 0; i < isovalue.size(); i++) {

    marching_cubes_nep
      (scalar_grid, isotable, isovalue[i], nep_num_dup,
       slist, vertex_coord, merge_data, nep_info);
    mcube_time.Add(nep_info.time);

    if (flag_subsample) { rescale_coord(subsample_resolution, vertex_coord); };
    output_isosurface(i, isotable, vertex_coord, slist, nep_info, io_time);
  }

}


void apply_nep
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  ELAPSED_CPU_TIME cpu_time;

  cpu_time.getElapsed();  // initialize cpu time

  if (use_octree) {
    // create octree
    IJKOCTREE::OCTREE octree(dimension, axis_size);
    octree.SetMinMax(scalar_grid.Scalar());
    check_octree(octree);

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;

    nep_using_datastruct
      (scalar_grid, isotable, merge_data, mcube_time, io_time, octree);
  } 
  else if (use_minmax_regions) {
    MINMAX_DATA minmax;
    set_minmax(scalar_grid, region_length, minmax);

    if (use_snapmc) {
      cerr << "Error.  Minmax regions not implemented with snapMC." << endl;
      exit(31);
    }

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;

    nep_using_datastruct
      (scalar_grid, isotable, merge_data, mcube_time, io_time, minmax);
  }
  else {
    // Marching Cubes algorithm with no preprocessing
    nep_no_datastruct
      (scalar_grid, isotable, merge_data, mcube_time, io_time);
  }

}

// **************************************************
// snapMC
// **************************************************

void snapMC_no_datastruct
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, SCALAR_TYPE isovalue,
 vector<VERTEX_INDEX> & simplex_vert, vector<COORD_TYPE> & vertex_coord, 
 NEP_INFO & nep_info)
{
  PROCEDURE_ERROR error("extract_isosurface_snap");

  clock_t t_start = clock();

  SNAP_GRID_CREATOR snap_grid
    (scalar_grid, isovalue, snap_value, nep_info.time.snap);

  if (!use_list) {
    snapMC(scalar_grid, snap_grid, isotable, isovalue,
	   nep_num_dup, simplex_vert, vertex_coord, merge_data, nep_info);
  }
  else {
    clock_t t2 = clock();

    std::vector<VERTEX_INDEX> snap_vlist;
    get_nonempty_snap_cubes(scalar_grid, isotable, isovalue, snap_vlist);

    clock_t t3 = clock();


    snapMC_from_list(scalar_grid, snap_grid, isotable, isovalue,
		     nep_num_dup, snap_vlist, true, simplex_vert,
		     vertex_coord, merge_data, nep_info);

    nep_info.time.preprocessing = float(t3-t2)/CLOCKS_PER_SEC;
  }

  clock_t t_end = clock();

  // store times
  nep_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


void snapMC_no_datastruct
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  NEP_INFO nep_info;
  nep_info.grid.num_cubes = scalar_grid.ComputeNumCubes();

  io_time.write_time = 0;
  for (int i = 0; i < isovalue.size(); i++) {

    snapMC_no_datastruct
      (scalar_grid, isotable, merge_data, isovalue[i],
       slist, vertex_coord, nep_info);
    mcube_time.Add(nep_info.time);

    if (flag_subsample) { rescale_coord(subsample_resolution, vertex_coord); };
    output_isosurface(i, isotable, vertex_coord, slist, nep_info, io_time);
  }

}

void apply_snapMC
(const MC_GRID & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  ELAPSED_CPU_TIME cpu_time;

  cpu_time.getElapsed();  // initialize cpu time

  if (use_octree) {
    cerr << "Error.  Octree not yet implemented with SnapMC." << endl;
    exit(36);
  } 

  if (use_minmax_regions) {
    cerr << "Error.  Minmax regions not yet implemented with SnapMC." << endl;
    exit(36);
  }


  // snapMC algorithm with no preprocessing
  snapMC_no_datastruct
    (scalar_grid, isotable, merge_data, mcube_time, io_time);
}

// **************************************************
// INTERVAL VOLUME EXTRACTION
// **************************************************

void extract_interval_volumes
(const MC_GRID & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract interval volumes from scalar and output to file
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  IVOL_MERGE_DATA merge_data(dimension, axis_size);
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  if (use_octree) {
    cerr << "Error.  Octree not yet implemented with MCVol." << endl;
    exit(36);
  }

  if (use_minmax_regions) {
    cerr << "Error.  Minmax regions not yet implemented with MCVol." << endl;
    exit(37);
  }

  if (use_nep_isotable) {
    cerr << "Error.  NEP isotable not yet implemented with MCVol." << endl;
    exit(38);
  }

  MCUBE_INFO mcube_info;
  mcube_info.grid.num_cubes = compute_num_grid_cubes(dimension, axis_size);

  io_time.write_time = 0;
  for (int i = 0; i+1 < isovalue.size(); i++) {

    MCVol(scalar_grid, isotable, isovalue[i], isovalue[i+1], 
	  slist, vertex_coord, merge_data, mcube_info);
    mcube_time.Add(mcube_info.time);

    if (!nowrite_flag) {
      if (report_time_flag)  
	wall_time.getElapsed();    // set for next call

      if (flag_subsample) 
	{ rescale_coord(subsample_resolution, vertex_coord); };
      write_mesh(i, isotable, vertex_coord, slist);

      if (report_time_flag) {
	io_time.write_time += wall_time.getElapsed();
      };
    };

    if (!use_stdout && !flag_silent) {
      int numv = (vertex_coord.size())/dimension;
      int nums = (slist.size())/isotable.NumVerticesPerSimplex();
      VERTEX_INDEX num_mixed_cubes = mcube_info.grid.num_mixed_cubes;
      VERTEX_INDEX num_grid_cubes = mcube_info.grid.num_cubes;
      float percent = float(num_mixed_cubes)/float(num_grid_cubes);
      int ipercent = int(100*percent);
      cerr << "  Interval [" << isovalue[i] 
	   << ", " << isovalue[i+1] << "].  " 
	   << numv << " volume vertices.  "
	   << nums << " volume simplices." << endl;
      cerr << "    " << num_mixed_cubes 
	   << " (" << ipercent << "%) non-empty cubes." << endl;
    }
  }

}

// **************************************************
// DATA STRUCTURES
// **************************************************

void check_octree(const IJKOCTREE::OCTREE & octree)
{
  if (use_snapmc) {
    cerr << "Error.  Octree not yet implemented with SnapMC." << endl;
    exit(101);
  } 

  if (genivol) {
    cerr << "Error.  Octree not yet implemented with MCVol." << endl;
    exit(102);
  }

  ERROR error;
  if (!octree.Check(error)) { 
    cerr << "Error in constructing octree." << endl;
    error.Print(cerr);
    exit(25);
  }
}

// **************************************************
// SnapMC
// **************************************************

// *** DEBUG ***
void write_snap_nrrd
(const MC_GRID & scalar_grid, const SCALAR_TYPE isovalue,  
 const SNAP_TYPE snap_value)
{
  char * output_filename = "snap.nrrd";

  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  SCALAR_TYPE * scalar_snap = NULL;  // snapped scalar values
  VERTEX_INDEX * snap_toward = NULL; // snap grid vertex iv toward 

  size_t * nrrd_size = new size_t[dimension];
  for (int i = 0; i < dimension; i++)
    nrrd_size[i] = axis_size[i];

  // create grid of snapped scalar values
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  scalar_snap = new SCALAR_TYPE[num_grid_vertices];
  snap_toward = new VERTEX_INDEX[num_grid_vertices];
  snap_scalar_values
    (scalar_grid, isovalue, snap_value, scalar_snap, snap_toward);

  // create nrrd data structure
  Nrrd *nval;

  nval = nrrdNew();
  nrrdWrap_nva(nval, scalar_snap, nrrdTypeFloat, dimension, nrrd_size);
  nrrdSave(output_filename, nval, NULL);
  nrrdNix(nval);

  delete [] scalar_snap;
  scalar_snap = NULL;
  delete [] snap_toward;
  snap_toward = NULL;
  delete [] nrrd_size;
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
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  AXIS_SIZE_TYPE num_region_edges[dimension];

  for (int d = 0; d < dimension; d++) {
    num_region_edges[d] = region_length;
  }

  VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();

  minmax.Create(dimension, scalar_grid.Scalar(), axis_size, num_region_edges);
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

