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
#include <assert.h>
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
void extract_mesh(const MC_GRID_BASE & scalar_grid, 
		  string & isotable_filename, MCUBE_TIME & mcube_time,
		  IO_TIME & io_time);
void apply_mcube
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);
void apply_nep
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);
void apply_snapMC
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);


// interval volume extraction routines
void extract_interval_volumes
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable, 
 MCUBE_TIME & mcube_time, IO_TIME & io_time);

// routines
void parse_command_line(int argc, char **argv);
void read_isosurface_table
(const int dimension, ISOSURFACE_TABLE & isotable, 
 string & isotable_filename, IO_TIME & io_time);
void check_octree(const IJKOCTREE::OCTREE & octree);
void report_num_cubes(const MCUBE_INFO & mcube_info);
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
void output_isosurface
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const SNAP_INFO & snap_info, IO_TIME & io_time);



// *** DEBUG ***
void write_snap_nrrd
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,  
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

  ELAPSED_TIME wall_time;

  MCUBE_TIME mcube_time;
  IO_TIME io_time = {0.0, 0.0, 0.0};

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
    AXIS_SIZE_TYPE * axis_size = NULL;
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

    axis_size = new AXIS_SIZE_TYPE[dimension];
    size_t size[NRRD_DIM_MAX];

    nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size); 

    for (int d = 0; d < dimension; d++) { axis_size[d] = size[d]; }

    MC_GRID full_scalar_grid(dimension, axis_size);
    nrrd2scalar(nin, full_scalar_grid.ScalarPtr());

    if (report_time_flag) { 
      io_time.read_nrrd_time = wall_time.getElapsed(); 
    };

    VERTEX_INDEX num_grid_cubes = full_scalar_grid.ComputeNumCubes();

    if (!flag_subsample) {
      if (!use_stdout && !flag_silent) 
	{ cerr << num_grid_cubes << " grid cubes." << endl; }

      extract_mesh(full_scalar_grid, isotable_filename, mcube_time, io_time);
    }
    else {
      // subsample grid
      MC_GRID subsampled_grid;
      subsampled_grid.Subsample(full_scalar_grid, subsample_resolution);

      if (!use_stdout && !flag_silent) {
	VERTEX_INDEX num_subsampled_cubes = subsampled_grid.ComputeNumCubes();
	cerr << num_grid_cubes << " grid cubes.  "
	     << num_subsampled_cubes << " subsampled grid cubes." << endl;
      }

      extract_mesh(subsampled_grid, isotable_filename, mcube_time, io_time);
    };

    delete [] axis_size;
    axis_size = NULL;

    nrrdNuke(nin);

    if (report_time_flag) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cerr << endl;
      report_time(io_time, mcube_time, total_elapsed_time,
		  input_filename, isotable_filename.c_str());
    };

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

void extract_mesh(const MC_GRID_BASE & scalar_grid, 
		  string & isotable_filename, MCUBE_TIME & mcube_time,
		  IO_TIME & io_time)
// extract isosurface or interval volumes from scalar grid
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  if (!genivol) {
    if (use_nep_isotable || use_snapmc) {
      NEP_ISOSURFACE_TABLE isotable(dimension);
      read_isosurface_table(dimension, isotable, isotable_filename, io_time);
      set_in_facet_iso_patches(isotable);
      
      NEP_ISO_MERGE_DATA merge_data(dimension, axis_size);
      if (!use_snapmc) {
	apply_nep(scalar_grid, isotable, merge_data, mcube_time, io_time);
      }
      else {
	apply_snapMC(scalar_grid, isotable, merge_data, 
		     mcube_time, io_time);
      }
    }
    else {
      ISOSURFACE_TABLE isotable(dimension);
      read_isosurface_table(dimension, isotable, isotable_filename, io_time);

      ISO_MERGE_DATA merge_data(dimension, axis_size);
      apply_mcube(scalar_grid, isotable, merge_data, mcube_time, io_time);
    };
  }
  else {
    ISOSURFACE_TABLE isotable(dimension);
    read_isosurface_table(dimension, isotable, isotable_filename, io_time);

    extract_interval_volumes
      (scalar_grid, isotable, mcube_time, io_time);
  };
}

// **************************************************
// MARCHING CUBES
// **************************************************

template <class DATASTRUCT> void mcube_using_datastruct
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable, 
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
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable, 
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
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable, 
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
    octree.SetMinMax(scalar_grid.ScalarPtrConst());
    check_octree(octree);

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;

    mcube_using_datastruct(scalar_grid, isotable, merge_data, 
			   mcube_time, io_time, octree);
  } 
  else if (use_minmax_regions) {
    IJKMCUBE::MINMAX_REGIONS minmax;
    minmax.ComputeMinMax(scalar_grid, region_length);

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
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
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
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
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
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
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
    octree.SetMinMax(scalar_grid.ScalarPtrConst());
    check_octree(octree);

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;

    nep_using_datastruct
      (scalar_grid, isotable, merge_data, mcube_time, io_time, octree);
  } 
  else if (use_minmax_regions) {
    IJKMCUBE::MINMAX_REGIONS minmax;
    minmax.ComputeMinMax(scalar_grid, region_length);


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

template <class DATASTRUCT> void snapMC_using_datastruct
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time,
 const DATASTRUCT & datastruct)
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  SNAP_INFO snap_info;
  snap_info.grid.num_cubes = scalar_grid.ComputeNumCubes();

  io_time.write_time = 0;
  for (int i = 0; i < isovalue.size(); i++) {

    snapMC(scalar_grid, isotable, isovalue[i], snap_value, nep_num_dup,
	   slist, vertex_coord, merge_data, datastruct, snap_info);
    mcube_time.Add(snap_info.time);

    if (flag_subsample) 
      { rescale_coord(subsample_resolution, vertex_coord); };
    output_isosurface(i, isotable, vertex_coord, slist, snap_info, io_time);
  }
}

void snapMC_no_datastruct
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, SCALAR_TYPE isovalue,
 vector<VERTEX_INDEX> & simplex_vert, vector<COORD_TYPE> & vertex_coord, 
 SNAP_INFO & snap_info)
{
  clock_t t_start = clock();

  SNAP_GRID snap_grid
    (scalar_grid, isovalue, snap_value, snap_info.time.snap);

  if (!use_list) {
    snapMC(scalar_grid, snap_grid, isotable, isovalue,
	   nep_num_dup, simplex_vert, vertex_coord, merge_data, snap_info);
  }
  else {
    clock_t t2 = clock();

    std::vector<VERTEX_INDEX> snap_vlist;
    get_nonempty_snap_cubes(scalar_grid, isotable, isovalue, snap_vlist);

    clock_t t3 = clock();


    snapMC_from_list(scalar_grid, snap_grid, isotable, isovalue,
		     nep_num_dup, snap_vlist, true, simplex_vert,
		     vertex_coord, merge_data, snap_info);

    snap_info.time.preprocessing = float(t3-t2)/CLOCKS_PER_SEC;
  }

  clock_t t_end = clock();

  // store times
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


void snapMC_no_datastruct
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> slist;
  ELAPSED_TIME wall_time;

  SNAP_INFO snap_info;
  snap_info.grid.num_cubes = scalar_grid.ComputeNumCubes();

  io_time.write_time = 0;
  for (int i = 0; i < isovalue.size(); i++) {

    snapMC_no_datastruct
      (scalar_grid, isotable, merge_data, isovalue[i],
       slist, vertex_coord, snap_info);
    mcube_time.Add(snap_info.time);

    if (flag_subsample) { rescale_coord(subsample_resolution, vertex_coord); };
    output_isosurface(i, isotable, vertex_coord, slist, snap_info, io_time);
  }

}

void apply_snapMC
(const MC_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable, 
 MERGE_DATA & merge_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
// extract isosurfaces from scalar and output to file
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  ELAPSED_CPU_TIME cpu_time;

  cpu_time.getElapsed();  // initialize cpu time

  if (use_octree) {
    SNAP_OCTREE octree(dimension, axis_size);
    octree.SetMinMax(scalar_grid.ScalarPtrConst());

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;
    
    snapMC_using_datastruct(scalar_grid, isotable, merge_data, 
			    mcube_time, io_time, octree);
  }
  else if (use_minmax_regions) {

    if (region_length < 2) {
      cerr << "Region length must be at least 2 with SnapMC." << endl;
      exit(38);
    }

    SNAP_MINMAX minmax;
    minmax.ComputeMinMax(scalar_grid, region_length);

    mcube_time.preprocessing = 
      float(cpu_time.getElapsed())/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;

    snapMC_using_datastruct(scalar_grid, isotable, merge_data, 
			    mcube_time, io_time, minmax);
  }
  else {
    // snapMC algorithm with no preprocessing
    snapMC_no_datastruct
      (scalar_grid, isotable, merge_data, mcube_time, io_time);
  }
}

// **************************************************
// INTERVAL VOLUME EXTRACTION
// **************************************************

void extract_interval_volumes
(const MC_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable, 
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
(const MC_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,  
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


void output_nep_isosurface
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

void output_isosurface
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const NEP_INFO & nep_info, IO_TIME & io_time)
{
  output_nep_isosurface(i, isotable, vertex_coord, slist, nep_info, io_time);
}

void output_isosurface
(const int i, const ISOSURFACE_TABLE & isotable, 
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const SNAP_INFO & snap_info, IO_TIME & io_time)
{
  output_nep_isosurface(i, isotable, vertex_coord, slist, snap_info, io_time);

  if (!use_stdout && !flag_silent) {
    cerr << "    " << snap_info.num_snapped_iso_vertices
	 << " snapped isosurface vertices." << endl;
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
  if (argc == 1) { usage_error(); };

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {
    PARAMETER param = get_parameter_token(argv[iarg]);
    if (param == UNKNOWN_PARAM) break;

    switch(param) {
    case REGION_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &region_length);
      use_minmax_regions = true;
      break;

    case OCTREE_PARAM:
      use_octree = true;
      break;

    case LIST_PARAM:
      use_list = true;
      break;

    case NEP_PARAM:
      // use isotable which differentiates between scalar values less
      //   than the isovalue (negative), equals, and greater than
      //   the isovalue (positive)
      use_nep_isotable = true;  
      break;

    case NEP_DUP_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &nep_num_dup);
      break;

    case IVOL_PARAM:
      genivol = true;
      break;

    case SNAP_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%f", &snap_value);
      use_snapmc = true;
      break;

    case SUBSAMPLE_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &subsample_resolution);
      flag_subsample = true;
      break;

    case DIR_PARAM: 
      iarg++;
      if (iarg >= argc) usage_error();
      isotable_directory = argv[iarg];
      break;

    case OFF_PARAM:
      output_format = OFF;
      break;

    case IV_PARAM:
      output_format = IV;
      break;

    case OUTPUT_FILENAME_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      output_filename = argv[iarg];
      break;

    case STDOUT_PARAM:
      use_stdout = true;
      break;

    case NOWRITE_PARAM:
      nowrite_flag = true;
      break;

    case SILENT_PARAM:
      flag_silent = true;
      break;

    case TIME_PARAM:
      report_time_flag = true;
      break;

    case HELP_PARAM:
      help();
      break;
    };

    iarg++;
  };

  // remaining parameters should be list of isovalues followed
  // by input file name

  // check for more parameter tokens
  for (int j = iarg; j < argc; j++) {
    if (get_parameter_token(argv[j]) != UNKNOWN_PARAM) {
      // argv[iarg] is not an isovalue
      cerr << "Error. Illegal parameter: " << argv[iarg] << endl;
      usage_error();
    }
  }

  if (iarg+2 > argc) {
    cerr << "Error.  Missing input isovalue or input file name." << endl;
    usage_error();
  };

  // store isovalues
  for (int j = iarg; j+1 < argc; j++) {
    isovalue_string.push_back(argv[j]);
    SCALAR_TYPE value;

    istringstream input_string(argv[j]);
    input_string >> value;

    if (input_string.fail()) {
      cerr << "Error. \"" << argv[j] << "\" is not a valid input isovalue." 
	   << endl;
      usage_error();
    };

    isovalue.push_back(value);
  }

  input_filename = argv[argc-1];

  if (use_octree && use_minmax_regions) {
    cerr << "Error.  Can't use both -region and -octree parameters.";
    usage_error();
  }

  if (genivol && isovalue.size() < 2) {
    cerr << "Error.  Need at least two isovalues for interval volume generation." << endl;
    usage_error();
  }

  if (flag_subsample && subsample_resolution <= 1) {
    cerr << "Error.  Subsample resolution must be an integer greater than 1."
	 << endl;
    usage_error();
  };

  if (output_filename != NULL && use_stdout) {
    cerr << "Error.  Can't use both -o and -stdout parameters."
	 << endl;
    usage_error();
  };

  if (use_snapmc && (snap_value < 0.0 || snap_value > 0.5)) {
    cerr << "Error.  Illegal snap value " << snap_value << "."
	 << endl;
    cerr << "        Snap value must be in range [0.0, 0.5]." << endl;
    usage_error();
  };
}

void usage_msg()
{
  cerr << "Usage: ijkmcube [-octree] [-region L] [-nep] [-snap D] [-ivol] [-list] [-subsample S] [-help] [-off|-iv] [-dir {isotable_directory}] [-o {output_filename}] [-stdout] [-nowrite] [-time] {isovalue1 isovalue2 ...} {input filename}" 
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
  cerr << endl;
  cerr << "ijkmcube - Marching cubes isosurface generation algorithm." << endl;
  cerr << endl;
  cerr << "Options:" << endl;

  cerr << "  -octree: Create octree for faster isosurface extraction."
       << endl;
  cerr << "  -region L: Preprocess into regions of dimension LxLx... for faster"
       << endl;
  cerr << "              isosurface extraction.  Each region has a minimum and maximum" << endl;
  cerr << "              scalar value."
       << endl;
  cerr << "  -nep:  Use isotable which differentiates scalar values less than (negative)," << endl;
  cerr << "         equal to, and greater than (positive) the isovalue." << endl;
  cerr << "  -snap D: Snap isosurface vertices within distance D of grid vertices." << endl;
  cerr << "           Value D must be in range [0.0, 0.5]." << endl;
  cerr << "  -ivol: Generate interval volume." << endl;
  cerr << "  -subsample S: Subsample grid at every S vertices." << endl;
  cerr << "                S must be an integer greater than 1." << endl;
  cerr << "  -list:  Preprocess by creating list of mixed cubes." << endl;
  cerr << "  -off: Output in geomview OFF format. (Default.)" << endl;
  cerr << "  -iv: Output in OpenInventor .iv format." << endl;
  cerr << "  -dir {isotable_directory}: Directory containing appropriate isosurface table." << endl;
  cerr << "  -o {output_filename}: Write isosurface to file {output_filename}." << endl;
  cerr << "  -stdout: Write isosurface to standard output." << endl;
  cerr << "  -nowrite: Don't write isosurface." << endl;
  cerr << "  -time: Output running time." << endl;
  cerr << "  -s: Silent mode." << endl;
  cerr << "  -help: Print this help message." << endl;
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

