/// \file ijkgenpatch.h
/// generate isosurface patch

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006, 2003, 2001 Rephael Wenger

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

#ifndef _GENISOPATCH_
#define _GENISOPATCH_

#include <vector>

#include "isotable.h"

void genisopatch
  (const ISOTABLE::ISOSURFACE_TABLE_POLYHEDRON & polyhedron,
   const int * const vertex_sign, 
   vector<int> & edge_list, int & num_simplices);

void ijkgenpatch_nep
  (const ISOSURFACE_TABLE_POLYHEDRON & polyhedron, 
   const int * const vertex_sign,
   vector<ISOSURFACE_VERTEX_INDEX> & isov_list, 
   vector<ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE> & isov_type,
   int & num_simplices);

#endif
