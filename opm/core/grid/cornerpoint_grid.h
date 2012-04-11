/*===========================================================================
//
// File: preprocess.h
//
// Created: Fri Jun 19 08:43:04 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date: 2010-08-27 19:12:16 +0200 (Fri, 27 Aug 2010) $
//
// $Revision: 930 $
//
//==========================================================================*/

/*
  Copyright 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_CORNERPOINT_GRID_HEADER_INCLUDED
#define OPM_CORNERPOINT_GRID_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/grid/cpgpreprocess/preprocess.h>

#ifdef __cplusplus
extern "C" {
#endif

    struct UnstructuredGrid *
    preprocess (const struct grdecl *in, double tol);

    void compute_geometry     (struct UnstructuredGrid *g);
    
    void free_grid(struct UnstructuredGrid *g);
#ifdef __cplusplus
}
#endif

#endif /* OPM_CORNERPOINT_GRID_HEADER_INCLUDED */
