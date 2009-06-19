//===========================================================================
//
// File: uniquepoints.h
//
// Created: Fri Jun 19 08:46:30 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
Copyright 2009 SINTEF ICT, Applied Mathematics.
Copyright 2009 Statoil ASA.

This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_UNIQUEPOINTS_HEADER
#define OPENRS_UNIQUEPOINTS_HEADER

int finduniquepoints(const struct grdecl *g,  /* input */
		     int                 *p,  /* for each z0 in zcorn, z0 = z[p0] */
		     sparse_table_t      *z,  /* list of uniq zcorn valules for each pillar*/
		     double               t); /* tolerance*/

#endif // OPENRS_UNIQUEPOINTS_HEADER

