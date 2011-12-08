/*===========================================================================
//
// File: cart_grid.h
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_CART_GRID_H_HEADER
#define OPM_CART_GRID_H_HEADER

#ifdef __cplusplus
extern "C" {
#endif

#include <grid.h>

void
destroy_cart_grid(grid_t *G);

grid_t *
create_cart_grid(int nx, int ny, int nz);

#ifdef __cplusplus
}
#endif
#endif  /* OPM_CART_GRID_H_HEADER */
