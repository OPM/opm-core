/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_TRANS_TPFA_HEADER_INCLUDED
#define OPM_TRANS_TPFA_HEADER_INCLUDED

#include "grid.h"

#ifdef __cplusplus
extern "C" {
#endif

void
tpfa_htrans_compute(grid_t *G, const double *perm, double *htrans);

void
tpfa_trans_compute(grid_t *G, const double *htrans, double *trans);

void
tpfa_eff_trans_compute(grid_t       *G,
                       const double *totmob,
                       const double *htrans,
                       double       *trans);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_TRANS_TPFA_HEADER_INCLUDED */
