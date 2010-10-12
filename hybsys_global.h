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

#ifndef HYBSYS_GLOBAL_H_INCLUDED
#define HYBSYS_GLOBAL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include "grid.h"
#include "well.h"
#include "sparse_sys.h"

struct CSRMatrix *
hybsys_define_globconn(grid_t *G, well_t *W);


void
hybsys_global_assemble_cell(int nconn, int *l2g,
                            const double     *S,
                            const double     *r,
                            struct CSRMatrix *A,
                            double           *b);

void
hybsys_global_assemble_well_sym(int ngconn_tot,
                                int ngconn, const int *gconn,
                                int nwconn, const int *wconn,
                                const double     *r2w,
                                const double     *w2w,
                                const double     *r,
                                struct CSRMatrix *A,
                                double           *b);



#ifdef __cplusplus
}
#endif

#endif  /* HYBSYS_GLOBAL_H_INCLUDED */
