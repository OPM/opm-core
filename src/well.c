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

#include <stdlib.h>

#include "well.h"


/* Release memory resources for cell->well mapping. */
/* ---------------------------------------------------------------------- */
void
deallocate_cell_wells(int *cwpos, int *cwells)
/* ---------------------------------------------------------------------- */
{
    free(cwells);
    free(cwpos);
}


/* Allocate memory resources for cell->well mapping.
 *
 * Returns 1 if successful and 0 if not.  CSR array pair set to NULL
 * unless allocation succeeds. */
/* ---------------------------------------------------------------------- */
int
allocate_cell_wells(int nc, well_t *W, int **cwpos, int **cwells)
/* ---------------------------------------------------------------------- */
{
    int i, totwconn;

    totwconn = W->well_connpos[W->number_of_wells];

    *cwpos  = malloc((nc + 1)     * sizeof **cwpos );
    *cwells = malloc(2 * totwconn * sizeof **cwells);

    if ((*cwpos == NULL) || (*cwells == NULL)) {
        deallocate_cell_wells(*cwpos, *cwells);

        *cwpos   = NULL;
        *cwells  = NULL;

        totwconn = 0;
    } else {
        for (i = 0; i < nc + 1; i++) {
            (*cwpos)[i] = 0;
        }
    }

    return totwconn;
}


/* Derive cell->well mapping from well->cell (connection) mapping. */
/* ---------------------------------------------------------------------- */
void
derive_cell_wells(int nc, well_t *W, int *cwpos, int *cwells)
/* ---------------------------------------------------------------------- */
{
    int i, w, *c, *connpos;

    connpos = W->well_connpos;

    c = W->well_cells;
    for (w = 0; w < W->number_of_wells; w++) {
        for (; c != W->well_cells + connpos[w + 1]; c++) {
            cwpos[*c + 1] += 1;
        }
    }

    for (i = 1; i <= nc; i++) {
        cwpos[0] += cwpos[i];
        cwpos[i]  = cwpos[0] - cwpos[i];
    }

    cwpos[0] = 0;
    c        = W->well_cells;
    for (w = 0; w < W->number_of_wells; w++) {
        for (; c != W->well_cells + connpos[w + 1]; c++) {
            cwells[ 2*cwpos[*c + 1] + 0 ] = w;
            cwells[ 2*cwpos[*c + 1] + 1 ] = c - W->well_cells;

            cwpos[*c + 1] += 1;
        }
    }
}
