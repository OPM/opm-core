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
    int totwconn;

    totwconn = W->well_connpos[W->number_of_wells];

    *cwpos  = calloc(nc + 1       , sizeof **cwpos );
    *cwells = malloc(2 * totwconn * sizeof **cwells);

    if ((*cwpos == NULL) || (*cwells == NULL)) {
        deallocate_cell_wells(*cwpos, *cwells);

        *cwpos   = NULL;
        *cwells  = NULL;

        totwconn = 0;
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
