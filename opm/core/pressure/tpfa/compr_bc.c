/*===========================================================================
//
// File: compr_bc.c
//
// Created: 2011-10-24 16:07:17+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/pressure/tpfa/compr_bc.h>

static int
expand_source_tables(int alloc, struct compr_bc *bc)
{
    enum compr_bc_type *t;
    int                *f;
    double             *p, *v, *s;

    t = realloc(bc->type      , alloc * 1           * sizeof *t);
    f = realloc(bc->face      , alloc * 1           * sizeof *f);
    p = realloc(bc->press     , alloc * 1           * sizeof *p);
    v = realloc(bc->flux      , alloc * 1           * sizeof *v);
    s = realloc(bc->saturation, alloc * bc->nphases * sizeof *s);

    if ((t == NULL) || (f == NULL) ||
        (p == NULL) || (v == NULL) || (s == NULL)) {

        free(s);  free(v);  free(p);  free(f);  free(t);

        alloc = 0;
    } else {
        bc->type       = t;  bc->face = f;
        bc->press      = p;  bc->flux = v;
        bc->saturation = s;
    }

    return alloc;

}


/* ======================================================================
 * Public methods below separator.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
struct compr_bc *
compr_bc_allocate(int np, int nbc)
/* ---------------------------------------------------------------------- */
{
    int               status;
    struct compr_bc *bc;

    if (np <= 0) {
        bc = NULL;
    } else {
        bc = malloc(1 * sizeof *bc);

        if (bc != NULL) {
            bc->nbc     = 0 ;
            bc->cpty    = 0 ;
            bc->nphases = np;

            bc->type       = NULL;
            bc->face       = NULL;
            bc->press      = NULL;
            bc->flux       = NULL;
            bc->saturation = NULL;

            if (nbc > 0) {
                status = expand_source_tables(nbc, bc);
            } else {
                status = 1;
            }

            if (status <= 0) {
                compr_bc_deallocate(bc);
                bc = NULL;
            }
        }
    }

    return bc;
}


/* ---------------------------------------------------------------------- */
void
compr_bc_deallocate(struct compr_bc *bc)
/* ---------------------------------------------------------------------- */
{
    if (bc != NULL) {
        free(bc->saturation);
        free(bc->flux      );
        free(bc->press     );
        free(bc->face      );
        free(bc->type      );
    }

    free(bc);
}


/* ---------------------------------------------------------------------- */
int
compr_bc_append(enum compr_bc_type  t  ,
                int                 f  ,
                int                 np ,
                double              p  ,
                double              v  ,
                const double       *sat,
                struct compr_bc    *bc )
/* ---------------------------------------------------------------------- */
{
    int alloc, status;

    if (np != bc->nphases) {
        return -1;
    }

    if (bc->nbc == bc->cpty) {
        alloc  = (bc->cpty > 0) ? 2 * bc->cpty : 1;
        status = expand_source_tables(alloc, bc);
    } else {
        assert (bc->nbc < bc->cpty);
        status = 1;
    }

    if (status > 0) {
        bc->type [ bc->nbc ] = t;
        bc->face [ bc->nbc ] = f;
        bc->press[ bc->nbc ] = p;
        bc->flux [ bc->nbc ] = v;

        memcpy(bc->saturation + (bc->nbc * np), sat,
               np * sizeof *bc->saturation);

        bc->nbc += 1;
    }

    return status > 0;
}


/* ---------------------------------------------------------------------- */
void
compr_bc_clear(struct compr_bc *bc)
/* ---------------------------------------------------------------------- */
{
    bc->nbc = 0;
}
