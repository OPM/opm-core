/*===========================================================================
//
// File: compr_source.c
//
// Created: 2011-10-19 19:22:24+0200
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

#include "compr_source.h"

static int
expand_source_tables(int alloc, struct compr_src *src)
{
    int    *c;
    double *v, *s;

    c = realloc(src->cell      , alloc * 1            * sizeof *c);
    v = realloc(src->flux      , alloc * 1            * sizeof *v);
    s = realloc(src->saturation, alloc * src->nphases * sizeof *s);

    if ((c == NULL) || (v == NULL) || (s == NULL)) {

        free(s);  free(v);   free(c);

        alloc = 0;
    } else {
        src->cell = c;  src->cpty       = alloc;
        src->flux = v;  src->saturation = s;
    }

    return alloc;

}


/* ======================================================================
 * Public methods below separator.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
struct compr_src *
compr_src_allocate(int np, int nsrc)
/* ---------------------------------------------------------------------- */
{
    int               status;
    struct compr_src *src;

    if (np <= 0) {
        src = NULL;
    } else {
        src = malloc(1 * sizeof *src);

        if (src != NULL) {
            src->nsrc    = 0 ;
            src->cpty    = 0 ;
            src->nphases = np;

            src->cell       = NULL;
            src->flux       = NULL;
            src->saturation = NULL;

            if (nsrc > 0) {
                status = expand_source_tables(nsrc, src);
            } else {
                status = 1;
            }

            if (status <= 0) {
                compr_src_deallocate(src);
                src = NULL;
            }
        }
    }

    return src;
}


/* ---------------------------------------------------------------------- */
void
compr_src_deallocate(struct compr_src *src)
/* ---------------------------------------------------------------------- */
{
    if (src != NULL) {
        free(src->saturation);
        free(src->flux      );
        free(src->cell      );
    }

    free(src);
}


/* ---------------------------------------------------------------------- */
int
append_compr_source_term(int               c  ,
                         int               np ,
                         double            v  ,
                         const double     *sat,
                         struct compr_src *src)
/* ---------------------------------------------------------------------- */
{
    int alloc, status;

    if (np != src->nphases) {
        return -1;
    }

    if (src->nsrc == src->cpty) {
        alloc  = (src->cpty > 0) ? 2 * src->cpty : 1;
        status = expand_source_tables(alloc, src);
    } else {
        assert (src->nsrc < src->cpty);
        status = 1;
    }

    if (status > 0) {
        src->cell[ src->nsrc ] = c;
        src->flux[ src->nsrc ] = v;

        memcpy(src->saturation + (src->nsrc * np), sat,
               np * sizeof *src->saturation);

        src->nsrc += 1;
    }

    return status > 0;
}


/* ---------------------------------------------------------------------- */
void
clear_compr_source_term(struct compr_src *src)
/* ---------------------------------------------------------------------- */
{
    src->nsrc = 0;
}
