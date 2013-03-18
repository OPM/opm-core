/*===========================================================================
//
// File: transport_source.c
//
// Created: 2011-10-05 19:58:38+0200
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

#include <stdlib.h>
#include <string.h>

#include <opm/core/transport/implicit/transport_source.h>


/* ---------------------------------------------------------------------- */
static int
expand_source_tables(int alloc, struct TransportSource *src)
/* ---------------------------------------------------------------------- */
{
    int    *c;
    double *p, *v, *s, *z;

    c = realloc(src->cell      , alloc * 1           * sizeof *c);
    p = realloc(src->pressure  , alloc * 1           * sizeof *p);
    v = realloc(src->flux      , alloc * 1           * sizeof *v);
    s = realloc(src->saturation, alloc * src->nphase * sizeof *s);
    z = realloc(src->surfvolume, alloc * src->nphase * sizeof *z);

    if ((c == NULL) ||
        (p == NULL) || (v == NULL) ||
        (s == NULL) || (z == NULL)) {

        free(z);  free(s);  free(v);  free(p);   free(c);

        alloc = 0;
    } else {
        src->cell       = c;  src->cpty       = alloc;
        src->pressure   = p;  src->flux       = v;
        src->saturation = s;  src->surfvolume = z;
    }

    return alloc;
}


/* ======================================================================
 * Public methods below separator.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
struct TransportSource *
create_transport_source(int nsrc, int nphase)
/* ---------------------------------------------------------------------- */
{
    int                     status;
    struct TransportSource *src;

    if (nphase <= 0) {
        src = NULL;
    } else {
        src = malloc(1 * sizeof *src);

        if (src != NULL) {
            src->nsrc = src->cpty = 0;

            src->nphase     = nphase;
            src->cell       = NULL;
            src->pressure   = NULL;  src->flux       = NULL;
            src->saturation = NULL;  src->surfvolume = NULL;

            if (nsrc > 0) {
                status = expand_source_tables(nsrc, src);
            } else {
                status = 1;
            }

            if (status <= 0) {
                destroy_transport_source(src);
                src = NULL;
            }
        }
    }

    return src;
}


/* ---------------------------------------------------------------------- */
void
destroy_transport_source(struct TransportSource *src)
/* ---------------------------------------------------------------------- */
{
    if (src != NULL) {
        free(src->surfvolume);  free(src->saturation);
        free(src->flux)      ;  free(src->pressure)  ;
        free(src->cell)      ;
    }

    free(src);
}


/* ---------------------------------------------------------------------- */
int
append_transport_source(int                     c,
                        int                     nphase,
                        double                  p,
                        double                  v,
                        const double           *sat,
                        const double           *z,
                        struct TransportSource *src)
/* ---------------------------------------------------------------------- */
{
    int alloc, status;

    if (nphase != src->nphase) {
        return -1;
    }

    if (src->nsrc == src->cpty) {
        alloc  = (src->cpty > 0) ? 2 * src->cpty : 1;
        status = expand_source_tables(alloc, src);
    } else {
        status = 1;
    }

    if (status > 0) {
        src->cell    [ src->nsrc ] = c;
        src->pressure[ src->nsrc ] = p;
        src->flux    [ src->nsrc ] = v;

        memcpy(src->saturation + (src->nsrc * nphase), sat,
               nphase * sizeof *src->saturation);

        memcpy(src->surfvolume + (src->nsrc * nphase),  z ,
               nphase * sizeof *src->surfvolume);

        src->nsrc += 1;
    }

    return status > 0;
}


/* ---------------------------------------------------------------------- */
void
clear_transport_source(struct TransportSource *src)
/* ---------------------------------------------------------------------- */
{
    src->nsrc = 0;
}
