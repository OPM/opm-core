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

#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/pressure/tpfa/compr_quant_general.h>


void
compr_quantities_gen_deallocate(struct compr_quantities_gen *cq)
{
    if (cq != NULL) {
        free(cq->Ac);
    }

    free(cq);
}


struct compr_quantities_gen *
compr_quantities_gen_allocate(size_t nc, size_t nf, int np)
{
    size_t                       alloc_sz, np2;
    struct compr_quantities_gen *cq;

    cq = malloc(1 * sizeof *cq);

    if (cq != NULL) {
        np2 = np * np;

        alloc_sz  = np2 * nc;   /* Ac */
        alloc_sz += np2 * nc;   /* dAc */
        alloc_sz += np2 * nf;   /* Af */
        alloc_sz += np  * nf;   /* phasemobf */
        alloc_sz += 1   * nc;   /* voldiscr */

        cq->Ac = malloc(alloc_sz * sizeof *cq->Ac);

        if (cq->Ac == NULL) {
            compr_quantities_gen_deallocate(cq);
            cq = NULL;
        } else {
            cq->dAc       = cq->Ac        + (np2 * nc);
            cq->Af        = cq->dAc       + (np2 * nc);
            cq->phasemobf = cq->Af        + (np2 * nf);
            cq->voldiscr  = cq->phasemobf + (np  * nf);

            cq->nphases   = np;

            vector_zero(alloc_sz, cq->Ac);
        }
    }

    return cq;
}
