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

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/pressure/fsh.h>
#include <opm/core/pressure/fsh_common_impl.h>
#include <opm/core/pressure/mimetic/hybsys.h>
#include <opm/core/pressure/mimetic/hybsys_global.h>




/* ---------------------------------------------------------------------- */
/* Compute cell pressures (cpress) and interface fluxes (fflux) from
 * current solution of system of linear equations, h->x.  Back
 * substitution process, projected half-contact fluxes. */
/* ---------------------------------------------------------------------- */
void
fsh_press_flux(grid_t *G,
               const double *Binv, const double *gpress,
               struct fsh_data *h,
               double *cpress, double *fflux,
               double *wpress, double *wflux)
/* ---------------------------------------------------------------------- */
{
    int c, f, i;
    double s;

    hybsys_compute_press_flux(G->number_of_cells,
                              G->cell_facepos,
                              G->cell_faces,
                              gpress, Binv,
                              h->pimpl->sys,
                              h->x, cpress, h->pimpl->cflux,
                              h->pimpl->work);

    if (h->pimpl->nw > 0) {
        assert ((wpress != NULL) && (wflux != NULL));
        hybsys_compute_press_flux_well(G->number_of_cells, G->cell_facepos,
                                       G->number_of_faces, h->pimpl->nw,
                                       h->pimpl->cwell_pos, h->pimpl->cwells,
                                       Binv, h->pimpl->WI,
                                       h->pimpl->wdp, h->pimpl->sys,
                                       h->pimpl->wsys, h->x, cpress,
                                       h->pimpl->cflux, wpress, wflux,
                                       h->pimpl->work);
    }

    for (f = 0; f < G->number_of_faces; f++) { fflux[f] = 0.0; }

    i = 0;
    for (c = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];
            s = 2.0*(G->face_cells[2*f + 0] == c) - 1.0;

            fflux[f] += s * h->pimpl->cflux[i];
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        i = (G->face_cells[2*f + 0] >= 0) +
            (G->face_cells[2*f + 1] >= 0);

        fflux[f] /= i;
    }
}
