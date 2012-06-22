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

#include <opm/core/pressure/tpfa/compr_quant.h>


/* ---------------------------------------------------------------------- */
/* Compute B \ (V') == zeta(cellNo) .* faceFlux2CellFlux(fflux) */
/* ---------------------------------------------------------------------- */
void
compr_flux_term(struct UnstructuredGrid       *G,
                const double *fflux,
                const double *zeta,
                double       *Biv)
/* ---------------------------------------------------------------------- */
{
    int    c, i, f;
    double s;

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];
            s = 2.0*(c == G->face_cells[2*f + 0]) - 1.0;

            Biv[i] = zeta[c] * s * fflux[f];
        }
    }
}


/* ---------------------------------------------------------------------- */
/* Compute P == ct .* pv ./ dt */
/* ---------------------------------------------------------------------- */
void
compr_accum_term(size_t        nc,
                 double        dt,
                 const double *porevol,
                 const double *totcompr,
                 double       *P)
/* ---------------------------------------------------------------------- */
{
    size_t c;

    for (c = 0; c < nc; c++) {
        P[c] = totcompr[c] * porevol[c] / dt;
    }
}


/* ---------------------------------------------------------------------- */
/* Add pressure accumulation term (P*p_0) to cell sources */
/* ---------------------------------------------------------------------- */
void
compr_src_add_press_accum(size_t        nc,
                          const double *p0,
                          const double *P,
                          double       *src)
/* ---------------------------------------------------------------------- */
{
    size_t c;

    for (c = 0; c < nc; c++) {
        src[c] += P[c] * p0[c];
    }
}
