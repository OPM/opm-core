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

#ifndef OPM_FSH_HEADER_INCLUDED
#define OPM_FHS_HEADER_INCLUDED

#include "fsh_common.h"
#include "grid.h"
#include "well.h"
#include "flow_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

struct fsh_data;

/** Constructs incompressible hybrid flow-solver data object for a
 *  given grid and well pattern.
 */
struct fsh_data *
fsh_construct(grid_t *G, well_t *W);



/** Assembles the hybridized linear system for face pressures.
 */
void
fsh_assemble(flowbc_t        *bc,
             const double    *src,
             const double    *Binv,
             const double    *Biv,
             const double    *P,
             const double    *gpress,
             well_control_t  *wctrl,
             const double    *WI,
             const double    *BivW,
             const double    *wdp,
             struct fsh_data *h);

/** Computes cell pressures, face fluxes, well pressures and well
 * fluxes from face pressures.
 */
void
fsh_press_flux(grid_t *G,
               const double *Binv, const double *gpress,
               struct fsh_data *h,
               double *cpress, double *fflux,
               double *wpress, double *wflux);

#ifdef __cplusplus
}
#endif


#endif  /* OPM_FSH_HEADER_INCLUDED */
