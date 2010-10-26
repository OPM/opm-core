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

#ifndef OPM_FSH_COMMON_IMPL_HEADER_INCLUDED
#define OPM_FSH_COMMON_IMPL_HEADER_INCLUDED

/* Internal header.  Don't install. */

struct fsh_impl {
    int nc, nf, nw;             /* Number of cells, faces, wells */

    /* Topology */
    int           *gdof_pos;    /* Pointers, grid DOFs (== cell_facepos) */
    int           *gdof;        /* Grid DOFs           (== cell_faces) */

    int           *cwell_pos;   /* Start pointers, well DOFs (c->w) */
    int           *cwells;      /* Well DOFs (c->w) */

    /* Discretisation data */
    double        *WI;          /* Permuted well production indices */
    double        *wdp;         /* Permuted well gravity pressures */

    double        *cflux;       /* Cell (half-contact) fluxes */

    struct hybsys      *sys;    /* Hybrid cell contribs */
    struct hybsys_well *wsys;   /* Hybrid cell contribs from wells */

    double *work;               /* Scratch array, floating point */
    int    *iwork;              /* Scratch array, integers */

    /* Linear storage goes here... */
    int    *idata;              /* Actual storage array, integers */
    double *ddata;              /* Actual storage array, floating point */
};


struct fsh_impl *
fsh_impl_allocate_basic(size_t idata_sz, size_t ddata_sz);

void
fsh_count_grid_dof(grid_t *G, int *max_ngdof, size_t *sum_ngdof2);

int
fsh_impose_bc(int              ndof,
              int             *dof,
              flowbc_t        *bc,
              struct fsh_impl *pimpl);

void
fsh_define_impl_arrays(size_t           nc,
                       size_t           nnu,
                       size_t           nhf,
                       size_t           max_ncf,
                       well_t          *W,
                       struct fsh_impl *pimpl);

void
fsh_define_cell_wells(size_t nc, well_t *W, struct fsh_impl *pimpl);

void
fsh_define_linsys_arrays(struct fsh_data *h);

#endif /* OPM_FSH_COMMON_IMPL_HEADER_INCLUDED */
