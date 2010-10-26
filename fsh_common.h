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

#ifndef OPM_FSH_COMMON_HEADER_INCLUDED
#define OPM_FSH_COMMON_HEADER_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

struct CSRMatrix;
struct fsh_impl;

/** Contains the linear system for assembly, as well as internal data
 *  for the assembly routines.
 */
struct fsh_data {
    /* Let \f$n_i\f$ be the number of connections/faces of grid cell
     * number \f$i\f$. Then max_ngconn = \f$\max_i n_i\f$
     */
    int               max_ngconn;
    /* With n_i as above, sum_ngconn2 = \f$\sum_i n_i^2\f$ */
    size_t            sum_ngconn2;

    /* Linear system */
    struct CSRMatrix *A;        /* Coefficient matrix */
    double           *b;        /* System RHS */
    double           *x;        /* Solution */

    /* Private implementational details. */
    struct fsh_impl  *pimpl;
};



/** Destroys the fsh data object */
void
fsh_destroy(struct fsh_data *h);

#ifdef __cplusplus
}
#endif


#endif  /* OPM_FSH_COMMON_HEADER_INCLUDED */
