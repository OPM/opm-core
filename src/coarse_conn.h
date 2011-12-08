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

#ifndef OPM_COARSE_CONN_HEADER_INCLUDED
#define OPM_COARSE_CONN_HEADER_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

struct coarse_topology {
    int nblocks;                /* Number of blocks */
    int nfaces;                 /* Number of faces */

    int *neighbours;            /* Neighbourship definition */

    int *blkfacepos;            /* Index into blkfaces */
    int *blkfaces;              /* Faces per block */

    int *subfacepos;            /* Index into subfaces */
    int *subfaces;              /* FS faces per coarse face */
};


struct coarse_topology *
coarse_topology_create(int nc, int nf, int expct_nconn,
                       const int *p, const int *neighbours);


void
coarse_topology_destroy(struct coarse_topology *t);


#ifdef __cplusplus
}
#endif

#endif  /* OPM_COARSE_CONN_HEADER_INCLUDED */
