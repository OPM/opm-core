/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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


#include "config.h"
#include <opm/core/grid/cart_grid.h>
#include <opm/core/grid.h>
#include <cstdio>

int main(void)
{
    using namespace std;
    struct UnstructuredGrid *g = create_grid_cart2d(2, 2, 1., 1.);
    int i;
    int k;
    for (i = 0; i < g->number_of_cells; ++i) {
        fprintf(stderr, "%d: ", i);
        for (k = g->cell_facepos[i]; k < g->cell_facepos[i + 1]; ++k) {
            fprintf(stderr, "%d ", g->cell_faces[k]);
        }
        fprintf(stderr, "\n");
    }
    destroy_grid(g);
    return 0;
}
