/*
 * test_read_grid.c
 *
 *  Created on: Aug 28, 2012
 *      Author: bska
 */

#include <stdio.h>

#include <opm/core/grid.h>
#include <opm/core/grid/cart_grid.h>

int
main(void)
{
    struct UnstructuredGrid *G1, *G2;

    G1 = read_grid("cart_grid_2d.txt");
    G2 = create_grid_cart2d(2, 2);

    destroy_grid(G2);
    destroy_grid(G1);

    return 0;
}
