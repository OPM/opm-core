#include <opm/core/grid.h>
#include <vector>

/**
 * Extract each column of the grid.
 * \note Assumes the pillars of the grid are all vertically aligned.
 * \param grid The grid from which to extract the columns.
 * \param columns will for each i + cartgrim[0]*j  contain the k values contained in the column
 *        centered at (i, j).
 */
void extractColumn( const UnstructuredGrid& grid, std::vector<std::vector<int> >& columns ) {

    columns.resize(grid.cartdims[0]*grid.cartdims[1]);

    // This is used for sorting (and discarded afterwards)
    std::vector<std::vector<double> > z_values;
    z_values.resize(columns.size());
    for(int i = 0; i < grid.cartdims[0]; i++) {
        for(int j = 0; j < grid.cartdims[1]; j++) {
            int plane_index = i + j*grid.cartdims[0];
            for(int k = 0; k < grid.cartdims[2]; k++) {
                columns[plane_index].push_back(plane_index + k*grid.cartdims[0]*grid.cartdims[1]);
                z_values[plane_index].push_back(grid.cell_volumes[3]);
            }
        }
    }


}
