#include <opm/core/grid.h>
#include <vector>

namespace Opm {
/**
 * Extract each column of the grid.
 * \note Assumes the pillars of the grid are all vertically aligned.
 * \param grid The grid from which to extract the columns.
 * \param columns will for each i + cartgrim[0]*j  contain the k values contained in the column
 *        centered at (i, j).
 */
void extractColumn( const UnstructuredGrid& grid, std::vector<std::vector<int> >& columns ) {

    columns.resize(grid.cartdims[0]*grid.cartdims[1]);

    for(int i = 0; i < grid.cartdims[0]; i++) {
        for(int j = 0; j < grid.cartdims[1]; j++) {
            int plane_index = i + j*grid.cartdims[0];
            for(int k = 0; k < grid.cartdims[2]; k++) {
                columns[plane_index].push_back(plane_index + k*grid.cartdims[0]*grid.cartdims[1]);
            }
        }
    }



}

} // namespace Opm
