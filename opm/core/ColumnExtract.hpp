#include <opm/core/grid.h>
#include <vector>
#include <map>
#include <algorithm>

namespace Opm {

/// Helper struct for extractColumn
/// Compares the underlying k-index
struct ExtractColumnCompare
{
    ExtractColumnCompare(const UnstructuredGrid& grid)
        : grid(grid)
    {
        // empty
    }

    bool operator()(const int i, const int j)
    {
        // Extract k-index
        int k_i = grid.global_cell[i] / grid.cartdims[0] / grid.cartdims[1];
        int k_j = grid.global_cell[j] / grid.cartdims[0] / grid.cartdims[1];

        return k_i <= k_j;
    }

    const UnstructuredGrid& grid;
};

/// Extract each column of the grid.
///  \note Assumes the pillars of the grid are all vertically aligned.
///  \param grid The grid from which to extract the columns.
///  \param columns will for each i + cartgrim[0]*j where (i, j) represents a non-empty column,
////        contain the cell indices contained in the column
///         centered at (i, j).
void extractColumn( const UnstructuredGrid& grid, std::map<int, std::vector<int> >& columns )
{
    const int* dims = grid.cartdims;
    const int* global = grid.global_cell;
    for (int i = 0; i < grid.number_of_cells; ++i) {
        // Extract Cartesian coordinates
        int index = grid.global_cell[i];
        int i_cart = index % dims[0];
        int k_cart = index / dims[0] / dims[1];
        int j_cart = (index - k_cart*dims[0]*dims[1])/ dims[0];

        // Finally add it to the map
        columns[i_cart + j_cart*dims[0]].push_back(i);
    }

    for (std::map<int, std::vector<int> >::iterator it = columns.begin(); it != columns.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), ExtractColumnCompare(grid));
    }
}

} // namespace Opm
