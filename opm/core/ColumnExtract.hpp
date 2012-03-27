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
	int index_i = grid.global_cell ? grid.global_cell[i] : i;
        int k_i = index_i / grid.cartdims[0] / grid.cartdims[1];
	int index_j = grid.global_cell ? grid.global_cell[j] : j;
        int k_j = index_j / grid.cartdims[0] / grid.cartdims[1];

        return k_i < k_j;
    }

    const UnstructuredGrid& grid;
};

/// Extract each column of the grid.
///  \note Assumes the pillars of the grid are all vertically aligned.
///  \param grid The grid from which to extract the columns.
///  \param columns will for each (i, j) where (i, j) represents a non-empty column,
////        contain the cell indices contained in the column
///         centered at (i, j) in the second variable, and i+jN in the first variable.
void extractColumn( const UnstructuredGrid& grid, std::pair<std::vector<int>, std::vector<std::vector<int> > >& columns )
{
    const int* dims = grid.cartdims;
    
    // Keeps track of column_index ---> index of vector
    std::map<int, int> global_to_local;
    for (int i = 0; i < grid.number_of_cells; ++i) {
        // Extract Cartesian coordinates
        int index = grid.global_cell ? grid.global_cell[i] : i; // If null, assume mapping is identity.
        int i_cart = index % dims[0];
        int k_cart = index / dims[0] / dims[1];
        int j_cart = (index - k_cart*dims[0]*dims[1])/ dims[0];

        int local_index;
        std::map<int, int>::iterator local_index_iterator = global_to_local.find(i_cart+j_cart*dims[0]);
        if( local_index_iterator != global_to_local.end() ) {
            local_index = local_index_iterator->second;
        } else {
            local_index = columns.second.size();
            global_to_local[i_cart+j_cart*dims[0]] = local_index;            
            columns.second.push_back(std::vector<int>());
        }
        
        columns.first.push_back(i_cart + j_cart*dims[0]);
        columns.second[local_index].push_back(i);
    }

    int num_cols = columns.second.size();
    for(int i = 0; i < num_cols; ++i) {
        std::sort(columns.second[i].begin(), columns.second[i].end(), ExtractColumnCompare(grid));
    }
}

} // namespace Opm
