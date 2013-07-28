#include <opm/core/grid.h>
#include <vector>
#include <map>
#include <algorithm>

namespace Opm {

    namespace {

        /// Helper struct for extractColumn
        /// Compares the underlying k-index
        struct ExtractColumnCompare
        {
            ExtractColumnCompare(const UnstructuredGrid& g)
            : grid(g)
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


        /// Neighbourhood query.
        /// \return true if two cells are neighbours.
        bool neighbours(const UnstructuredGrid& grid, const int c0, const int c1)
        {
            for (int hf = grid.cell_facepos[c0]; hf < grid.cell_facepos[c0 + 1]; ++hf) {
                const int f = grid.cell_faces[hf];
                if (grid.face_cells[2*f] == c1 || grid.face_cells[2*f+1] == c1) {
                    return true;
                }
            }
            return false;
        }

    } // anonymous namespace


/// Extract each column of the grid.
///  \note Assumes the pillars of the grid are all vertically aligned.
///  \param grid The grid from which to extract the columns.
///  \param columns will for each (i, j) where (i, j) represents a non-empty column,
////        contain the cell indices contained in the column
///         centered at (i, j) in the second variable, and i+jN in the first variable.
inline void extractColumn( const UnstructuredGrid& grid, std::vector<std::vector<int> >& columns )
{
    const int* dims = grid.cartdims;

    // Keeps track of column_index ---> index of vector
    std::map<int, int> global_to_local;
    for (int cell = 0; cell < grid.number_of_cells; ++cell) {
        // Extract Cartesian coordinates
        int index = grid.global_cell ? grid.global_cell[cell] : cell; // If null, assume mapping is identity.
        int i_cart = index % dims[0];
        int k_cart = index / dims[0] / dims[1];
        int j_cart = (index - k_cart*dims[0]*dims[1])/ dims[0];

        int local_index;
        std::map<int, int>::iterator local_index_iterator = global_to_local.find(i_cart+j_cart*dims[0]);
        if (local_index_iterator != global_to_local.end()) {
            local_index = local_index_iterator->second;
        } else {
            local_index = columns.size();
            global_to_local[i_cart+j_cart*dims[0]] = local_index;
            columns.push_back(std::vector<int>());
        }
        columns[local_index].push_back(cell);
    }

    int num_cols = columns.size();
    for (int col = 0; col < num_cols; ++col) {
        std::sort(columns[col].begin(), columns[col].end(), ExtractColumnCompare(grid));
    }

    // At this point, a column may contain multiple disjoint sets of cells.
    // We must split these columns into connected parts.
    std::vector< std::vector<int> > new_columns;
    for (int col = 0; col < num_cols; ++col) {
        const int colsz = columns[col].size();
        int first_of_col = 0;
        for (int k = 1; k < colsz; ++k) {
            const int c0 = columns[col][k-1];
            const int c1 = columns[col][k];
            if (!neighbours(grid, c0, c1)) {
                // Must split. Move the cells [first_of_col, ... , k-1] to
                // a new column, known to be connected.
                new_columns.push_back(std::vector<int>());
                new_columns.back().assign(columns[col].begin() + first_of_col, columns[col].begin() + k);
                // The working column now starts with index k.
                first_of_col = k;
            }
        }
        if (first_of_col != 0) {
            // The column was split, the working part should be
            // the entire column. We erase the cells before first_of_col.
            // (Could be more efficient if we instead chop off end.)
            columns[col].erase(columns[col].begin(), columns[col].begin() + first_of_col);
        }
    }

    // Must tack on the new columns to complete the set.
    const int num_cols_all = num_cols + new_columns.size();
    columns.resize(num_cols_all);
    for (int col = num_cols; col < num_cols_all; ++col) {
        columns[col].swap(new_columns[col - num_cols]);
    }

}

} // namespace Opm
