#include <opm/core/grid.h>
#include <vector>

namespace Opm {

    /// Extract each column of the grid.
    ///  \note Assumes the pillars of the grid are all vertically aligned.
    ///  \param grid The grid from which to extract the columns.
    ///  \param columns will for each i + cartgrim[0]*j  contain the k values contained in the column
    ///         centered at (i, j).
    void extractColumn( const UnstructuredGrid& grid, std::vector<std::vector<int> >& columns )
    {
	const int* dims = grid.cartdims;
	columns.resize(dims[0]*dims[1]);
	for (int i = 0; i < dims[0]; ++i) {
	    for (int j = 0; j < dims[1]; ++j) {
		int plane_index = i + j*dims[0];
		for (int k = 0; k < dims[2]; ++k) {
		    columns[plane_index].push_back(plane_index + k*dims[0]*dims[1]);
		}
	    }
	}
    }

} // namespace Opm
