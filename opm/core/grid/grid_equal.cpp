#include <string.h>  // C string.h to get memcmp()

#include <opm/common/util/numeric/cmp.hpp>

#include <opm/core/grid.h>

/*
   The grid_equal() function is separated out into a separate file to
   ensure that it is compiled with a C++ compiler, so that the
   namespace features used in the opm/common/util/numeric/cmp.cpp
   implementation compiles.
*/


bool
grid_equal(const struct UnstructuredGrid * grid1 , const struct UnstructuredGrid * grid2) {
    if ((grid1->dimensions      == grid2->dimensions)      &&
        (grid1->number_of_cells == grid2->number_of_cells) &&
        (grid1->number_of_faces == grid2->number_of_faces) &&
        (grid1->number_of_nodes == grid2->number_of_nodes)) {

        // Exact integer comparisons
        {
            if (memcmp(grid1->face_nodepos , grid2->face_nodepos , (grid1->number_of_faces + 1) * sizeof * grid1->face_nodepos) != 0)
            return false;

            if (memcmp(grid1->face_nodes , grid2->face_nodes , grid1->face_nodepos[grid1->number_of_faces] * sizeof * grid1->face_nodes) != 0)
                return false;

            if (memcmp(grid1->face_cells , grid2->face_cells , 2 * grid1->number_of_faces * sizeof * grid1->face_cells) != 0)
                return false;

            if (memcmp(grid1->cell_faces , grid2->cell_faces ,  grid1->cell_facepos[grid1->number_of_cells] * sizeof * grid1->cell_faces) != 0)
                return false;

            if (memcmp(grid1->cell_facepos , grid2->cell_facepos , (grid1->number_of_cells + 1) * sizeof * grid1->cell_facepos) != 0)
                return false;

            if (grid1->global_cell && grid2->global_cell) {
                if (memcmp(grid1->global_cell , grid2->global_cell , grid1->number_of_cells * sizeof * grid1->global_cell) != 0)
                    return false;
            } else if (grid1->global_cell != grid2->global_cell)
                return false;

            if (grid1->cell_facetag && grid2->cell_facetag) {
                if (memcmp(grid1->cell_facetag , grid2->cell_facetag , grid1->cell_facepos[grid1->number_of_cells] * sizeof * grid1->cell_facetag) != 0)
                    return false;
            } else if (grid1->cell_facetag != grid2->cell_facetag)
                return false;
        }


        // Floating point comparisons.
        {
	    if (!Opm::cmp::array_equal<double>( grid1->node_coordinates , grid2->node_coordinates , static_cast<size_t>(grid1->dimensions * grid1->number_of_nodes)))
                return false;

  	    if (!Opm::cmp::array_equal<double>( grid1->face_centroids , grid2->face_centroids , static_cast<size_t>(grid1->dimensions * grid1->number_of_faces)))
                return false;

            if (!Opm::cmp::array_equal<double>( grid1->face_areas , grid2->face_areas , static_cast<size_t>(grid1->number_of_faces)))
                return false;

            if (!Opm::cmp::array_equal<double>( grid1->face_normals , grid2->face_normals , static_cast<size_t>(grid1->dimensions * grid1->number_of_faces)))
                return false;

            if (!Opm::cmp::array_equal<double>( grid1->cell_centroids , grid2->cell_centroids , static_cast<size_t>(grid1->dimensions * grid1->number_of_cells)))
                return false;

            if (!Opm::cmp::array_equal<double>( grid1->cell_volumes , grid2->cell_volumes , static_cast<size_t>(grid1->number_of_cells)))
                return false;
        }
        return true;
    } else
        return false;
}
