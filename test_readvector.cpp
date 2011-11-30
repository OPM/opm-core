#include <iostream>
#include <vector>
#include <cpgpreprocess/preprocess.h>
#include <cpgpreprocess/cgridinterface.h>

#include "readvector.hpp"

int main()
{
    std::vector<double> zcorn;
    read_vector_from_file("example/zcorn.txt", zcorn);
    std::vector<int> actnum;
    read_vector_from_file("example/actnum.txt", actnum);

    ::std::vector<double> coord;
    read_vector_from_file("example/coord.txt", coord);

    ::std::vector<double> dimens;
    read_vector_from_file("example/dimens.txt", dimens);

    struct grdecl grdecl;
    grdecl.zcorn = &zcorn[0];
    grdecl.coord = &coord[0];
    grdecl.actnum = &actnum[0];

    grdecl.dims[0] = dimens[0];
    grdecl.dims[1] = dimens[1];
    grdecl.dims[2] = dimens[2];

    struct CornerpointGrid G;

    preprocess(&grdecl, 0.0, &G);

    compute_geometry(&G);

    struct UnstructuredGrid *g = &G.grid;

    double vol = 0.0;
    for (int c = 0; c < g->number_of_cells; c++) {
        vol += g->cell_volumes[c];
    }
    std::cout << "Sum volumes = " << vol << '\n';

    for (int c = 0, i = 0; c < g->number_of_cells; c++) {
        for (; i < g->cell_facepos[c + 1]; i++) {
            std::cout << "(c,i) = (" << c << "," << G.cface_tag[i] << ")\n";
        }
    }

    free_cornerpoint_grid(&G);

    return 0;
}
