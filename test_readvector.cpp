#include <iostream>
#include <string>
#include <vector>

#include "grid.h"
#include "preprocess.h"
#include "cgridinterface.h"

#include "readvector.hpp"

struct CornerpointGrid
read_grid(const std::string& dir)
{
    std::string fn;
    fn = dir + '/' + "zcorn.txt";

    std::vector<double> zcorn;
    read_vector_from_file(fn, zcorn);

    fn = dir + '/' + "coord.txt";
    ::std::vector<double> coord;
    read_vector_from_file(fn, coord);

    fn = dir + '/' + "actnum.txt";
    std::vector<int> actnum;
    read_vector_from_file(fn, actnum);

    fn = dir + '/' + "dimens.txt";
    ::std::vector<double> dimens;
    read_vector_from_file(fn, dimens);

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

    return G;
}

int main()
{
    struct CornerpointGrid G;

    G = read_grid(std::string("example"));

    free_cornerpoint_grid(&G);
    
    return 0;
}
