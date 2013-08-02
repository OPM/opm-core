/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/grid.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    using namespace std;
    // Parameters.
    Opm::parameter::ParameterGroup param(argc, argv);

    // Parser.
    std::string ecl_file = param.get<std::string>("deck_filename");
    std::string input_file = param.get<std::string>("input_filename");
    std::string matrix_output = param.get<std::string>("matrix_output");
    std::string rs_output = param.get<std::string>("rs_output");
    std::string b_output = param.get<std::string>("b_output");
    std::string mu_output = param.get<std::string>("mu_output");
    Opm::EclipseGridParser deck(ecl_file);
    UnstructuredGrid grid;
    grid.number_of_cells = 1;
    grid.global_cell = NULL;
    grid.dimensions = 3;
    grid.cartdims[0] = 1;
    grid.cartdims[1] = 1;
    grid.cartdims[2] = 1;
    Opm::BlackoilPropertiesFromDeck props(deck, grid, param);
    Opm::BlackoilPvtProperties pvt;
    int samples = param.getDefault("dead_tab_size", 1025);
    pvt.init(deck, samples);


    const int n = 1;
    std::fstream inos(input_file.c_str());
    if(!inos.good()){
        std::cout << "Could not open :" << input_file << std::endl;
        exit(3);
    }
    std::fstream aos(matrix_output.c_str(), std::fstream::out | std::fstream::trunc);
    aos << setiosflags(ios::scientific) << setprecision(12);
    if(!(aos.good())){
        std::cout << "Could not open"<< matrix_output << std::endl;
        exit(3);
    }
    std::fstream bos(b_output.c_str(), std::fstream::out | std::fstream::trunc);
    bos << setiosflags(ios::scientific) << setprecision(12);
    if(!(bos.good())){
        std::cout << "Could not open"<< b_output << std::endl;
        exit(3);
    }
    std::fstream rsos(rs_output.c_str(), std::fstream::out | std::fstream::trunc);
    rsos << setiosflags(ios::scientific) << setprecision(12);
    if(!(rsos.good())){
        std::cout << "Could not open"<< rs_output << std::endl;
        exit(3);
    }
    std::fstream muos(mu_output.c_str(), std::fstream::out | std::fstream::trunc);
    if(!(muos.good())){
        std::cout << "Could not open"<< rs_output << std::endl;
        exit(3);
    }



    const int np = props.numPhases();
    const int max_np = 3;
    if (np > max_np) {
        THROW("Max #phases is 3.");
    }
    while((inos.good()) && (!inos.eof())){
        double p[n];
        double z[max_np*n];
        int cells[n] = { 0 };
        inos >> p[0];
        for(int i=0; i < np; ++i){
            inos >> z[i];
        }

        if(inos.good()){
            double A[max_np*max_np*n];
            double dA[max_np*max_np*n];
            props.matrix(n, p, z, cells, A, dA);
            std::copy(A, A + np*np*n, std::ostream_iterator<double>(aos, " "));
            std::copy(dA, dA + np*np*n, std::ostream_iterator<double>(aos, " "));
            aos << std::endl;
            double mu[max_np];
            //double dmu[max_np];//not implemented
            props.viscosity(n, p, z, cells, mu, 0);
            std::copy(mu, mu + np*n, std::ostream_iterator<double>(muos, " "));
            //std::copy(dmu, dmu + np*n, std::ostream_iterator<double>(muos, " "));
            aos << std::endl;

            double b[max_np];
            double dbdp[max_np];
            pvt.dBdp(n, p, z, b, dbdp);
            std::copy(b, b + np*n, std::ostream_iterator<double>(bos, " "));
            std::copy(dbdp, dbdp + np*n, std::ostream_iterator<double>(bos, " "));
            bos << std::endl;

            double rs[max_np];
            double drs[max_np];
            //pvt.R(n, p, z, rs);
            pvt.dRdp(n, p, z, rs,drs);
            std::copy(rs, rs + np*n, std::ostream_iterator<double>(rsos, " "));
            std::copy(drs, drs + np*n, std::ostream_iterator<double>(rsos, " "));
            rsos << std::endl;
        }
    }
    if (param.anyUnused()) {
        std::cout << "--------------------   Unused parameters:   --------------------\n";
        param.displayUsage();
        std::cout << "----------------------------------------------------------------" << std::endl;
    }
}
