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


#include "config.h"

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/eclipse/EclipseGridInspector.hpp>
#include <opm/core/fluid/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/grid.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    // Parameters.
    Opm::parameter::ParameterGroup param(argc, argv);

    // Parser.
    std::string ecl_file = param.get<std::string>("deck_filename");
    std::string input_file = param.get<std::string>("input_filename");
    std::string matrix_output = param.get<std::string>("matrix_outout");
    std::string rs_output = param.get<std::string>("rs_output");
    std::string b_output = param.get<std::string>("b_output");
    Opm::EclipseGridParser deck(ecl_file);
    UnstructuredGrid grid;
    grid.number_of_cells = 1;
    grid.global_cell = NULL;
    Opm::BlackoilPropertiesFromDeck props(deck, grid);
    Opm::BlackoilPvtProperties pvt;
    pvt.init(deck);
    
  
    const int n = 1;
    std::fstream inos(input_file.c_str());
    int np=props.numPhases();
    while(inos.good()){
        double p[n];
        double z[np*n];
        int cells[n] = { 0 };
        inos >> p[0];      
        for(int i=0; i < np; ++i){
            inos >> z[i];
        }
        double A[np*np*n];
        props.matrix(n, p, z, cells, A, 0);
        std::fstream aos(matrix_output.c_str());
        std::copy(A, A + np*np, std::ostream_iterator<double>(aos, " "));
        std::cout << std::endl;
        double b[np];
        double dbdp[np];            
        pvt.dBdp(n, p, z, b, dbdp);
        std::fstream bos(b_output.c_str());
        std::copy(b, b + np, std::ostream_iterator<double>(bos, " "));
        std::copy(b, dbdp + np, std::ostream_iterator<double>(bos, " "));            
        std::cout << std::endl;
        double rs[np];
        double drs[np];
        //pvt.R(n, p, z, rs);
        pvt.dRdp(n, p, z, rs,drs);            
        std::fstream rsos(rs_output.c_str());
        std::copy(b, rs + np, std::ostream_iterator<double>(rsos, " "));
        std::copy(b, drs + np, std::ostream_iterator<double>(rsos, " "));
        std::cout << std::endl;
    }
}
