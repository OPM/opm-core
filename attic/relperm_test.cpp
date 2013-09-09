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
#include <iterator>
#include <string>
#include <vector>
#include <fstream>

int main(int argc, char** argv)
{
    // Parameters.
    Opm::parameter::ParameterGroup param(argc, argv);

    // Parser.
    std::string ecl_file = param.get<std::string>("deck_filename");
    std::string input_file = param.get<std::string>("input_filename");
    std::string relperm_output = param.get<std::string>("relperm_output");
    //std::string relperm_output = param.get<std::string>("relperm_outout");
    Opm::EclipseGridParser deck(ecl_file);
    UnstructuredGrid grid;
    grid.number_of_cells = 1;
    grid.global_cell = NULL;
    grid.dimensions = 3;
    grid.cartdims[0] = 1;
    grid.cartdims[1] = 1;
    grid.cartdims[2] = 1;
    Opm::BlackoilPropertiesFromDeck props(deck, grid, param);

    std::fstream inos(input_file.c_str());//, std::fstream::in);
    if(!inos.good()){
        std::cout << "Could not open :" << input_file << std::endl;
        exit(3);
    }
    std::fstream kros(relperm_output.c_str(), std::fstream::out | std::fstream::trunc);
    if(!kros.good()){
        std::cout << "Could not open :" << input_file << std::endl;
        exit(3);
    }
    const int np = props.numPhases();
    const int max_np = 3;
    if (np > max_np) {
        OPM_THROW(std::runtime_error, "Max #phases is 3.");
    }
    while((inos.good()) && (!inos.eof())){
      double s[max_np];
      for(int i=0; i < np; ++i){
        inos >> s[i];
      }
      if(inos.good()){
          double kr[max_np];
          double dkr[max_np*max_np];
          int cell[1];
          cell[0]=1;
          props.relperm(1,s, cell, kr, dkr);
          std::copy(s, s + np, std::ostream_iterator<double>(kros, " "));
          kros << " ";
          std::copy(kr, kr + np, std::ostream_iterator<double>(kros, " "));
          kros << " ";
          std::copy(dkr, dkr + np*np, std::ostream_iterator<double>(kros, " "));
          kros << "\n";
      }
    }
    if (param.anyUnused()) {
        std::cout << "--------------------   Unused parameters:   --------------------\n";
        param.displayUsage();
        std::cout << "----------------------------------------------------------------" << std::endl;
    }
}
