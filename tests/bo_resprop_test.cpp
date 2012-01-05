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
#include <opm/core/fluid/BlackoilPropertiesFromDeck.hpp>

#include <iterator>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    // Parameters.
    Dune::parameter::ParameterGroup param(argc, argv);

    // Parser.
    std::string ecl_file = param.get<std::string>("filename");
    Dune::EclipseGridParser deck(ecl_file);

    Opm::BlackoilPropertiesFromDeck props(deck, ::std::vector<int>(1));

    const int n = 1;
    double p[n] = { 150e5 };
    double z[3*n] = { 1.0, 1.0, 1.0 };
    int cells[n] = { 0 };
    double A[9*n];
    props.matrix(n, p, z, cells, A, 0);
    std::copy(A, A + 9, std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
}
