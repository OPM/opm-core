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
#include <opm/core/utility/writeVtkData.hpp>
#include <fstream>
#include <array>

int main()
{
    std::array<int, 3> dims = {{ 2, 2, 2 }};
    std::array<double, 3> cell_size = {{ 1.0, 2.0, 3.0 }};
    Opm::DataMap m;
    double foov[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
    std::vector<double> foo(foov, foov + sizeof foov/sizeof foov[0]);
    m["foo"] = &foo;
    double barv[16] = { 1, 10, 2, 20, 3, 30, 4, 40, 5, 50, 6, 60, 7, 70, 8, 80 };
    std::vector<double> bar(barv, barv + sizeof barv/sizeof barv[0]);
    m["bar"] = &bar;
    std::ofstream os("testoutput.vtk");
    Opm::writeVtkData(dims, cell_size, m, os);
}
