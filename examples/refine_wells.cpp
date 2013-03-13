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

#include <opm/core/io/eclipse/EclipseGridParser.hpp>

// Double I and J coordinates of wells and completions.
// Do not change any well productivity indices or any other data.
int main(int argc, char** argv)
{
    using namespace Opm;
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " deck\n";
        return 1;
    }
    EclipseGridParser deck(argv[1], false);

    WELSPECS ws = deck.getWELSPECS();
    const int nw = ws.welspecs.size();
    for (int w = 0; w < nw; ++w) {
        ws.welspecs[w].I_ *= 2;
        ws.welspecs[w].J_ *= 2;
    }

    COMPDAT cd = deck.getCOMPDAT();
    const int nc = cd.compdat.size();
    for (int c = 0; c < nc; ++c) {
        cd.compdat[c].grid_ind_[0] *= 2;
        cd.compdat[c].grid_ind_[1] *= 2;
    }

    ws.write(std::cout);
    std::cout << '\n';
    cd.write(std::cout);
}
