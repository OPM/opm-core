//===========================================================================
//
// File: unit_test.cpp
//
// Created: Fri Jul 17 11:02:33 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include <boost/bind.hpp>

#include <opm/core/utility/Units.hpp>

using namespace Dune::prefix;
using namespace Dune::unit;

int main()
{
    std::cout << "m  = " << meter << '\n';
    std::cout << "kg = " << kilogram << '\n';
    std::cout << "s  = " << second << '\n';

    std::cout << "mD = " << milli*darcy << '\n';
    std::cout << "MD = " << mega*darcy << '\n';

    std::cout << "MD = " << convert::to(mega*darcy, milli*darcy) << " mD\n";

    std::cout << "1 bar = " << convert::to(convert::from(1.0, barsa), psia) << " psi\n";

    std::cout << "1 atm = " << convert::to(1*atm, barsa) << " bar\n";

    std::vector<double> flux(10, 10000*cubic(meter)/year);

    std::cout << "flux = [";
    std::copy(flux.begin(), flux.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\b] (m^3/s)\n";

    std::transform(flux.begin(), flux.end(), flux.begin(),
                   boost::bind(convert::to, _1, cubic(meter)/year));
    std::cout << "     = [";
    std::copy(flux.begin(), flux.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\b] (m^3/year)\n";

    return 0;
}
