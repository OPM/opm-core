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

#include "config.h"
/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE UnitsTest
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

/* --- our own headers --- */
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include <boost/bind.hpp>

#include <opm/core/utility/Units.hpp>

using namespace Opm::prefix;
using namespace Opm::unit;

BOOST_AUTO_TEST_SUITE ()

BOOST_AUTO_TEST_CASE (units)
{
    BOOST_REQUIRE_EQUAL (meter, 1);
    BOOST_REQUIRE_EQUAL (kilogram, 1);
    BOOST_REQUIRE_EQUAL (second, 1);

    BOOST_REQUIRE_CLOSE (milli*darcy, 9.86923667e-16, 0.01);
    BOOST_REQUIRE_CLOSE (mega*darcy, 9.86923e-7, 0.01);
    BOOST_REQUIRE_CLOSE (convert::to(mega*darcy, milli*darcy), 1e9, 0.01);

    BOOST_REQUIRE_CLOSE (convert::to(convert::from(1.0, barsa), psia), 14.5038, 0.01);
    BOOST_REQUIRE_CLOSE (convert::to(1*atm, barsa), 1.01325, 0.01);

    std::vector<double> flux(10, 10000*cubic(meter)/year);
	for (int i = 0; i < 10; ++i) {
        BOOST_REQUIRE_CLOSE (flux[i], 3.17098e-4, 0.01); 
	}

    std::transform(flux.begin(), flux.end(), flux.begin(),
                   boost::bind(convert::to, _1, cubic(meter)/year));
	for (int i = 0; i < 10; ++i) {
        BOOST_REQUIRE_CLOSE (flux[i], 1e4, 0.01);
	}
}

BOOST_AUTO_TEST_SUITE_END()
