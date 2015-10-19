/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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


#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE MinpvProcessorTest
#include <boost/test/unit_test.hpp>

#include <opm/core/grid/MinpvProcessor.hpp>

BOOST_AUTO_TEST_CASE(Processing)
{
    std::vector<double> zcorn = { 0, 0, 0, 0,
                                  1, 1, 1, 1,
                                  1, 1, 1, 1, 
                                  3, 3, 3, 3,
                                  3, 3, 3, 3,
                                  6, 6, 6, 6 };
    std::vector<double> zcorn2after = { 0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        0, 0, 0, 0, 
                                        3, 3, 3, 3,
                                        3, 3, 3, 3,
                                        6, 6, 6, 6 };
    std::vector<double> zcorn3after = { 0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        0, 0, 0, 0, 
                                        0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        6, 6, 6, 6  };
    std::vector<double> zcorn4after = { 0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        1, 1, 1, 1,
                                        3, 3, 3, 3,
                                        3, 3, 3, 3,
                                        6, 6, 6, 6 };
    std::vector<double> zcorn5after = { 0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        1, 1, 1, 1,
                                        1, 1, 1, 1,
                                        3, 3, 3, 3,
                                        6, 6, 6, 6 };

    std::vector<double> pv = { 1, 2, 3 };
    std::vector<int> actnum = { 1, 1, 1 };

    Opm::MinpvProcessor mp1(1, 1, 3);
    auto z1 = zcorn;
    bool fill_removed_cells = true;
    mp1.process(pv, 0.5, actnum, fill_removed_cells, z1.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcorn.begin(), zcorn.end());

    Opm::MinpvProcessor mp2(1, 1, 3);
    auto z2 = zcorn;
    mp2.process(pv, 1.5, actnum, fill_removed_cells, z2.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z2.begin(), z2.end(), zcorn2after.begin(), zcorn2after.end());

    Opm::MinpvProcessor mp3(1, 1, 3);
    auto z3 = zcorn;
    mp3.process(pv, 2.5, actnum, fill_removed_cells, z3.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z3.begin(), z3.end(), zcorn3after.begin(), zcorn3after.end());

    Opm::MinpvProcessor mp4(1, 1, 3);
    auto z4 = zcorn;
    mp4.process(pv, 1.5, actnum, !fill_removed_cells, z4.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z4.begin(), z4.end(), zcorn4after.begin(), zcorn4after.end());

    Opm::MinpvProcessor mp5(1, 1, 3);
    auto z5 = zcorn;
    mp5.process(pv, 2.5, actnum, !fill_removed_cells, z5.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z5.begin(), z5.end(), zcorn5after.begin(), zcorn5after.end());
}
