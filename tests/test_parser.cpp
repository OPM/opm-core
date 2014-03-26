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

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE OPM-ParserTest
#include <boost/test/unit_test.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckDoubleItem.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <memory>

BOOST_AUTO_TEST_CASE(CreateParser)
{
    const std::string filename1 = "testBlackoilState1.DATA";
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile( filename1 );

    BOOST_CHECK_EQUAL( 6U , deck->size() );
    Opm::DeckItemConstPtr actnum = deck->getKeyword("ACTNUM")->getRecord(0)->getItem(0);
    const std::vector<int>& actnum_data = actnum->getIntData();

    BOOST_CHECK_EQUAL( 1000U , actnum->size() );
    BOOST_CHECK_EQUAL( 1, actnum_data[0] );
    BOOST_CHECK_EQUAL( 2, actnum_data[400] );
    BOOST_CHECK_EQUAL( 3, actnum_data[999] );
}
