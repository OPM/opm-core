/* Copyright 2014 Statoil ASA
 * This file is licensed under GPL3, see http://www.opm-project.org/
*/

#include <config.h>

/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE TEST_UG
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

/* --- our own headers --- */
#include <algorithm>
#include <vector>
#include <opm/core/grid.h>
#include <opm/core/grid/cornerpoint_grid.h>  /* compute_geometry */
#include <opm/core/grid/GridManager.hpp>  /* compute_geometry */

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

using namespace std;



BOOST_AUTO_TEST_CASE(Equal) {
    const std::string filename1 = "CORNERPOINT_ACTNUM.DATA";
    const char *deck2Data =
        "RUNSPEC\n"
        "\n"
        "DIMENS\n"
        " 10 10 10 /\n"
        "GRID\n"
        "DXV\n"
        "10*0.25 /\n"
        "DYV\n"
        "10*0.25 /\n"
        "DZV\n"
        "10*0.25 /\n"
        "TOPS\n"
        "100*0.25 /\n"
        "EDIT\n"
        "\n";
    
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck1 = parser->parseFile( filename1 );
    Opm::DeckConstPtr deck2 = parser->parseString( deck2Data );
    
    Opm::GridManager grid1(deck1);
    Opm::GridManager grid2(deck2);
    
    const UnstructuredGrid* cgrid1 = grid1.c_grid();
    const UnstructuredGrid* cgrid2 = grid2.c_grid();
    

    BOOST_CHECK( grid_equal( cgrid1 , cgrid1 ));
    BOOST_CHECK( grid_equal( cgrid2 , cgrid2 ));
    BOOST_CHECK( !grid_equal( cgrid1 , cgrid2 ));
}

