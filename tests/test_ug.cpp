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
#include <opm/core/utility/opm_memcmp_double.h>
#include <opm/core/grid/cornerpoint_grid.h>  /* compute_geometry */
#include <opm/core/grid/GridManager.hpp>  /* compute_geometry */
#include <opm/core/grid/cpgpreprocess/preprocess.h>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseMode.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckItem.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>

using namespace std;



BOOST_AUTO_TEST_CASE(Equal) {
    Opm::ParseMode parseMode;
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
    Opm::DeckConstPtr deck1 = parser->parseFile( filename1 , parseMode);
    Opm::DeckConstPtr deck2 = parser->parseString( deck2Data , parseMode);
    
    BOOST_CHECK( deck1->hasKeyword("ZCORN") );
    BOOST_CHECK( deck1->hasKeyword("COORD") );
    
    Opm::GridManager grid1(deck1);
    Opm::GridManager grid2(deck2);
    
    const UnstructuredGrid* cgrid1 = grid1.c_grid();
    const UnstructuredGrid* cgrid2 = grid2.c_grid();
    


    BOOST_CHECK( grid_equal( cgrid1 , cgrid1 ));
    BOOST_CHECK( grid_equal( cgrid2 , cgrid2 ));
    BOOST_CHECK( !grid_equal( cgrid1 , cgrid2 ));
}



BOOST_AUTO_TEST_CASE(EqualEclipseGrid) {
    const std::string filename = "CORNERPOINT_ACTNUM.DATA";
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::ParseMode parseMode;
    Opm::DeckConstPtr deck = parser->parseFile( filename , parseMode);

    std::shared_ptr<const Opm::EclipseGrid> grid(new Opm::EclipseGrid(deck));

    Opm::GridManager gridM(grid);
    const UnstructuredGrid* cgrid1 = gridM.c_grid();
    struct UnstructuredGrid * cgrid2;
    {
        struct grdecl g;
        Opm::DeckKeywordConstPtr dimens = deck->getKeyword("DIMENS");
        Opm::DeckKeywordConstPtr coord = deck->getKeyword("COORD");
        Opm::DeckKeywordConstPtr zcorn = deck->getKeyword("ZCORN");
        Opm::DeckKeywordConstPtr actnum = deck->getKeyword("ACTNUM");
        
        g.dims[0] = dimens->getRecord(0)->getItem("NX")->getInt(0);
        g.dims[1] = dimens->getRecord(0)->getItem("NY")->getInt(0);
        g.dims[2] = dimens->getRecord(0)->getItem("NZ")->getInt(0);

        g.coord  = coord->getSIDoubleData().data();
        g.zcorn  = zcorn->getSIDoubleData().data();
        g.actnum = actnum->getIntData().data();
        g.mapaxes = NULL;
    
        
        cgrid2 = create_grid_cornerpoint(&g , 0.0);
        if (!cgrid2) 
            throw std::runtime_error("Failed to construct grid.");
    }
    
    
    BOOST_CHECK( grid_equal( cgrid1 , cgrid2 ));
    destroy_grid( cgrid2 );
}


BOOST_AUTO_TEST_CASE(TOPS_Fully_Specified) {
    const char *deck1Data =
        "RUNSPEC\n"
        "\n"
        "DIMENS\n"
        " 10 10 3 /\n"
        "GRID\n"
        "DX\n"
        "300*1000 /\n"
        "DY\n"
        "300*1000 /\n"
        "DZ\n"
        "100*20 100*30  100*50 /\n"
        "TOPS\n"
        "100*8325 /\n"
        "EDIT\n"
        "\n";


    const char *deck2Data =
        "RUNSPEC\n"
        "\n"
        "DIMENS\n"
        " 10 10 3 /\n"
        "GRID\n"
        "DX\n"
        "300*1000 /\n"
        "DY\n"
        "300*1000 /\n"
        "DZ\n"
        "100*20 100*30  100*50 /\n"
        "TOPS\n"
        "100*8325 100*8345  100*8375/\n"
        "EDIT\n"
        "\n";

    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::ParseMode parseMode;
    Opm::DeckConstPtr deck1 = parser->parseString( deck1Data , parseMode);
    Opm::DeckConstPtr deck2 = parser->parseString( deck2Data , parseMode);

    std::shared_ptr<const Opm::EclipseGrid> grid1(new Opm::EclipseGrid(deck1));
    std::shared_ptr<const Opm::EclipseGrid> grid2(new Opm::EclipseGrid(deck2));

    Opm::GridManager gridM1(grid1);
    Opm::GridManager gridM2(grid2);

    const UnstructuredGrid* cgrid1 = gridM1.c_grid();
    const UnstructuredGrid* cgrid2 = gridM2.c_grid();

    BOOST_CHECK( grid_equal( cgrid1 , cgrid2 ));
}



BOOST_AUTO_TEST_CASE(compare_double) {
    const double abs_epsilon = 1e-8;
    const double rel_epsilon = 1e-5;

    double v1,v2;
    /* Should be equal: */
    {
        v1 = 0.0;
        v2 = 0.0;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 0 );

        v1 = 1e-12;
        v2 = v1 + 0.5*abs_epsilon;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 0 );

        v1 = 7.0;
        v2 = 7.0;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 0 );

        v1 = -7.0;
        v2 = -7.0;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 0 );

        v1 = 0;
        v2 = 0.5 * abs_epsilon;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 0 );


        v1 = 1e7;
        v2 = 1e7 + 2*abs_epsilon;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 0 );

        v1 = 1e7;
        v2 = 1e7*(1 + 2*rel_epsilon);
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 0 );
    }

    /* Should be different: */
    {
        v1 = 0;
        v2 = 1.5 * abs_epsilon;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 1 );

        v1 = 1e-8;
        v2 = v1 + 1.5*abs_epsilon;
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 1 );

        v1 = 1;
        v2 = v1*(1 + 2*rel_epsilon + abs_epsilon);
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 1 );

        v1 = 10;
        v2 = v1*(1 + 2*rel_epsilon + abs_epsilon);
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 1 );

        v1 = 1e7;
        v2 = 1e7*(1 + 2*rel_epsilon + abs_epsilon);
        BOOST_CHECK_EQUAL( opm_memcmp_double( &v1 , &v2 , 1 ) , 1 );
    }
}
