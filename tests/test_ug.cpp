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
        const auto& dimens = deck->getKeyword("DIMENS");
        const auto& coord = deck->getKeyword("COORD");
        const auto& zcorn = deck->getKeyword("ZCORN");
        const auto& actnum = deck->getKeyword("ACTNUM");
        
        g.dims[0] = dimens.getRecord(0).getItem("NX").get< int >(0);
        g.dims[1] = dimens.getRecord(0).getItem("NY").get< int >(0);
        g.dims[2] = dimens.getRecord(0).getItem("NZ").get< int >(0);

        g.coord  = coord.getSIDoubleData().data();
        g.zcorn  = zcorn.getSIDoubleData().data();
        g.actnum = actnum.getIntData().data();
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


