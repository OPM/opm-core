/*
  Copyright 2015 Statoil ASA.

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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE EclipseRFTWriter
#include <boost/test/unit_test.hpp>

#include <opm/core/io/eclipse/EclipseWriteRFTHandler.hpp>
#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <ert/ecl/ecl_rft_file.h>
#include <ert/util/test_work_area.h>
#include <ert/util/util.h>

#include <vector>

namespace {

void verifyRFTFile(const std::string& rft_filename) {

    ecl_rft_file_type * new_rft_file = ecl_rft_file_alloc(rft_filename.c_str());
    std::shared_ptr<ecl_rft_file_type> rft_file;
    rft_file.reset(new_rft_file, ecl_rft_file_free);

    //Get RFT node for well/time OP_1/10 OKT 2008
    time_t recording_time = util_make_datetime(0, 0, 0, 10, 10, 2008);
    ecl_rft_node_type * ecl_rft_node = ecl_rft_file_get_well_time_rft(rft_file.get() , "OP_1" , recording_time);
    BOOST_CHECK(ecl_rft_node_is_RFT(ecl_rft_node));

    //Verify RFT data for completions (ijk) 9 9 1, 9 9 2 and 9 9 3 for OP_1
    const ecl_rft_cell_type * ecl_rft_cell1 = ecl_rft_node_lookup_ijk(ecl_rft_node, 8, 8, 0);
    const ecl_rft_cell_type * ecl_rft_cell2 = ecl_rft_node_lookup_ijk(ecl_rft_node, 8, 8, 1);
    const ecl_rft_cell_type * ecl_rft_cell3 = ecl_rft_node_lookup_ijk(ecl_rft_node, 8, 8, 2);

    BOOST_CHECK_CLOSE(ecl_rft_cell_get_pressure(ecl_rft_cell1), 210088*0.00001, 0.00001);
    BOOST_CHECK_CLOSE(ecl_rft_cell_get_pressure(ecl_rft_cell2), 210188*0.00001, 0.00001);
    BOOST_CHECK_CLOSE(ecl_rft_cell_get_pressure(ecl_rft_cell3), 210288*0.00001, 0.00001);

    BOOST_CHECK_EQUAL(ecl_rft_cell_get_sgas(ecl_rft_cell1), 0.0);
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_sgas(ecl_rft_cell2), 0.0);
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_sgas(ecl_rft_cell3), 0.0);

    BOOST_CHECK_EQUAL(ecl_rft_cell_get_swat(ecl_rft_cell1), 0.0);
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_swat(ecl_rft_cell2), 0.0);
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_swat(ecl_rft_cell3), 0.0);

    BOOST_CHECK_EQUAL(ecl_rft_cell_get_soil(ecl_rft_cell1), 1.0);
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_soil(ecl_rft_cell2), 1.0);
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_soil(ecl_rft_cell3), 1.0);

    BOOST_CHECK_EQUAL(ecl_rft_cell_get_depth(ecl_rft_cell1), (0.250 + (0.250/2)));
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_depth(ecl_rft_cell2), (2*0.250 + (0.250/2)));
    BOOST_CHECK_EQUAL(ecl_rft_cell_get_depth(ecl_rft_cell3), (3*0.250 + (0.250/2)));
}




Opm::DeckConstPtr createDeck(const std::string& input_str) {
    Opm::ParserPtr parser = std::make_shared<Opm::Parser>();
    Opm::DeckConstPtr deck = parser->parseString(input_str);
    return deck;
}


std::shared_ptr<Opm::WellState> createWellState(std::shared_ptr<Opm::BlackoilState> blackoilState)
{
    std::shared_ptr<Opm::WellState> wellState = std::make_shared<Opm::WellState>();
    wellState->init(0, *blackoilState);
    return wellState;
}



std::shared_ptr<Opm::BlackoilState> createBlackoilState(int timeStepIdx, std::shared_ptr<Opm::GridManager> ourFineGridManagerPtr)
{
    const UnstructuredGrid &ourFinerUnstructuredGrid = *ourFineGridManagerPtr->c_grid();

    std::shared_ptr<Opm::BlackoilState> blackoilState = std::make_shared<Opm::BlackoilState>();
    blackoilState->init(ourFinerUnstructuredGrid, 3);

    size_t numCells = ourFinerUnstructuredGrid.number_of_cells;

    auto &pressure = blackoilState->pressure();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx) {
        pressure[cellIdx] = timeStepIdx*1e5 + 1e4 + cellIdx;
    }
    return blackoilState;
}



std::shared_ptr<Opm::EclipseWriter> createEclipseWriter(std::shared_ptr<const Opm::Deck> deck,
                                                        std::shared_ptr<Opm::EclipseState> eclipseState,
                                                        std::shared_ptr<Opm::GridManager> ourFineGridManagerPtr,
                                                        const int * compressedToCartesianCellIdx)
{
    Opm::parameter::ParameterGroup params;
    params.insertParameter("deck_filename", "testcase.data");

    Opm::PhaseUsage phaseUsage = Opm::phaseUsageFromDeck(deck);

    const UnstructuredGrid &ourFinerUnstructuredGrid = *ourFineGridManagerPtr->c_grid();

    std::shared_ptr<Opm::EclipseWriter> eclipseWriter = std::make_shared<Opm::EclipseWriter>(params,
                                                                                             deck,
                                                                                             eclipseState,
                                                                                             phaseUsage,
                                                                                             ourFinerUnstructuredGrid.number_of_cells,
                                                                                             compressedToCartesianCellIdx);

    return eclipseWriter;
}

}

BOOST_AUTO_TEST_CASE(test_EclipseWriterRFTHandler)
{
    const std::string& deckString =
                                    "RUNSPEC\n"
                                    "OIL\n"
                                    "GAS\n"
                                    "WATER\n"
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
                                    "\n"
                                     "START             -- 0 \n"
                                    "1 NOV 1979 / \n"
                                    "SCHEDULE\n"
                                    "DATES             -- 1\n"
                                    " 1 DES 1979/ \n"
                                    "/\n"
                                    "WELSPECS\n"
                                    "    'OP_1'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
                                    "    'OP_2'       'OP'   4   4 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
                                    "/\n"
                                    "COMPDAT\n"
                                    " 'OP_1'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
                                    " 'OP_1'  9  9   2   2 'OPEN' 1*   46.825   0.311  4332.346 1*  1*  'X'  22.123 / \n"
                                    " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
                                    " 'OP_2'  4  4   4  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
                                    "/\n"
                                    "DATES             -- 2\n"
                                    " 10  OKT 2008 / \n"
                                    "/\n"
                                    "WRFT \n"
                                    "/ \n"
                                    "WELOPEN\n"
                                    " 'OP_1' OPEN / \n"
                                    " 'OP_2' OPEN / \n"
                                    "/\n"
                                    "DATES             -- 3\n"
                                    " 10  NOV 2008 / \n"
                                    "/\n";



    test_work_area_type * new_ptr = test_work_area_alloc("test_EclipseWriterRFTHandler");
    std::shared_ptr<test_work_area_type> test_area;
    test_area.reset(new_ptr, test_work_area_free);

    std::shared_ptr<const Opm::Deck>   deck         = createDeck(deckString);
    std::shared_ptr<Opm::EclipseState> eclipseState = std::make_shared<Opm::EclipseState>(deck);

    std::shared_ptr<Opm::SimulatorTimer> simulatorTimer = std::make_shared<Opm::SimulatorTimer>();
    simulatorTimer->init(eclipseState->getSchedule()->getTimeMap());

    std::shared_ptr<Opm::GridManager>  ourFineGridManagerPtr = std::make_shared<Opm::GridManager>(eclipseState->getEclipseGrid());
    const UnstructuredGrid &ourFinerUnstructuredGrid = *ourFineGridManagerPtr->c_grid();
    const int* compressedToCartesianCellIdx = Opm::UgGridHelpers::globalCell(ourFinerUnstructuredGrid);

    std::shared_ptr<Opm::EclipseWriter> eclipseWriter = createEclipseWriter(deck,
                                                                            eclipseState,
                                                                            ourFineGridManagerPtr,
                                                                            compressedToCartesianCellIdx);
    eclipseWriter->writeInit(*simulatorTimer);


    for (; simulatorTimer->currentStepNum() < simulatorTimer->numSteps(); ++ (*simulatorTimer)) {
        std::shared_ptr<Opm::BlackoilState> blackoilState2 = createBlackoilState(simulatorTimer->currentStepNum(),ourFineGridManagerPtr);
        std::shared_ptr<Opm::WellState> wellState = createWellState(blackoilState2);
        eclipseWriter->writeTimeStep(*simulatorTimer, *blackoilState2, *wellState);
    }

    std::string cwd(test_work_area_get_cwd(test_area.get()));
    std::string rft_filename = cwd + "/TESTCASE.RFT";
    verifyRFTFile(rft_filename);

}



