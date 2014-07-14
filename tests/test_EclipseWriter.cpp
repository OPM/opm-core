/*
  Copyright 2014 Andreas Lauser

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

#define BOOST_TEST_MODULE EclipseWriter
#include <boost/test/unit_test.hpp>

#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

// ERT stuff
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_endian_flip.h>
#include <ert/ecl/fortio.h>

#include <memory>

std::shared_ptr<Opm::EclipseWriter> eclWriter;
std::shared_ptr<Opm::SimulatorTimer> simTimer;
std::shared_ptr<const Opm::Deck> deck;
std::shared_ptr<Opm::EclipseGrid> eclGrid;
std::shared_ptr<Opm::GridManager> ourFineGridManagerPtr;
std::shared_ptr<Opm::BlackoilState> blackoilState;
std::shared_ptr<Opm::WellState> wellState;

void createEclipseWriter(const char *deckString)
{
    Opm::ParserConstPtr parser(new Opm::Parser());
    deck = parser->parseString(deckString);

    Opm::parameter::ParameterGroup params;
    params.insertParameter("deck_filename", "foo.data");

    auto runspecSection = std::make_shared<Opm::RUNSPECSection>(deck);
    auto gridSection = std::make_shared<Opm::GRIDSection>(deck);
    eclGrid.reset(new Opm::EclipseGrid(runspecSection, gridSection));

    BOOST_CHECK(eclGrid->getNX() == 3);
    BOOST_CHECK(eclGrid->getNY() == 3);
    BOOST_CHECK(eclGrid->getNZ() == 3);
    BOOST_CHECK(eclGrid->getCartesianSize() == 3*3*3);

    std::array<int, 3> cartSize;
    cartSize[0] = eclGrid->getNX();
    cartSize[1] = eclGrid->getNY();
    cartSize[2] = eclGrid->getNZ();

    simTimer.reset(new Opm::SimulatorTimer());
    Opm::TimeMapConstPtr timeMap(new Opm::TimeMap(deck));
    simTimer->init(timeMap);

    eclWriter.reset(new Opm::EclipseWriter(params,
                                           deck,
                                           eclGrid->getCartesianSize(),
                                           0,
                                           &cartSize[0]));

    // also create an UnstructuredGrid (required to create a BlackoilState)
    Opm::EclipseGridConstPtr constEclGrid(eclGrid);
    ourFineGridManagerPtr.reset(new Opm::GridManager(constEclGrid));

    const UnstructuredGrid &ourFinerUnstructuredGrid = *ourFineGridManagerPtr->c_grid();
    BOOST_CHECK(ourFinerUnstructuredGrid.cartdims[0] == 3);
    BOOST_CHECK(ourFinerUnstructuredGrid.cartdims[1] == 3);
    BOOST_CHECK(ourFinerUnstructuredGrid.cartdims[2] == 3);

    BOOST_CHECK(ourFinerUnstructuredGrid.number_of_cells == 3*3*3);

    // this check is disabled so far, because UnstructuredGrid uses some weird definition
    // of the term "face". For this grid, "number_of_faces" is 108 which is
    // 2*6*numCells...
    //BOOST_CHECK(ourFinerUnstructuredGrid.number_of_faces == 4*4*4);

    int numCells = ourFinerUnstructuredGrid.number_of_cells;
    for (int cellIdx = 0; cellIdx < numCells; ++cellIdx)
        BOOST_CHECK(ourFinerUnstructuredGrid.global_cell[cellIdx] == cellIdx);
}

void createBlackoilState(int timeStepIdx)
{
    // allocate a new BlackoilState object
    const UnstructuredGrid &ourFinerUnstructuredGrid = *ourFineGridManagerPtr->c_grid();
    blackoilState.reset(new Opm::BlackoilState);
    blackoilState->init(ourFinerUnstructuredGrid, 3);

    size_t numCells = ourFinerUnstructuredGrid.number_of_cells;
    size_t numFaces = ourFinerUnstructuredGrid.number_of_faces;

    BOOST_CHECK(blackoilState->pressure().size() == numCells);
    BOOST_CHECK(blackoilState->facepressure().size() == numFaces);
    BOOST_CHECK(blackoilState->faceflux().size() == numFaces);
    BOOST_CHECK(blackoilState->saturation().size() == numCells*3);
    BOOST_CHECK(blackoilState->gasoilratio().size() == numCells);
    BOOST_CHECK(blackoilState->rv().size() == numCells);

    // this check is disabled because BlackoilState does not seem to allocate memory for
    // this field. This means that it is probably unused and unneeded.
    //BOOST_CHECK(blackoilState->surfacevol().size() == numCells*3);

    // fill the state object with some data. The fun with this class is that it does not
    // exhibit a proper c++ way to do this (i.e., getter + setter methods). Instead
    // references to the arrays must be retrieved from the object and manipulated
    // directly. Don't try to call resize() or anything else which is not politically
    // correct on them!
    auto &pressure = blackoilState->pressure();
    auto &facepressure = blackoilState->facepressure();
    auto &faceflux = blackoilState->faceflux();
    auto &saturation = blackoilState->saturation();
    auto &gasoilratio = blackoilState->gasoilratio();
    auto &rv = blackoilState->rv();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx) {
        pressure[cellIdx] = timeStepIdx*1e5 + 1e4 + cellIdx;

        // set the phase saturations. Some fun with direct index manipulation is to be
        // had...
        saturation[3*cellIdx + 0] = timeStepIdx*1e5 +2.1e4 + cellIdx; // oil
        saturation[3*cellIdx + 1] = timeStepIdx*1e5 +2.2e4 + cellIdx; // gas
        saturation[3*cellIdx + 2] = timeStepIdx*1e5 +2.3e4 + cellIdx; // water

        // oil vaporization factor
        rv[cellIdx] = timeStepIdx*1e5 +3e4 + cellIdx;

        // gas dissolution factor
        gasoilratio[cellIdx] = timeStepIdx*1e5 + 4e4 + cellIdx;
    }

    // face specific data
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx) {
        facepressure[faceIdx] = timeStepIdx*1e5 + 5e4 + faceIdx;
        faceflux[faceIdx] = timeStepIdx*1e5 + 6e4 + faceIdx;
    }
}

void createWellState(int timeStepIdx)
{
    // allocate a new BlackoilState object
    wellState.reset(new Opm::WellState);
    wellState->init(0, *blackoilState);
}

void getErtData(ecl_kw_type *eclKeyword, std::vector<double> &data)
{
    size_t kwSize = ecl_kw_get_size(eclKeyword);
    float* ertData = static_cast<float*>(ecl_kw_iget_ptr(eclKeyword, 0));

    data.resize(kwSize);
    std::copy(ertData, ertData + kwSize, data.begin());
}

void getErtData(ecl_kw_type *eclKeyword, std::vector<int> &data)
{
    size_t kwSize = ecl_kw_get_size(eclKeyword);
    int* ertData = static_cast<int*>(ecl_kw_iget_ptr(eclKeyword, 0));

    data.resize(kwSize);
    std::copy(ertData, ertData + kwSize, data.begin());
}

void compareErtData(const std::vector<double> &src, const std::vector<double> &dst, double tolerance)
{
    BOOST_CHECK_EQUAL(src.size(), dst.size());
    if (src.size() != dst.size())
        return;

    for (size_t i = 0; i < src.size(); ++i)
        BOOST_CHECK_CLOSE(src[i], dst[i], tolerance);
}

void compareErtData(const std::vector<int> &src, const std::vector<int> &dst)
{
    BOOST_CHECK_EQUAL(src.size(), dst.size());
    if (src.size() != dst.size())
        return;

    for (size_t i = 0; i < src.size(); ++i)
        BOOST_CHECK_EQUAL(src[i], dst[i]);
}

void checkEgridFile()
{
    size_t numCells = ourFineGridManagerPtr->c_grid()->number_of_cells;

    // use ERT directly to inspect the EGRID file produced by EclipseWriter
    auto egridFile = fortio_open_reader("FOO.EGRID", /*isFormated=*/0, ECL_ENDIAN_FLIP);

    ecl_kw_type *eclKeyword;
    // yes, that's an assignment!
    while ((eclKeyword = ecl_kw_fread_alloc(egridFile))) {
        std::string keywordName(ecl_kw_get_header(eclKeyword));
        if (keywordName == "COORD") {
            std::vector<double> sourceData, resultData;
            eclGrid->exportCOORD(sourceData);
            getErtData(eclKeyword, resultData);
            compareErtData(sourceData, resultData, /*percentTolerance=*/1e-6);
        }
        else if (keywordName == "ZCORN") {
            std::vector<double> sourceData, resultData;
            eclGrid->exportZCORN(sourceData);
            getErtData(eclKeyword, resultData);
            compareErtData(sourceData, resultData, /*percentTolerance=*/1e-6);
        }
        else if (keywordName == "ACTNUM") {
            std::vector<int> sourceData, resultData;
            eclGrid->exportACTNUM(sourceData);
            getErtData(eclKeyword, resultData);

            if (resultData.size() == numCells && sourceData.size() == 0) {
                sourceData.resize(numCells);
                std::fill(sourceData.begin(), sourceData.end(), 1);
            }

            compareErtData(sourceData, resultData);
        }

        ecl_kw_free(eclKeyword);
    }

    fortio_fclose(egridFile);
}

void checkInitFile()
{
    // use ERT directly to inspect the INIT file produced by EclipseWriter
    auto initFile = fortio_open_reader("FOO.INIT", /*isFormated=*/0, ECL_ENDIAN_FLIP);

    ecl_kw_type *eclKeyword;
    // yes, that's an assignment!
    while ((eclKeyword = ecl_kw_fread_alloc(initFile))) {
        std::string keywordName(ecl_kw_get_header(eclKeyword));

        if (keywordName == "PORO") {
            const std::vector<double> &sourceData = deck->getKeyword("PORO")->getSIDoubleData();
            std::vector<double> resultData;
            getErtData(eclKeyword, resultData);

            compareErtData(sourceData, resultData, /*percentTolerance=*/1e-4);
        }

        if (keywordName == "PERMX") {
            std::vector<double> sourceData = deck->getKeyword("PERMX")->getSIDoubleData();
            std::vector<double> resultData;
            getErtData(eclKeyword, resultData);

            // convert the data from ERT from Field to SI units (mD to m^2)
            for (size_t i = 0; i < resultData.size(); ++i) {
                resultData[i] *= 9.869233e-16;
            }

            compareErtData(sourceData, resultData, /*percentTolerance=*/1e-4);
        }

        ecl_kw_free(eclKeyword);
    }

    fortio_fclose(initFile);
}

void checkRestartFile(int timeStepIdx)
{
    size_t numCells = ourFineGridManagerPtr->c_grid()->number_of_cells;

    Opm::PhaseUsage phaseUsage = Opm::phaseUsageFromDeck(deck);
    int numActivePhases = phaseUsage.num_phases;
    int waterPhaseIdx = phaseUsage.phase_pos[Opm::BlackoilPhases::Aqua];
    int gasPhaseIdx = phaseUsage.phase_pos[Opm::BlackoilPhases::Vapour];

    for (int i = 0; i <= timeStepIdx; ++i) {
        createBlackoilState(i);

        // use ERT directly to inspect the restart file produced by EclipseWriter
        auto rstFile = fortio_open_reader("FOO.UNRST", /*isFormated=*/0, ECL_ENDIAN_FLIP);

        int curSeqnum = -1;
        ecl_kw_type *eclKeyword;
        // yes, that's an assignment!
        while ((eclKeyword = ecl_kw_fread_alloc(rstFile))) {
            std::string keywordName(ecl_kw_get_header(eclKeyword));

            if (keywordName == "SEQNUM") {
                curSeqnum = *static_cast<int*>(ecl_kw_iget_ptr(eclKeyword, 0));
            }
            if (curSeqnum != i)
                continue;

            if (keywordName == "PRESSURE") {
                std::vector<double> sourceData = blackoilState->pressure();
                std::vector<double> resultData;
                getErtData(eclKeyword, resultData);

                // convert the data from ERT from Metric to SI units (bar to Pa)
                for (size_t i = 0; i < resultData.size(); ++i) {
                    resultData[i] *= 1e5;
                }

                compareErtData(sourceData, resultData, /*percentTolerance=*/1e-4);
            }

            if (keywordName == "SWAT") {
                std::vector<double> sourceData;
                std::vector<double> resultData;
                getErtData(eclKeyword, resultData);

                // extract the water saturation from the black-oil state
                sourceData.resize(numCells);
                for (size_t i = 0; i < sourceData.size(); ++i) {
                    // again, fun with direct index manipulation...
                    sourceData[i] = blackoilState->saturation()[i*numActivePhases + waterPhaseIdx];
                }

                compareErtData(sourceData, resultData, /*percentTolerance=*/1e-4);
            }

            if (keywordName == "SGAS") {
                std::vector<double> sourceData;
                std::vector<double> resultData;
                getErtData(eclKeyword, resultData);

                // extract the water saturation from the black-oil state
                sourceData.resize(numCells);
                for (size_t i = 0; i < sourceData.size(); ++i) {
                    // again, fun with direct index manipulation...
                    sourceData[i] = blackoilState->saturation()[i*numActivePhases + gasPhaseIdx];
                }

                compareErtData(sourceData, resultData, /*percentTolerance=*/1e-4);
            }
        }

        fortio_fclose(rstFile);
    }
}

void checkSummaryFile(int timeStepIdx)
{
    // TODO
}

BOOST_AUTO_TEST_CASE(EclipseWriterIntegration)
{
    const char *deckString =
        "RUNSPEC\n"
        "OIL\n"
        "GAS\n"
        "WATER\n"
        "METRIC\n"
        "DIMENS\n"
        "3 3 3/\n"
        "GRID\n"
        "DXV\n"
        "1.0 2.0 3.0 /\n"
        "DYV\n"
        "4.0 5.0 6.0 /\n"
        "DZV\n"
        "7.0 8.0 9.0 /\n"
        "TOPS\n"
        "9*100 /\n"
        "PROPS\n"
        "PORO\n"
        "27*0.3 /\n"
        "PERMX\n"
        "27*1 /\n"
        "SCHEDULE\n"
        "TSTEP\n"
        "1.0 2.0 3.0 4.0 /\n"
        "WELSPECS\n"
        "'INJ' 'G' 1 1 2000 'GAS' /\n"
        "'PROD' 'G' 3 3 1000 'OIL' /\n"
        "/\n";

    createEclipseWriter(deckString);

    eclWriter->writeInit(*simTimer);

    checkEgridFile();
    checkInitFile();

    for (; simTimer->currentStepNum() < simTimer->numSteps(); ++ (*simTimer)) {
        createBlackoilState(simTimer->currentStepNum());
        createWellState(simTimer->currentStepNum());
        eclWriter->writeTimeStep(*simTimer, *blackoilState, *wellState);
        checkRestartFile(simTimer->currentStepNum());
        checkSummaryFile(simTimer->currentStepNum());
    }
}
