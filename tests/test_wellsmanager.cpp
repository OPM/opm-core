/*
  Copyright 2013 Statoil ASA.

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

#define BOOST_TEST_MODULE WellsManagerTests
#include <boost/test/unit_test.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>

#include <opm/core/grid/GridManager.hpp>

void wells_static_check(const Wells* wells) {
    BOOST_CHECK_EQUAL(2, wells->number_of_wells);
    BOOST_CHECK_EQUAL(3, wells->number_of_phases);

    BOOST_CHECK_EQUAL("INJ1", wells->name[0]);
    BOOST_CHECK_EQUAL("PROD1", wells->name[1]);

    /* The mapping from well number into the wells->WI and wells->well_cells arrays. */
    BOOST_CHECK_EQUAL(0, wells->well_connpos[0]);
    BOOST_CHECK_EQUAL(1, wells->well_connpos[1]);
    BOOST_CHECK_EQUAL(2, wells->well_connpos[2]);

    /* Connection factor */
    BOOST_CHECK_CLOSE(1.2279166666666664e-12, wells->WI[0], 0.001);
    BOOST_CHECK_CLOSE(1.2279166666666664e-12, wells->WI[1], 0.001);

    /* Completed cells */
    BOOST_CHECK_EQUAL(0, wells->well_cells[0]);
    BOOST_CHECK_EQUAL(9 + 2 * 10 + 2 * 10 * 10, wells->well_cells[1]);
}

/*
 The number of controls is determined by looking at which elements
 have been given explicit - non-default - values in the WCONxxxx
 keyword. Is that at all interesting?
 */

void check_controls_epoch0(struct WellControls ** ctrls) {
    // The injector
    {
        const struct WellControls * ctrls0 = ctrls[0];
        BOOST_CHECK_EQUAL(3, well_controls_get_num(ctrls0));   // The number of controls for the injector == 3??

        BOOST_CHECK_EQUAL(SURFACE_RATE, well_controls_iget_type(ctrls0, 0));
        BOOST_CHECK_EQUAL(RESERVOIR_RATE, well_controls_iget_type(ctrls0, 1));
        BOOST_CHECK_EQUAL(BHP, well_controls_iget_type(ctrls0, 2));

        // The different targets
        BOOST_CHECK_EQUAL(100.0 / 86400, well_controls_iget_target(ctrls0, 0));
        BOOST_CHECK_EQUAL(200.0 / 86400, well_controls_iget_target(ctrls0, 1));
        BOOST_CHECK_EQUAL(400 * 100000, well_controls_iget_target(ctrls0, 2));

        // Which control is active
        BOOST_CHECK_EQUAL(0, well_controls_get_current(ctrls0));

        // The phase distribution in the active target
        {
            const double * distr = well_controls_iget_distr(ctrls0, 0);
            BOOST_CHECK_EQUAL(0, distr[0]);  // Water
            BOOST_CHECK_EQUAL(0, distr[1]);  // Oil
            BOOST_CHECK_EQUAL(1, distr[2]);  // Gas
        }
    }

    // The producer
    {
        const struct WellControls * ctrls1 = ctrls[1];
        BOOST_CHECK_EQUAL(2, well_controls_get_num(ctrls1));   // The number of controls for the producer == 2??
        BOOST_CHECK_EQUAL(SURFACE_RATE, well_controls_iget_type(ctrls1, 0));
        BOOST_CHECK_EQUAL(BHP, well_controls_iget_type(ctrls1, 1));

        // The different targets
        BOOST_CHECK_EQUAL(-20000.0 / 86400, well_controls_iget_target(ctrls1, 0));
        BOOST_CHECK_EQUAL(1000 * 100000, well_controls_iget_target(ctrls1, 1));

        // Which control is active
        BOOST_CHECK_EQUAL(0, well_controls_get_current(ctrls1));

        // The phase distribution in the active target
        {
            const double * distr = well_controls_iget_distr(ctrls1, 0);
            BOOST_CHECK_EQUAL(0, distr[0]);  // Water
            BOOST_CHECK_EQUAL(1, distr[1]);  // Oil
            BOOST_CHECK_EQUAL(0, distr[2]);  // Gas
        }
    }
}

void check_controls_epoch1(struct WellControls ** ctrls) {
    // The injector
    {
        const struct WellControls * ctrls0 = ctrls[0];
        BOOST_CHECK_EQUAL(3, well_controls_get_num(ctrls0));   // The number of controls for the injector == 3??

        BOOST_CHECK_EQUAL(SURFACE_RATE, well_controls_iget_type(ctrls0, 0));
        BOOST_CHECK_EQUAL(RESERVOIR_RATE, well_controls_iget_type(ctrls0, 1));
        BOOST_CHECK_EQUAL(BHP, well_controls_iget_type(ctrls0, 2));

        // The different targets
        BOOST_CHECK_CLOSE(10.0 / 86400, well_controls_iget_target(ctrls0, 0), 0.001);
        BOOST_CHECK_CLOSE(20.0 / 86400, well_controls_iget_target(ctrls0, 1), 0.001);
        BOOST_CHECK_CLOSE(40 * 100000, well_controls_iget_target(ctrls0, 2), 0.001);

        // Which control is active
        BOOST_CHECK_EQUAL(1, well_controls_get_current(ctrls0));

        {
            const double * distr = well_controls_iget_distr(ctrls0, 1);
            BOOST_CHECK_EQUAL(1, distr[0]);  // Water
            BOOST_CHECK_EQUAL(0, distr[1]);  // Oil
            BOOST_CHECK_EQUAL(0, distr[2]);  // Gas
        }
    }

    // The producer
    {
        const struct WellControls * ctrls1 = ctrls[1];
        BOOST_CHECK_EQUAL(3, well_controls_get_num(ctrls1));   // The number of controls for the producer - now 3.
        BOOST_CHECK_EQUAL(SURFACE_RATE, well_controls_iget_type(ctrls1, 0));
        BOOST_CHECK_EQUAL(RESERVOIR_RATE, well_controls_iget_type(ctrls1, 1));
        BOOST_CHECK_EQUAL(BHP, well_controls_iget_type(ctrls1, 2));

        // The different targets
        BOOST_CHECK_CLOSE(-999.0 / 86400, well_controls_iget_target(ctrls1, 0), 0.001);
        BOOST_CHECK_CLOSE(-123.0 / 86400, well_controls_iget_target(ctrls1, 1), 0.001);
        BOOST_CHECK_CLOSE(100 * 100000, well_controls_iget_target(ctrls1, 2), 0.001);

        // Which control is active
        BOOST_CHECK_EQUAL(1, well_controls_get_current(ctrls1));

        {
            const double * distr = well_controls_iget_distr(ctrls1, 1);
            BOOST_CHECK_EQUAL(1, distr[0]);  // Water
            BOOST_CHECK_EQUAL(1, distr[1]);  // Oil
            BOOST_CHECK_EQUAL(1, distr[2]);  // Gas
        }
    }
}

void check_controls_epoch3(struct WellControls ** ctrls) {
    // The new producer
    const struct WellControls * ctrls1 = ctrls[1];
    // Note: controls include default (1 atm) BHP control.
    // and also an ORAT control
    BOOST_CHECK_EQUAL(2, well_controls_get_num(ctrls1));
}

BOOST_AUTO_TEST_CASE(New_Constructor_Works) {

    const std::string filename = "wells_manager_data.data";
    Opm::Parser parser;
    Opm::ParseContext parseContext;
    Opm::Deck deck = parser.parseFile(filename, parseContext);

    Opm::EclipseState eclipseState(deck, parseContext);
    Opm::GridManager gridManager(eclipseState.getInputGrid());
    const auto& grid = eclipseState.getInputGrid();
    const Opm::TableManager table ( deck );
    const Opm::Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Schedule sched(deck, grid, eclipseProperties, Opm::Phases(true, true, true), parseContext );

 
    {
        Opm::WellsManager wellsManager(eclipseState, sched, 0, *gridManager.c_grid());
        wells_static_check(wellsManager.c_wells());
        check_controls_epoch0(wellsManager.c_wells()->ctrls);
    }

    {
        Opm::WellsManager wellsManager(eclipseState, sched, 1, *gridManager.c_grid());
        wells_static_check(wellsManager.c_wells());
        check_controls_epoch1(wellsManager.c_wells()->ctrls);
    }

    {
        Opm::WellsManager wellsManager(eclipseState, sched, 3, *gridManager.c_grid());
        const Wells* wells = wellsManager.c_wells();

        // There is 3 wells in total in the deck at the 3rd schedule step.
        // PROD1 is shut and should therefore not be counted.
        // The new well is therefore the secound well.
        BOOST_CHECK_EQUAL(2, wells->number_of_wells);
        BOOST_CHECK_EQUAL(wells->name[0], "INJ1");
        BOOST_CHECK_EQUAL(wells->name[1], "NEW");

        check_controls_epoch3(wellsManager.c_wells()->ctrls);
    }

}

BOOST_AUTO_TEST_CASE(WellsEqual) {
    const std::string filename = "wells_manager_data.data";
    Opm::ParseContext parseContext;
    Opm::Parser parser;
    Opm::Deck deck(parser.parseFile(filename, parseContext));
    Opm::EclipseState eclipseState(deck, parseContext);
    Opm::GridManager gridManager(eclipseState.getInputGrid());
    const auto& grid = eclipseState.getInputGrid();
    const Opm::TableManager table ( deck );
    const Opm::Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Schedule sched(deck, grid, eclipseProperties, Opm::Phases(true, true, true), parseContext );


    Opm::WellsManager wellsManager0(eclipseState, sched, 0, *gridManager.c_grid());
    Opm::WellsManager wellsManager1(eclipseState, sched, 1, *gridManager.c_grid());

    BOOST_CHECK(wells_equal( wellsManager0.c_wells() , wellsManager0.c_wells(),false));
    BOOST_CHECK(!wells_equal( wellsManager0.c_wells() , wellsManager1.c_wells(),false));
}

BOOST_AUTO_TEST_CASE(ControlsEqual) {
    const std::string filename = "wells_manager_data.data";
    Opm::ParseContext parseContext;
    Opm::Parser parser;
    Opm::Deck deck(parser.parseFile(filename, parseContext));
    Opm::EclipseState eclipseState(deck, parseContext);
    Opm::GridManager gridManager(eclipseState.getInputGrid());
    const auto& grid = eclipseState.getInputGrid();
    const Opm::TableManager table ( deck );
    const Opm::Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Schedule sched(deck, grid, eclipseProperties, Opm::Phases(true, true, true), parseContext );



    Opm::WellsManager wellsManager0(eclipseState, sched, 0, *gridManager.c_grid());
    Opm::WellsManager wellsManager1(eclipseState, sched, 1, *gridManager.c_grid());

    BOOST_CHECK(well_controls_equal( wellsManager0.c_wells()->ctrls[0] , wellsManager0.c_wells()->ctrls[0] , false));
    BOOST_CHECK(well_controls_equal( wellsManager0.c_wells()->ctrls[1] , wellsManager0.c_wells()->ctrls[1] , false));
    BOOST_CHECK(well_controls_equal( wellsManager1.c_wells()->ctrls[0] , wellsManager1.c_wells()->ctrls[0] , false));
    BOOST_CHECK(well_controls_equal( wellsManager1.c_wells()->ctrls[1] , wellsManager1.c_wells()->ctrls[1] , false));

    BOOST_CHECK(!well_controls_equal( wellsManager0.c_wells()->ctrls[0] , wellsManager0.c_wells()->ctrls[1] , false));
    BOOST_CHECK(!well_controls_equal( wellsManager0.c_wells()->ctrls[1] , wellsManager0.c_wells()->ctrls[0] , false));
    BOOST_CHECK(!well_controls_equal( wellsManager1.c_wells()->ctrls[0] , wellsManager0.c_wells()->ctrls[0] , false));
    BOOST_CHECK(!well_controls_equal( wellsManager1.c_wells()->ctrls[1] , wellsManager0.c_wells()->ctrls[1] , false));
}

BOOST_AUTO_TEST_CASE(WellShutOK) {
    const std::string filename = "wells_manager_data.data";
    Opm::ParseContext parseContext;
    Opm::Parser parser;
    Opm::Deck deck(parser.parseFile(filename, parseContext));
    Opm::EclipseState eclipseState(deck, parseContext);
    Opm::GridManager gridManager(eclipseState.getInputGrid());
    const auto& grid = eclipseState.getInputGrid();
    const Opm::TableManager table ( deck );
    const Opm::Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Schedule sched(deck, grid, eclipseProperties, Opm::Phases(true, true, true), parseContext );


    Opm::WellsManager wellsManager2(eclipseState, sched, 2, *gridManager.c_grid());

    // Shut wells are not added to the deck. i.e number of wells should be 2-1
    BOOST_CHECK(wellsManager2.c_wells()->number_of_wells == 1);

    //BOOST_CHECK_NO_THROW( Opm::WellsManager wellsManager2(eclipseState , 2 , *gridManager.c_grid(), NULL));
}

BOOST_AUTO_TEST_CASE(WellSTOPOK) {
    const std::string filename = "wells_manager_data_wellSTOP.data";
    Opm::ParseContext parseContext;
    Opm::Parser parser;
    Opm::Deck deck(parser.parseFile(filename, parseContext));
    Opm::EclipseState eclipseState(deck, parseContext);
    Opm::GridManager gridManager(eclipseState.getInputGrid());
    const auto& grid = eclipseState.getInputGrid();
    const Opm::TableManager table ( deck );
    const Opm::Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Schedule sched(deck, grid, eclipseProperties, Opm::Phases(true, true, true), parseContext );



    Opm::WellsManager wellsManager(eclipseState, sched, 0, *gridManager.c_grid());

    const Wells* wells = wellsManager.c_wells();
    const struct WellControls* ctrls0 = wells->ctrls[0];
    const struct WellControls* ctrls1 = wells->ctrls[1];

    BOOST_CHECK(well_controls_well_is_stopped(ctrls0)); // The first well is closed
    BOOST_CHECK(well_controls_well_is_open(ctrls1));  // The second well is open
}
