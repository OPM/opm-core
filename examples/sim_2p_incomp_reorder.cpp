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


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/initState.hpp>
#include <opm/core/utility/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorTwophase.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <tr1/array>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <vector>
#include <numeric>




// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    using namespace Opm;

    std::cout << "\n================    Test program for incompressible two-phase flow     ===============\n\n";
    Opm::parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<Opm::GridManager> grid;
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
    boost::scoped_ptr<Opm::WellsManager> wells;
    boost::scoped_ptr<Opm::RockCompressibility> rock_comp;
    Opm::SimulatorTimer simtimer;
    Opm::TwophaseState state;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    if (use_deck) {
        std::string deck_filename = param.get<std::string>("deck_filename");
        Opm::EclipseGridParser deck(deck_filename);
        // Grid init
        grid.reset(new Opm::GridManager(deck));
        // Rock and fluid init
        const int* gc = grid->c_grid()->global_cell;
        std::vector<int> global_cell(gc, gc + grid->c_grid()->number_of_cells);
        props.reset(new Opm::IncompPropertiesFromDeck(deck, global_cell));
        // Wells init.
        wells.reset(new Opm::WellsManager(deck, *grid->c_grid(), props->permeability()));
        // check_well_controls = param.getDefault("check_well_controls", false);
        // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
        // Timer init.
        if (deck.hasField("TSTEP")) {
            simtimer.init(deck);
        } else {
            simtimer.init(param);
        }
        // Rock compressibility.
        rock_comp.reset(new Opm::RockCompressibility(deck));
        // Gravity.
        gravity[2] = deck.hasField("NOGRAV") ? 0.0 : Opm::unit::gravity;
        // Init state variables (saturation and pressure).
        if (param.has("init_saturation")) {
            initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
        } else {
            initStateFromDeck(*grid->c_grid(), *props, deck, gravity[2], state);
        }
    } else {
        // Grid init.
        const int nx = param.getDefault("nx", 100);
        const int ny = param.getDefault("ny", 100);
        const int nz = param.getDefault("nz", 1);
        const double dx = param.getDefault("dx", 1.0);
        const double dy = param.getDefault("dy", 1.0);
        const double dz = param.getDefault("dz", 1.0);
        grid.reset(new Opm::GridManager(nx, ny, nz, dx, dy, dz));
        // Rock and fluid init.
        props.reset(new Opm::IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        // Wells init.
        wells.reset(new Opm::WellsManager());
        // Timer init.
        simtimer.init(param);
        // Rock compressibility.
        rock_comp.reset(new Opm::RockCompressibility(param));
        // Gravity.
        gravity[2] = param.getDefault("gravity", 0.0);
        // Init state variables (saturation and pressure).
        initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
    }

    // Warn if gravity but no density difference.
    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    if (use_gravity) {
        if (props->density()[0] == props->density()[1]) {
            std::cout << "**** Warning: nonzero gravity, but zero density difference." << std::endl;
        }
    }

    // Source-related variables init.
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> totmob;
    std::vector<double> omega; // Will remain empty if no gravity.
    std::vector<double> rc; // Will remain empty if no rock compressibility.

    // Extra rock init.
    std::vector<double> porevol;
    if (rock_comp->isActive()) {
        computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, state.pressure(), porevol);
    } else {
        computePorevolume(*grid->c_grid(), props->porosity(), porevol);
    }
    double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Initialising src
    std::vector<double> src(num_cells, 0.0);
    if (wells->c_wells()) {
        // Do nothing, wells will be the driving force, not source terms.
        // Opm::wellsToSrc(*wells->c_wells(), num_cells, src);
    } else {
        const double default_injection = use_gravity ? 0.0 : 0.1;
        const double flow_per_sec = param.getDefault<double>("injected_porevolumes_per_day", default_injection)
            *tot_porevol_init/Opm::unit::day;
        src[0] = flow_per_sec;
        src[num_cells - 1] = -flow_per_sec;
    }

    // Boundary conditions.
    Opm::FlowBCManager bcs;
    if (param.getDefault("use_pside", false)) {
        int pside = param.get<int>("pside");
        double pside_pressure = param.get<double>("pside_pressure");
        bcs.pressureSide(*grid->c_grid(), Opm::FlowBCManager::Side(pside), pside_pressure);
    }

    // Linear solver.
    Opm::LinearSolverFactory linsolver(param);

    const double *grav = use_gravity ? &gravity[0] : 0;


    Opm::SimulatorTwophase simulator(param,
                                     *grid->c_grid(),
                                     *props,
                                     rock_comp->isActive() ? rock_comp.get() : 0,
                                     wells->c_wells(),
                                     src,
                                     bcs.c_bcs(),
                                     linsolver,
                                     grav);

    // Warn if any parameters are unused.
    if (param.anyUnused()) {
        std::cout << "--------------------   Unused parameters:   --------------------\n";
        param.displayUsage();
        std::cout << "----------------------------------------------------------------" << std::endl;
    }

    // Write parameters used for later reference.
    // if (output) {
    //     param.writeParam(output_dir + "/spu_2p.param");
    // }

    WellState well_state;
    well_state.init(wells->c_wells(), state);

    simulator.run(simtimer, state, well_state);
}
