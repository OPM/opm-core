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

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/initState.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/transport/reorder/TransportModelTracerTof.hpp>
#include <opm/core/transport/reorder/TransportModelTracerTofDiscGal.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/filesystem.hpp>

#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>


namespace
{
    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }
} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    using namespace Opm;

    std::cout << "\n================    Test program for incompressible tof computations     ===============\n\n";
    parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<EclipseGridParser> deck;
    boost::scoped_ptr<GridManager> grid;
    boost::scoped_ptr<IncompPropertiesInterface> props;
    boost::scoped_ptr<Opm::WellsManager> wells;
    TwophaseState state;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    if (use_deck) {
        std::string deck_filename = param.get<std::string>("deck_filename");
        deck.reset(new EclipseGridParser(deck_filename));
        // Grid init
        grid.reset(new GridManager(*deck));
        // Rock and fluid init
        props.reset(new IncompPropertiesFromDeck(*deck, *grid->c_grid()));
        // Wells init.
        wells.reset(new Opm::WellsManager(*deck, *grid->c_grid(), props->permeability()));
        // Gravity.
        gravity[2] = deck->hasField("NOGRAV") ? 0.0 : unit::gravity;
        // Init state variables (saturation and pressure).
        if (param.has("init_saturation")) {
            initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
        } else {
            initStateFromDeck(*grid->c_grid(), *props, *deck, gravity[2], state);
        }
    } else {
        // Grid init.
        const int nx = param.getDefault("nx", 100);
        const int ny = param.getDefault("ny", 100);
        const int nz = param.getDefault("nz", 1);
        const double dx = param.getDefault("dx", 1.0);
        const double dy = param.getDefault("dy", 1.0);
        const double dz = param.getDefault("dz", 1.0);
        grid.reset(new GridManager(nx, ny, nz, dx, dy, dz));
        // Rock and fluid init.
        props.reset(new IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        // Wells init.
        wells.reset(new Opm::WellsManager());
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
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Initialising src
    std::vector<double> porevol;
    computePorevolume(*grid->c_grid(), props->porosity(), porevol);
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> src(num_cells, 0.0);
    if (use_deck) {
        // Do nothing, wells will be the driving force, not source terms.
    } else {
        const double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        const double default_injection = use_gravity ? 0.0 : 0.1;
        const double flow_per_sec = param.getDefault<double>("injected_porevolumes_per_day", default_injection)
            *tot_porevol_init/unit::day;
        src[0] = flow_per_sec;
        src[num_cells - 1] = -flow_per_sec;
    }

    // Boundary conditions.
    FlowBCManager bcs;
    if (param.getDefault("use_pside", false)) {
        int pside = param.get<int>("pside");
        double pside_pressure = param.get<double>("pside_pressure");
        bcs.pressureSide(*grid->c_grid(), FlowBCManager::Side(pside), pside_pressure);
    }

    // Linear solver.
    LinearSolverFactory linsolver(param);

    // Pressure solver.
    Opm::IncompTpfa psolver(*grid->c_grid(), *props, 0, linsolver,
                            0.0, 0.0, 0,
                            grav, wells->c_wells(), src, bcs.c_bcs());

    // Choice of tof solver.
    bool use_dg = param.getDefault("use_dg", false);
    int dg_degree = -1;
    bool use_cvi = false;
    bool use_limiter = false;
    bool use_multidim_upwind = false;
    if (use_dg) {
        dg_degree = param.getDefault("dg_degree", 0);
        use_cvi = param.getDefault("use_cvi", false);
        use_limiter = param.getDefault("use_limiter", false);
    } else {
        use_multidim_upwind = param.getDefault("use_multidim_upwind", false);
    }

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    std::ofstream epoch_os;
    std::string output_dir;
    if (output) {
        output_dir =
            param.getDefault("output_dir", std::string("output"));
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            THROW("Creating directories failed: " << fpath);
        }
        std::string filename = output_dir + "/epoch_timing.param";
        epoch_os.open(filename.c_str(), std::fstream::trunc | std::fstream::out);
        // open file to clean it. The file is appended to in SimulatorTwophase
        filename = output_dir + "/step_timing.param";
        std::fstream step_os(filename.c_str(), std::fstream::trunc | std::fstream::out);
        step_os.close();
        param.writeParam(output_dir + "/simulation.param");
    }

    // Init wells.
    Opm::WellState well_state;
    well_state.init(wells->c_wells(), state);

    // Main solvers.
    Opm::time::StopWatch pressure_timer;
    double ptime = 0.0;
    Opm::time::StopWatch transport_timer;
    double ttime = 0.0;
    Opm::time::StopWatch total_timer;
    total_timer.start();
    std::cout << "\n\n================    Starting main solvers     ===============" << std::endl;

    // Solve pressure.
    pressure_timer.start();
    psolver.solve(1.0, state, well_state);
    pressure_timer.stop();
    double pt = pressure_timer.secsSinceStart();
    std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
    ptime += pt;

    // Process transport sources (to include bdy terms and well flows).
    std::vector<double> transport_src;
    Opm::computeTransportSource(*grid->c_grid(), src, state.faceflux(), 1.0,
                                wells->c_wells(), well_state.perfRates(), transport_src);

    // Solve time-of-flight.
    transport_timer.start();
    std::vector<double> tof;
    if (use_dg) {
        Opm::TransportModelTracerTofDiscGal tofsolver(*grid->c_grid(), use_cvi, use_limiter);
        tofsolver.solveTof(&state.faceflux()[0], &porevol[0], &transport_src[0], dg_degree, tof);
    } else {
        Opm::TransportModelTracerTof tofsolver(*grid->c_grid(), use_multidim_upwind);
        tofsolver.solveTof(&state.faceflux()[0], &porevol[0], &transport_src[0], tof);
    }
    transport_timer.stop();
    double tt = transport_timer.secsSinceStart();
    std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
    ttime += tt;
    total_timer.stop();

    // Output.
    if (output) {
        std::string tof_filename = output_dir + "/tof.txt";
        std::ofstream tof_stream(tof_filename.c_str());
        std::copy(tof.begin(), tof.end(), std::ostream_iterator<double>(tof_stream, "\n"));
    }

    std::cout << "\n\n================    End of simulation     ===============\n"
              << "Total time taken: " << total_timer.secsSinceStart()
              << "\n  Pressure time:  " << ptime
              << "\n  Transport time: " << ttime << std::endl;
}
