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

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/SparseTable.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/props/IncompPropertiesSinglePhase.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/pressure/IncompTpfaSinglePhase.hpp>
#include <opm/core/tof/TofReorder.hpp>
#include <opm/core/tof/TofDiscGalReorder.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <memory>
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

    void buildTracerheadsFromWells(const Wells& wells,
                                   Opm::SparseTable<int>& tracerheads)
    {
        tracerheads.clear();
        const int num_wells = wells.number_of_wells;
        for (int w = 0; w < num_wells; ++w) {
            if (wells.type[w] != INJECTOR) {
                continue;
            }
            tracerheads.appendRow(wells.well_cells + wells.well_connpos[w],
                                  wells.well_cells + wells.well_connpos[w + 1]);
        }
    }

} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;

    std::cout << "\n================    Test program for incompressible tof computations     ===============\n\n";
    parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Read the deck.
    std::string deck_filename = param.get<std::string>("deck_filename");
    Parser parser;
    DeckConstPtr deck = parser.parseFile(deck_filename);
    EclipseStateConstPtr eclipseState = std::make_shared<EclipseState>(deck);

    // Grid init
    GridManager grid_manager(deck);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    // Rock and fluid init
    IncompPropertiesSinglePhase props(deck, eclipseState, grid);
    // Wells init.
    WellsManager wells_manager(eclipseState , 0, grid, props.permeability());
    const Wells& wells = *wells_manager.c_wells();

    // Pore volume.
    std::vector<double> porevol;
    computePorevolume(grid, props.porosity(), porevol);
    int num_cells = grid.number_of_cells;

    // Linear solver.
    LinearSolverFactory linsolver(param);

    // Pressure solver.
    Opm::IncompTpfaSinglePhase psolver(grid, props, linsolver, wells);

    // Choice of tof solver.
    bool use_dg = param.getDefault("use_dg", false);
    bool use_multidim_upwind = false;
    // Need to initialize dg solver here, since it uses parameters now.
    std::unique_ptr<Opm::TofDiscGalReorder> dg_solver;
    if (use_dg) {
        dg_solver.reset(new Opm::TofDiscGalReorder(grid, param));
    } else {
        use_multidim_upwind = param.getDefault("use_multidim_upwind", false);
    }
    bool compute_tracer = param.getDefault("compute_tracer", false);

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
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        std::string filename = output_dir + "/epoch_timing.param";
        epoch_os.open(filename.c_str(), std::fstream::trunc | std::fstream::out);
        // open file to clean it. The file is appended to in SimulatorTwophase
        filename = output_dir + "/step_timing.param";
        std::fstream step_os(filename.c_str(), std::fstream::trunc | std::fstream::out);
        step_os.close();
        param.writeParam(output_dir + "/simulation.param");
    }

    // Check if we have misspelled anything
    warnIfUnusedParams(param);

    // Main solvers.
    Opm::time::StopWatch pressure_timer;
    double ptime = 0.0;
    Opm::time::StopWatch transport_timer;
    double ttime = 0.0;
    Opm::time::StopWatch total_timer;
    total_timer.start();
    std::cout << "\n\n================    Starting main solvers     ===============" << std::endl;

    // Solve pressure.
    std::vector<double> press;
    std::vector<double> flux;
    std::vector<double> bhp;
    std::vector<double> wellrates;
    pressure_timer.start();
    psolver.solve(press, flux, bhp, wellrates);
    pressure_timer.stop();
    double pt = pressure_timer.secsSinceStart();
    std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
    ptime += pt;

    // Process transport sources (to include bdy terms and well flows).
    std::vector<double> src(num_cells, 0.0);
    std::vector<double> transport_src;
    Opm::computeTransportSource(grid, src, flux, 1.0,
                                &wells, wellrates, transport_src);

    // Solve time-of-flight.
    transport_timer.start();
    std::vector<double> tof;
    std::vector<double> tracer;
    Opm::SparseTable<int> tracerheads;
    if (compute_tracer) {
        buildTracerheadsFromWells(wells, tracerheads);
    }
    if (use_dg) {
        if (compute_tracer) {
            dg_solver->solveTofTracer(flux.data(), porevol.data(), transport_src.data(), tracerheads, tof, tracer);
        } else {
            dg_solver->solveTof(flux.data(), porevol.data(), transport_src.data(), tof);
        }
    } else {
        Opm::TofReorder tofsolver(grid, use_multidim_upwind);
        if (compute_tracer) {
            tofsolver.solveTofTracer(flux.data(), porevol.data(), transport_src.data(), tracerheads, tof, tracer);
        } else {
            tofsolver.solveTof(flux.data(), porevol.data(), transport_src.data(), tof);
        }
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
        tof_stream.precision(16);
        std::copy(tof.begin(), tof.end(), std::ostream_iterator<double>(tof_stream, "\n"));
        if (compute_tracer) {
            std::string tracer_filename = output_dir + "/tracer.txt";
            std::ofstream tracer_stream(tracer_filename.c_str());
            tracer_stream.precision(16);
            const int nt = tracer.size()/num_cells;
            for (int i = 0; i < nt*num_cells; ++i) {
                tracer_stream << tracer[i] << (((i + 1) % nt == 0) ? '\n' : ' ');
            }
        }
    }

    std::cout << "\n\n================    End of simulation     ===============\n"
              << "Total time taken: " << total_timer.secsSinceStart()
              << "\n  Pressure time:  " << ptime
              << "\n  Transport time: " << ttime << std::endl;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
