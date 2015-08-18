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
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/props/IncompPropertiesSinglePhase.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/pressure/IncompTpfaSinglePhase.hpp>
#include <opm/core/flowdiagnostics/FlowDiagnostics.hpp>
#include <opm/core/flowdiagnostics/TofReorder.hpp>
#include <opm/core/flowdiagnostics/TofDiscGalReorder.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseMode.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <memory>
#include <boost/filesystem.hpp>

#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>


namespace
{
    const static double alq_invalid = -std::numeric_limits<double>::max();
    const static int vfp_invalid = -std::numeric_limits<int>::max();

    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }

    void buildTracerheadsFromWells(const Wells& wells,
                                   const bool trace_injectors,
                                   Opm::SparseTable<int>& tracerheads)
    {
        tracerheads.clear();
        const int num_wells = wells.number_of_wells;
        const WellType wanted_type = trace_injectors ? INJECTOR : PRODUCER;
        for (int w = 0; w < num_wells; ++w) {
            if (wells.type[w] != wanted_type) {
                continue;
            }
            tracerheads.appendRow(wells.well_cells + wells.well_connpos[w],
                                  wells.well_cells + wells.well_connpos[w + 1]);
        }
    }

    void setBhpWells(Wells& wells)
    {
        const int num_wells = wells.number_of_wells;
        for (int w = 0; w < num_wells; ++w) {
            WellControls* ctrl = wells.ctrls[w];
            const double target = (wells.type[w] == INJECTOR) ? 200*Opm::unit::barsa : 100*Opm::unit::barsa;
            const double distr[3] = { 1.0, 0.0, 0.0 }; // Large enough irrespective of #phases.
            well_controls_add_new(BHP, target,
                    alq_invalid, vfp_invalid,
                    distr, ctrl);
            well_controls_set_current(ctrl, well_controls_get_num(ctrl) - 1);
        }
    }

    void computeTransportSourceSinglePhase(const UnstructuredGrid& grid,
                                           const std::vector<double>& src,
                                           const std::vector<double>& faceflux,
                                           const double inflow_frac,
                                           const Wells* wells,
                                           const std::vector<double>& well_perfrates,
                                           std::vector<double>& transport_src)
    {
        using namespace Opm;
        int nc = grid.number_of_cells;
        transport_src.resize(nc);
        // Source term and boundary contributions.
        for (int c = 0; c < nc; ++c) {
            transport_src[c] = 0.0;
            transport_src[c] += src[c] > 0.0 ? inflow_frac*src[c] : src[c];
            for (int hf = grid.cell_facepos[c]; hf < grid.cell_facepos[c + 1]; ++hf) {
                int f = grid.cell_faces[hf];
                const int* f2c = &grid.face_cells[2*f];
                double bdy_influx = 0.0;
                if (f2c[0] == c && f2c[1] == -1) {
                    bdy_influx = -faceflux[f];
                } else if (f2c[0] == -1 && f2c[1] == c) {
                    bdy_influx = faceflux[f];
                }
                if (bdy_influx != 0.0) {
                    transport_src[c] += bdy_influx > 0.0 ? inflow_frac*bdy_influx : bdy_influx;
                }
            }
        }

        // Well contributions.
        if (wells) {
            const int nw = wells->number_of_wells;
            for (int w = 0; w < nw; ++w) {
                for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                    const int perf_cell = wells->well_cells[perf];
                    double perf_rate = well_perfrates[perf];
                    if (perf_rate > 0.0) {
                        // perf_rate is a total inflow rate, we want a water rate.
                        if (wells->type[w] != INJECTOR) {
                            std::cout << "**** Warning: crossflow in well "
                                      << w << " perf " << perf - wells->well_connpos[w]
                                      << " ignored. Rate was "
                                      << perf_rate/Opm::unit::day << " m^3/day." << std::endl;
                            perf_rate = 0.0;
                        } else {
                            perf_rate *= inflow_frac;
                        }
                    }
                    transport_src[perf_cell] += perf_rate;
                }
            }
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
    parameter::ParameterGroup param(argc, argv);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Read the deck.
    std::string deck_filename = param.get<std::string>("deck_filename");
    Parser parser;
    ParseMode parseMode;
    DeckConstPtr deck = parser.parseFile(deck_filename , parseMode);
    EclipseStateConstPtr eclipseState = std::make_shared<EclipseState>(deck , parseMode);

    // Grid init
    GridManager grid_manager(deck);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    // Rock and fluid init
    IncompPropertiesSinglePhase props(deck, eclipseState, grid);
    // Wells init.
    WellsManager wells_manager(eclipseState , 0, grid, props.permeability());

    std::shared_ptr<Wells> my_wells(clone_wells(wells_manager.c_wells()), destroy_wells);
    setBhpWells(*my_wells);
    const Wells& wells = *my_wells;

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
    computeTransportSourceSinglePhase(grid, src, flux, 1.0,
                                      &wells, wellrates, transport_src);

    std::string tof_filenames[2] = { output_dir + "/ftof.txt", output_dir + "/btof.txt" };
    std::string tracer_filenames[2] = { output_dir + "/ftracer.txt", output_dir + "/btracer.txt" };
    std::vector<double> tracers[2];

    // We compute tof twice, direction == 0 is from injectors, 1 is from producers.
    for (int direction = 0; direction < 2; ++direction) {
        // Turn direction of flux and flip source terms if starting from producers.
        if (direction == 1) {
            for (auto it = flux.begin(); it != flux.end(); ++it) {
                (*it) = -(*it);
            }
            for (auto it = transport_src.begin(); it != transport_src.end(); ++it) {
                (*it) = -(*it);
            }
        }

        // Solve time-of-flight.
        transport_timer.start();
        std::vector<double> tof;
        std::vector<double> tracer;
        Opm::SparseTable<int> tracerheads;
        if (compute_tracer) {
            buildTracerheadsFromWells(wells, direction == 0, tracerheads);
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
        if (direction == 0) {
            std::cout << "Forward ";
        } else {
            std::cout << "Backward ";
        }
        std::cout << "time-of-flight/tracer solve took: " << tt << " seconds." << std::endl;
        ttime += tt;

        // Output.
        if (output) {
            std::string tof_filename = tof_filenames[direction];
            std::ofstream tof_stream(tof_filename.c_str());
            tof_stream.precision(16);
            std::copy(tof.begin(), tof.end(), std::ostream_iterator<double>(tof_stream, "\n"));
            if (compute_tracer) {
                std::string tracer_filename = tracer_filenames[direction];
                std::ofstream tracer_stream(tracer_filename.c_str());
                tracer_stream.precision(16);
                const int nt = tracer.size()/num_cells;
                for (int i = 0; i < nt*num_cells; ++i) {
                    tracer_stream << tracer[i] << (((i + 1) % nt == 0) ? '\n' : ' ');
                }
                tracers[direction] = tracer;
            }
        }
    }

    // If we have tracers, compute well pairs.
    if (compute_tracer) {
        auto wp = Opm::computeWellPairs(wells, porevol, tracers[0], tracers[1]);
        std::string wellpair_filename = output_dir + "/wellpairs.txt";
        std::ofstream wellpair_stream(wellpair_filename.c_str());
        const int nwp = wp.size();
        for (int ii = 0; ii < nwp; ++ii) {
            wellpair_stream << std::get<0>(wp[ii]) << ' ' << std::get<1>(wp[ii]) << ' ' << std::get<2>(wp[ii]) << '\n';
        }
    }

    total_timer.stop();

    std::cout << "\n\n================    End of simulation     ===============\n"
              << "Total time taken:  " << total_timer.secsSinceStart()
              << "\n  Pressure time:   " << ptime
              << "\n  Tof/tracer time: " << ttime << std::endl;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
