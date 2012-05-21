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

#include <opm/core/pressure/CompressibleTpfa.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/initState.hpp>
#include <opm/core/utility/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/BlackoilPropertiesBasic.hpp>
#include <opm/core/fluid/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/ColumnExtract.hpp>
#include <opm/core/BlackoilState.hpp>
#include <opm/core/WellState.hpp>
#include <opm/core/transport/GravityColumnSolver.hpp>

#include <opm/core/transport/reorder/TransportModelTwophase.hpp>

#include <boost/filesystem/convenience.hpp>
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

#define TRANSPORT_SOLVER_FIXED 0


template <class State>
static void outputState(const UnstructuredGrid& grid,
                        const State& state,
                        const int step,
                        const std::string& output_dir)
{
    // Write data in VTK format.
    std::ostringstream vtkfilename;
    vtkfilename << output_dir << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
        THROW("Failed to open " << vtkfilename.str());
    }
    Opm::DataMap dm;
    dm["saturation"] = &state.saturation();
    dm["pressure"] = &state.pressure();
    std::vector<double> cell_velocity;
    Opm::estimateCellVelocity(grid, state.faceflux(), cell_velocity);
    dm["velocity"] = &cell_velocity;
    Opm::writeVtkData(grid, dm, vtkfile);

    // Write data (not grid) in Matlab format
    for (Opm::DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
        std::ostringstream fname;
        fname << output_dir << "/" << it->first << "-" << std::setw(3) << std::setfill('0') << step << ".dat";
        std::ofstream file(fname.str().c_str());
        if (!file) {
            THROW("Failed to open " << fname.str());
        }
        const std::vector<double>& d = *(it->second);
        std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
    }
}


static void outputWaterCut(const Opm::Watercut& watercut,
                           const std::string& output_dir)
{
    // Write water cut curve.
    std::string fname = output_dir  + "/watercut.txt";
    std::ofstream os(fname.c_str());
    if (!os) {
        THROW("Failed to open " << fname);
    }
    watercut.write(os);
}


static void outputWellReport(const Opm::WellReport& wellreport,
                             const std::string& output_dir)
{
    // Write well report.
    std::string fname = output_dir  + "/wellreport.txt";
    std::ofstream os(fname.c_str());
    if (!os) {
        THROW("Failed to open " << fname);
    }
    wellreport.write(os);
}


// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    using namespace Opm;

    std::cout << "\n================    Test program for weakly compressible two-phase flow     ===============\n\n";
    Opm::parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Reading various control parameters.
    const bool use_reorder = param.getDefault("use_reorder", true);
    const bool output = param.getDefault("output", true);
    std::string output_dir;
    int output_interval = 1;
    if (output) {
        output_dir = param.getDefault("output_dir", std::string("output"));
        // Ensure that output dir exists
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            THROW("Creating directories failed: " << fpath);
        }
        output_interval = param.getDefault("output_interval", output_interval);
    }
    // const int num_transport_substeps = param.getDefault("num_transport_substeps", 1);

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<Opm::GridManager> grid;
    boost::scoped_ptr<Opm::BlackoilPropertiesInterface> props;
    boost::scoped_ptr<Opm::WellsManager> wells;
    boost::scoped_ptr<Opm::RockCompressibility> rock_comp;
    Opm::SimulatorTimer simtimer;
    Opm::BlackoilState state;
    bool check_well_controls = false;
    int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    if (use_deck) {
        std::string deck_filename = param.get<std::string>("deck_filename");
        Opm::EclipseGridParser deck(deck_filename);
        // Grid init
        grid.reset(new Opm::GridManager(deck));
        // Rock and fluid init
        const int* gc = grid->c_grid()->global_cell;
        std::vector<int> global_cell(gc, gc + grid->c_grid()->number_of_cells);
        props.reset(new Opm::BlackoilPropertiesFromDeck(deck, global_cell));
        // Wells init.
        wells.reset(new Opm::WellsManager(deck, *grid->c_grid(), props->permeability()));
        check_well_controls = param.getDefault("check_well_controls", false);
        max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
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
        props.reset(new Opm::BlackoilPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
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
        if (props->surfaceDensity()[0] == props->surfaceDensity()[1]) {
            std::cout << "**** Warning: nonzero gravity, but zero density difference." << std::endl;
        }
    }
    bool use_segregation_split = false;
    bool use_column_solver = false;
    bool use_gauss_seidel_gravity = false;
    if (use_gravity && use_reorder) {
        use_segregation_split = param.getDefault("use_segregation_split", use_segregation_split);
        if (use_segregation_split) {
            use_column_solver = param.getDefault("use_column_solver", use_column_solver);
            if (use_column_solver) {
                use_gauss_seidel_gravity = param.getDefault("use_gauss_seidel_gravity", use_gauss_seidel_gravity);
            }
        }
    }

    // Check that rock compressibility is not used with solvers that do not handle it.
    // int nl_pressure_maxiter = 0;
    // double nl_pressure_tolerance = 0.0;
    if (rock_comp->isActive()) {
        THROW("No rock compressibility in comp. pressure solver yet.");
        if (!use_reorder) {
            THROW("Cannot run implicit (non-reordering) transport solver with rock compressibility yet.");
        }
        // nl_pressure_maxiter = param.getDefault("nl_pressure_maxiter", 10);
        // nl_pressure_tolerance = param.getDefault("nl_pressure_tolerance", 1.0); // in Pascal
    }

    // Source-related variables init.
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> totmob;
    std::vector<double> omega; // Will remain empty if no gravity.
    std::vector<double> rc; // Will remain empty if no rock compressibility.

    // Extra rock init.
    std::vector<double> porevol;
    if (rock_comp->isActive()) {
        THROW("CompressibleTpfa solver does not handle this.");
        computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, state.pressure(), porevol);
    } else {
        computePorevolume(*grid->c_grid(), props->porosity(), porevol);
    }
    double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(num_cells);
    std::vector<double> src(num_cells, 0.0);

    // Initialising src
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

    std::vector<double> reorder_src = src;

    // Solvers init.
    // Linear solver.
    Opm::LinearSolverFactory linsolver(param);
    // Pressure solver.
    const double nl_press_res_tol = param.getDefault("nl_press_res_tol", 1e-6);
    const double nl_press_change_tol = param.getDefault("nl_press_change_tol", 10.0);
    const int nl_press_maxiter = param.getDefault("nl_press_maxiter", 20);
    const double *grav = use_gravity ? &gravity[0] : 0;
    Opm::CompressibleTpfa psolver(*grid->c_grid(), *props, linsolver,
                                  nl_press_res_tol, nl_press_change_tol, nl_press_maxiter,
                                  grav, wells->c_wells());
    // Reordering solver.
#if TRANSPORT_SOLVER_FIXED
    const double nl_tolerance = param.getDefault("nl_tolerance", 1e-9);
    const int nl_maxiter = param.getDefault("nl_maxiter", 30);
    Opm::TransportModelTwophase reorder_model(*grid->c_grid(), *props, nl_tolerance, nl_maxiter);
    if (use_gauss_seidel_gravity) {
        reorder_model.initGravity(grav);
    }
#endif // TRANSPORT_SOLVER_FIXED
    // Column-based gravity segregation solver.
    std::vector<std::vector<int> > columns;
    if (use_column_solver) {
        Opm::extractColumn(*grid->c_grid(), columns);
    }

    // The allcells vector is used in calls to computeTotalMobility()
    // and computeTotalMobilityOmega().
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }

    // Warn if any parameters are unused.
    if (param.anyUnused()) {
        std::cout << "--------------------   Unused parameters:   --------------------\n";
        param.displayUsage();
        std::cout << "----------------------------------------------------------------" << std::endl;
    }

    // Write parameters used for later reference.
    if (output) {
        param.writeParam(output_dir + "/spu_2p.param");
    }

    // Main simulation loop.
    Opm::time::StopWatch pressure_timer;
    double ptime = 0.0;
    Opm::time::StopWatch transport_timer;
    double ttime = 0.0;
    Opm::time::StopWatch total_timer;
    total_timer.start();
    std::cout << "\n\n================    Starting main simulation loop     ===============" << std::endl;
    double init_satvol[2] = { 0.0 };
    double satvol[2] = { 0.0 };
    double injected[2] = { 0.0 };
    double produced[2] = { 0.0 };
    double tot_injected[2] = { 0.0 };
    double tot_produced[2] = { 0.0 };
    Opm::computeSaturatedVol(porevol, state.saturation(), init_satvol);
    std::cout << "\nInitial saturations are    " << init_satvol[0]/tot_porevol_init
              << "    " << init_satvol[1]/tot_porevol_init << std::endl;
    Opm::Watercut watercut;
    watercut.push(0.0, 0.0, 0.0);
    Opm::WellReport wellreport;
    Opm::WellState well_state;
    std::vector<double> fractional_flows;
    std::vector<double> well_resflows_phase;
    int num_wells = 0;
    if (wells->c_wells()) {
        num_wells = wells->c_wells()->number_of_wells;
        well_state.init(wells->c_wells());
        well_resflows_phase.resize((wells->c_wells()->number_of_phases)*(num_wells), 0.0);
        wellreport.push(*props, *wells->c_wells(),
                        state.pressure(), state.surfacevol(), state.saturation(),
                        0.0, well_state.bhp(), well_state.perfRates());
    }
    for (; !simtimer.done(); ++simtimer) {
        // Report timestep and (optionally) write state to disk.
        simtimer.report(std::cout);
        if (output && (simtimer.currentStepNum() % output_interval == 0)) {
            outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir);
        }

        // Solve pressure.
        if (check_well_controls) {
            computeFractionalFlow(*props, allcells, state.pressure(), state.surfacevol(), state.saturation(), fractional_flows);
        }
        if (check_well_controls) {
            wells->applyExplicitReinjectionControls(well_resflows_phase, well_resflows_phase);
        }
        bool well_control_passed = !check_well_controls;
        int well_control_iteration = 0;
        do { // Well control outer loop.
            pressure_timer.start();
            psolver.solve(simtimer.currentStepLength(), state, well_state);
            pressure_timer.stop();
            double pt = pressure_timer.secsSinceStart();
            std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
            ptime += pt;

            if (check_well_controls) {
                Opm::computePhaseFlowRatesPerWell(*wells->c_wells(),
                                                  fractional_flows,
                                                  well_state.perfRates(),
                                                  well_resflows_phase);
                std::cout << "Checking well conditions." << std::endl;
                // For testing we set surface := reservoir
                well_control_passed = wells->conditionsMet(well_state.bhp(), well_resflows_phase, well_resflows_phase);
                ++well_control_iteration;
                if (!well_control_passed && well_control_iteration > max_well_control_iterations) {
                    THROW("Could not satisfy well conditions in " << max_well_control_iterations << " tries.");
                }
                if (!well_control_passed) {
                    std::cout << "Well controls not passed, solving again." << std::endl;
                } else {
                    std::cout << "Well conditions met." << std::endl;
                }
            }
        } while (!well_control_passed);

        // Process transport sources (to include bdy terms and well flows).
        Opm::computeTransportSource(*grid->c_grid(), src, state.faceflux(), 1.0,
                                    wells->c_wells(), well_state.perfRates(), reorder_src);

        // Solve transport.
        transport_timer.start();
#if TRANSPORT_SOLVER_FIXED
        double stepsize = simtimer.currentStepLength();
        if (num_transport_substeps != 1) {
            stepsize /= double(num_transport_substeps);
            std::cout << "Making " << num_transport_substeps << " transport substeps." << std::endl;
        }
        for (int tr_substep = 0; tr_substep < num_transport_substeps; ++tr_substep) {
            Opm::toWaterSat(state.saturation(), reorder_sat);
            reorder_model.solve(&state.faceflux()[0], &porevol[0], &reorder_src[0],
                                stepsize, &reorder_sat[0]);
            Opm::toBothSat(reorder_sat, state.saturation());
            Opm::computeInjectedProduced(*props, state.saturation(), reorder_src, stepsize, injected, produced);
            if (use_segregation_split) {
                reorder_model.solveGravity(columns, &porevol[0], stepsize, reorder_sat);
                Opm::toBothSat(reorder_sat, state.saturation());
            }
        }
#endif // TRANSPORT_SOLVER_FIXED
        transport_timer.stop();
        double tt = transport_timer.secsSinceStart();
        std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
        ttime += tt;

        // Report volume balances.
        Opm::computeSaturatedVol(porevol, state.saturation(), satvol);
        tot_injected[0] += injected[0];
        tot_injected[1] += injected[1];
        tot_produced[0] += produced[0];
        tot_produced[1] += produced[1];
        std::cout.precision(5);
        const int width = 18;
        std::cout << "\nVolume balance report (all numbers relative to total pore volume).\n";
        std::cout << "    Saturated volumes:     "
                  << std::setw(width) << satvol[0]/tot_porevol_init
                  << std::setw(width) << satvol[1]/tot_porevol_init << std::endl;
        std::cout << "    Injected volumes:      "
                  << std::setw(width) << injected[0]/tot_porevol_init
                  << std::setw(width) << injected[1]/tot_porevol_init << std::endl;
        std::cout << "    Produced volumes:      "
                  << std::setw(width) << produced[0]/tot_porevol_init
                  << std::setw(width) << produced[1]/tot_porevol_init << std::endl;
        std::cout << "    Total inj volumes:     "
                  << std::setw(width) << tot_injected[0]/tot_porevol_init
                  << std::setw(width) << tot_injected[1]/tot_porevol_init << std::endl;
        std::cout << "    Total prod volumes:    "
                  << std::setw(width) << tot_produced[0]/tot_porevol_init
                  << std::setw(width) << tot_produced[1]/tot_porevol_init << std::endl;
        std::cout << "    In-place + prod - inj: "
                  << std::setw(width) << (satvol[0] + tot_produced[0] - tot_injected[0])/tot_porevol_init
                  << std::setw(width) << (satvol[1] + tot_produced[1] - tot_injected[1])/tot_porevol_init << std::endl;
        std::cout << "    Init - now - pr + inj: "
                  << std::setw(width) << (init_satvol[0] - satvol[0] - tot_produced[0] + tot_injected[0])/tot_porevol_init
                  << std::setw(width) << (init_satvol[1] - satvol[1] - tot_produced[1] + tot_injected[1])/tot_porevol_init
                  << std::endl;
        std::cout.precision(8);

        watercut.push(simtimer.currentTime() + simtimer.currentStepLength(),
                      produced[0]/(produced[0] + produced[1]),
                      tot_produced[0]/tot_porevol_init);
        if (wells->c_wells()) {
            wellreport.push(*props, *wells->c_wells(),
                            state.pressure(), state.surfacevol(), state.saturation(),
                            simtimer.currentTime() + simtimer.currentStepLength(),
                            well_state.bhp(), well_state.perfRates());
        }
    }
    total_timer.stop();

    std::cout << "\n\n================    End of simulation     ===============\n"
              << "Total time taken: " << total_timer.secsSinceStart()
              << "\n  Pressure time:  " << ptime
              << "\n  Transport time: " << ttime << std::endl;

    if (output) {
        outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir);
        outputWaterCut(watercut, output_dir);
        if (wells->c_wells()) {
            outputWellReport(wellreport, output_dir);
        }
    }
}
