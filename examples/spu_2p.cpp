/*===========================================================================
//
// File: spu_2p.cpp
//
// Created: 2011-10-05 10:29:01+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>

// #define EXPERIMENT_ISTL
#ifdef EXPERIMENT_ISTL
#include <opm/core/linalg/LinearSolverIstl.hpp>
#endif

#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/GridManager.hpp>
#include <opm/core/grid.h>
#include <opm/core/WellsManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/SimpleFluid2p.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

#include <opm/core/ColumnExtract.hpp>
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




class ReservoirState {
public:
    ReservoirState(const UnstructuredGrid* g, const double init_sat = 0.0)
    : press_ (g->number_of_cells, 0.0),
      fpress_(g->number_of_faces, 0.0),
      flux_  (g->number_of_faces, 0.0),
      sat_   (2 * g->number_of_cells, 0.0)
    {
	for (int cell = 0; cell < g->number_of_cells; ++cell) {
	    sat_[2*cell] = init_sat;
	    sat_[2*cell + 1] = 1.0 - init_sat;
	}
    }

    enum ExtremalSat { MinSat, MaxSat };

    void setToMinimumWaterSat(const Opm::IncompPropertiesInterface& props)
    {
	const int n = props.numCells();
	std::vector<int> cells(n);
	for (int i = 0; i < n; ++i) {
	    cells[i] = i;
	}
	setWaterSat(cells, props, MinSat);
    }

    void setWaterSat(const std::vector<int>& cells,
		     const Opm::IncompPropertiesInterface& props,
		     ExtremalSat es)
    {
	const int n = cells.size();
	std::vector<double> smin(2*n);
	std::vector<double> smax(2*n);
	props.satRange(n, &cells[0], &smin[0], &smax[0]);
	const double* svals = (es == MinSat) ? &smin[0] : &smax[0];
	for (int ci = 0; ci < n; ++ci) {
	    const int cell = cells[ci];
	    sat_[2*cell] = svals[2*ci];
	    sat_[2*cell + 1] = 1.0 - sat_[2*cell];
	}
    }

    int numPhases() const { return sat_.size()/press_.size(); }

    ::std::vector<double>& pressure    () { return press_ ; }
    ::std::vector<double>& facepressure() { return fpress_; }
    ::std::vector<double>& faceflux    () { return flux_  ; }
    ::std::vector<double>& saturation  () { return sat_   ; }

    const ::std::vector<double>& pressure    () const { return press_ ; }
    const ::std::vector<double>& facepressure() const { return fpress_; }
    const ::std::vector<double>& faceflux    () const { return flux_  ; }
    const ::std::vector<double>& saturation  () const { return sat_   ; }

private:
    ::std::vector<double> press_ ;
    ::std::vector<double> fpress_;
    ::std::vector<double> flux_  ;
    ::std::vector<double> sat_   ;
};



template <class State>
void outputState(const UnstructuredGrid* grid,
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
    Opm::estimateCellVelocity(*grid, state.faceflux(), cell_velocity);
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

/// Create a src vector equivalent to a wells structure.
/// For this to be valid, the wells must be all rate-controlled and
/// single-perforation.
void wellsToSrc(const Wells& wells, const int num_cells, std::vector<double>& src)
{
    src.resize(num_cells);
    for (int w = 0; w < wells.number_of_wells; ++w) {
	if (wells.ctrls[w]->num != 1) {
	    THROW("In wellsToSrc(): well has more than one control.");
	}
	if (wells.ctrls[w]->type[0] != RATE) {
	    THROW("In wellsToSrc(): well is BHP, not RATE.");
	}
	if (wells.well_connpos[w+1] - wells.well_connpos[w] != 1) {
	    THROW("In wellsToSrc(): well has multiple perforations.");
	}
	const double flow = wells.ctrls[w]->target[0];
	const double cell = wells.well_cells[wells.well_connpos[w]];
	src[cell] = (wells.type[w] == INJECTOR) ? flow : -flow;
    }
}


// --------------- Types needed to define transport solver ---------------

class SimpleFluid2pWrappingProps
{
public:
    SimpleFluid2pWrappingProps(const Opm::IncompPropertiesInterface& props)
    : props_(props),
      smin_(props.numCells()*props.numPhases()),
      smax_(props.numCells()*props.numPhases())
    {
        if (props.numPhases() != 2) {
            THROW("SimpleFluid2pWrapper requires 2 phases.");
        }
        const int num_cells = props.numCells();
        std::vector<int> cells(num_cells);
        for (int c = 0; c < num_cells; ++c) {
            cells[c] = c;
        }
        props.satRange(num_cells, &cells[0], &smin_[0], &smax_[0]);
    }

    double density(int phase) const
    {
        return props_.density()[phase];
    }

    template <class Sat,
    class Mob,
    class DMob>
    void mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const
    {
        props_.relperm(1, &s[0], &c, &mob[0], &dmob[0]);
        const double* mu = props_.viscosity();
        mob[0] /= mu[0];
        mob[1] /= mu[1];
        // Recall that we use Fortran ordering for kr derivatives,
        // therefore dmob[i*2 + j] is row j and column i of the
        // matrix.
        // Each row corresponds to a kr function, so which mu to
        // divide by also depends on the row, j.
        dmob[0*2 + 0] /= mu[0];
        dmob[0*2 + 1] /= mu[1];
        dmob[1*2 + 0] /= mu[0];
        dmob[1*2 + 1] /= mu[1];
    }

    template <class Sat,
    class Pcap,
    class DPcap>
    void pc(int c, const Sat& s, Pcap& pcap, DPcap& dpcap) const
    {
        double pc[2];
        double dpc[4];
        props_.capPress(1, &s[0], &c, pc, dpc);
        pcap = pc[0];
        ASSERT(pc[1] == 0.0);
        dpcap = dpc[0];
        ASSERT(dpc[1] == 0.0);
        ASSERT(dpc[2] == 0.0);
        ASSERT(dpc[3] == 0.0);
    }

    double s_min(int c) const
    {
        return smin_[2*c + 0];
    }

    double s_max(int c) const
    {
        return smax_[2*c + 0];
    }

private:
    const Opm::IncompPropertiesInterface& props_;
    std::vector<double> smin_;
    std::vector<double> smax_;
};

typedef SimpleFluid2pWrappingProps TwophaseFluid;
typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;

using namespace Opm::ImplicitTransportDefault;

typedef NewtonVectorCollection< ::std::vector<double> >      NVecColl;
typedef JacobianSystem        < struct CSRMatrix, NVecColl > JacSys;

template <class Vector>
class MaxNorm {
public:
    static double
    norm(const Vector& v) {
        return AccumulationNorm <Vector, MaxAbs>::norm(v);
    }
};

typedef Opm::ImplicitTransport<TransportModel,
        JacSys        ,
        MaxNorm       ,
        VectorNegater ,
        VectorZero    ,
        MatrixZero    ,
        VectorAssign  > TransportSolver;



// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    std::cout << "\n================    Test program for incompressible two-phase flow     ===============\n\n";
    Opm::parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Reading various control parameters.
    const int num_psteps = param.getDefault("num_psteps", 1);
    const double stepsize_days = param.getDefault("stepsize_days", 1.0);
    const double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);
    const bool guess_old_solution = param.getDefault("guess_old_solution", false);
    const bool use_reorder = param.getDefault("use_reorder", true);
    const bool output = param.getDefault("output", true);
    std::string output_dir;
    if (output) {
        output_dir = param.getDefault("output_dir", std::string("output"));
        // Ensure that output dir exists
        boost::filesystem::path fpath(output_dir);
        create_directories(fpath);
    }

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<Opm::GridManager> grid;
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
    boost::scoped_ptr<Opm::WellsManager> wells;
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
    }

    // Extra rock init.
    std::vector<double> porevol;
    computePorevolume(*grid->c_grid(), *props, porevol);
    double tot_porevol = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Extra fluid init for transport solver.
    TwophaseFluid fluid(*props);

    // Gravity init.
    double gravity[3] = { 0.0 };
    double g = param.getDefault("gravity", 0.0);
    bool use_gravity = g != 0.0;
    if (use_gravity) {
        gravity[grid->c_grid()->dimensions - 1] = g;
        if (props->density()[0] == props->density()[1]) {
            std::cout << "**** Warning: nonzero gravity, but zero density difference." << std::endl;
        }
    }
    bool use_segregation_split = false;
    bool use_column_solver = false;
    if (use_gravity && use_reorder) {
        use_segregation_split = param.getDefault("use_segregation_split", use_segregation_split);
        if (use_segregation_split) {
            use_column_solver = param.getDefault("use_column_solver", use_column_solver);
        }
    }

    // Solvers init.
    // Pressure solver.
#ifdef EXPERIMENT_ISTL
    Opm::LinearSolverIstl linsolver(param);
#else
    Opm::LinearSolverUmfpack linsolver;
#endif // EXPERIMENT_ISTL
    const double *grav = use_gravity ? &gravity[0] : 0;
    Opm::IncompTpfa psolver(*grid->c_grid(), props->permeability(), grav, linsolver);
    // Non-reordering solver.
    TransportModel  model  (fluid, *grid->c_grid(), porevol, grav, guess_old_solution);
    if (use_gravity) {
        model.initGravityTrans(*grid->c_grid(), psolver.getHalfTrans());
    }
    TransportSolver tsolver(model);
    // Reordering solver.
    const double nltol = param.getDefault("nl_tolerance", 1e-9);
    const int maxit = param.getDefault("nl_maxiter", 30);
    Opm::TransportModelTwophase reorder_model(*grid->c_grid(), &porevol[0], *props, nltol, maxit);
    // Column-based gravity segregation solver.
    typedef std::map<int, std::vector<int> > ColMap;
    ColMap columns;
    if (use_column_solver) {
        Opm::extractColumn(*grid->c_grid(), columns);
    }
    Opm::GravityColumnSolver<TransportModel> colsolver(model, *grid->c_grid(), nltol, maxit);

    // Boundary conditions.
    Opm::FlowBCManager bcs;

    // State-related and source-related variables init.
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> totmob;
    std::vector<double> omega; // Will remain empty if no gravity.
    double init_sat = param.getDefault("init_sat", 0.0);
    ReservoirState state(grid->c_grid(), init_sat);
    if (!param.has("init_sat")) {
	state.setToMinimumWaterSat(*props);
    }
    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(num_cells);
    std::vector<double> src(num_cells, 0.0);
    int scenario = param.getDefault("scenario", 0);
    switch (scenario) {
    case 0:
	{
	    std::cout << "==== Scenario 0: single-cell source and sink.\n";
	    if (wells->c_wells()) {
		wellsToSrc(*wells->c_wells(), num_cells, src);
	    } else {
		double flow_per_sec = 0.1*tot_porevol/Opm::unit::day;
		src[0] = flow_per_sec;
		src[grid->c_grid()->number_of_cells - 1] = -flow_per_sec;
	    }
	    break;
	}
    case 1:
	{
	    std::cout << "==== Scenario 1: half source, half sink.\n";
	    double flow_per_sec = 0.1*porevol[0]/Opm::unit::day;
	    std::fill(src.begin(), src.begin() + src.size()/2, flow_per_sec);
	    std::fill(src.begin() + src.size()/2, src.end(), -flow_per_sec);
	    break;
	}
    case 2:
	{
	    std::cout << "==== Scenario 2: gravity convection.\n";
	    if (!use_gravity) {
		std::cout << "**** Warning: running gravity convection scenario, but gravity is zero." << std::endl;
	    }
	    if (use_deck) {
		std::cout << "**** Warning: running gravity convection scenario, which expects a cartesian grid."
			  << std::endl;
	    }
	    std::vector<int> left_cells;
	    left_cells.reserve(num_cells/2);
	    const int *glob_cell = grid->c_grid()->global_cell;
	    for (int cell = 0; cell < num_cells; ++cell) {
		const int* cd = grid->c_grid()->cartdims;
		const int gc = glob_cell == 0 ? cell : glob_cell[cell];
		bool left = (gc % cd[0]) < cd[0]/2;
		if (left) {
		    left_cells.push_back(cell);
		}
	    }
	    state.setWaterSat(left_cells, *props, ReservoirState::MaxSat);
	    break;
	}
    case 3:
	{
	    std::cout << "==== Scenario 3: gravity segregation.\n";
	    if (!use_gravity) {
		std::cout << "**** Warning: running gravity segregation scenario, but gravity is zero." << std::endl;
	    }
	    if (use_deck) {
		std::cout << "**** Warning: running gravity segregation scenario, which expects a cartesian grid."
			  << std::endl;
	    }
	    std::vector<double>& sat = state.saturation();
	    const int *glob_cell = grid->c_grid()->global_cell;
	    // Water on top
	    for (int cell = 0; cell < num_cells; ++cell) {
		const int* cd = grid->c_grid()->cartdims;
		const int gc = glob_cell == 0 ? cell : glob_cell[cell];
		bool top = (gc / cd[0] / cd[1]) < cd[2]/2;
		sat[2*cell] = top ? 1.0 : 0.0;
		sat[2*cell + 1 ] = 1.0 - sat[2*cell];
	    }
	    break;
	}
    default:
	{
	    THROW("==== Scenario " << scenario << " is unknown.");
	}
    }
    TransportSource* tsrc = create_transport_source(2, 2);
    double ssrc[]   = { 1.0, 0.0 };
    double ssink[]  = { 0.0, 1.0 };
    double zdummy[] = { 0.0, 0.0 };
    for (int cell = 0; cell < num_cells; ++cell) {
        if (src[cell] > 0.0) {
            append_transport_source(cell, 2, 0, src[cell], ssrc, zdummy, tsrc);
        } else if (src[cell] < 0.0) {
            append_transport_source(cell, 2, 0, src[cell], ssink, zdummy, tsrc);
        }
    }
    std::vector<double> reorder_src = src;

    // Dirichlet boundary conditions.
    if (param.getDefault("use_pside", false)) {
	int pside = param.get<int>("pside");
	double pside_pressure = param.get<double>("pside_pressure");
	bcs.pressureSide(*grid->c_grid(), Opm::FlowBCManager::Side(pside), pside_pressure);
    }

    // Control init.
    Opm::ImplicitTransportDetails::NRReport  rpt;
    Opm::ImplicitTransportDetails::NRControl ctrl;
    double current_time = 0.0;
    double total_time = stepsize*num_psteps;
    if (!use_reorder || use_segregation_split) {
        ctrl.max_it = param.getDefault("max_it", 20);
        ctrl.verbosity = param.getDefault("verbosity", 0);
        ctrl.max_it_ls = param.getDefault("max_it_ls", 5);
    }

    // Linear solver init.
    using Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;
    CSRMatrixUmfpackSolver linsolve;

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
    double aver_sats[2] = { 0.0 };
    Opm::computeAverageSat(porevol, state.saturation(), aver_sats);
    std::cout << "\nInitial average saturations are    " << aver_sats[0] << "    " << aver_sats[1] << std::endl;
    for (int pstep = 0; pstep < num_psteps; ++pstep) {
        std::cout << "\n\n---------------    Simulation step number " << pstep
                << "    ---------------"
                << "\n      Current time (days)     " << Opm::unit::convert::to(current_time, Opm::unit::day)
        << "\n      Current stepsize (days) " << Opm::unit::convert::to(stepsize, Opm::unit::day)
        << "\n      Total time (days)       " << Opm::unit::convert::to(total_time, Opm::unit::day)
        << "\n" << std::endl;

        if (output) {
            outputState(grid->c_grid(), state, pstep, output_dir);
        }

        // Solve pressure.
        if (use_gravity) {
            computeTotalMobilityOmega(*props, allcells, state.saturation(), totmob, omega);
        } else {
            computeTotalMobility(*props, allcells, state.saturation(), totmob);
        }
        pressure_timer.start();
        psolver.solve(totmob, omega, src, bcs.c_bcs(), state.pressure(), state.faceflux());
        pressure_timer.stop();
        double pt = pressure_timer.secsSinceStart();
        std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
        ptime += pt;

	// Process transport sources (to include bdy terms).
	if (use_reorder) {
	    Opm::computeTransportSource(*grid->c_grid(), src, state.faceflux(), 1.0, reorder_src);
	} else {
	    clear_transport_source(tsrc);
	    for (int cell = 0; cell < num_cells; ++cell) {
		if (src[cell] > 0.0) {
		    append_transport_source(cell, 2, 0, src[cell], ssrc, zdummy, tsrc);
		} else if (src[cell] < 0.0) {
		    append_transport_source(cell, 2, 0, src[cell], ssink, zdummy, tsrc);
		}
	    }
	}

        // Solve transport
        transport_timer.start();
        if (use_reorder) {
            Opm::toWaterSat(state.saturation(), reorder_sat);
            reorder_model.solve(&state.faceflux()[0], &reorder_src[0], stepsize, &reorder_sat[0]);
            Opm::toBothSat(reorder_sat, state.saturation());
            if (use_segregation_split) {
                if (use_column_solver) {
                    colsolver.solve(columns, stepsize, state.saturation());
                } else {
                    std::vector<double> fluxes = state.faceflux();
                    std::fill(state.faceflux().begin(), state.faceflux().end(), 0.0);
                    tsolver.solve(*grid->c_grid(), tsrc, stepsize, ctrl, state, linsolve, rpt);
                    std::cout << rpt;
                    state.faceflux() = fluxes;
                }
            }
        } else {
            tsolver.solve(*grid->c_grid(), tsrc, stepsize, ctrl, state, linsolve, rpt);
            std::cout << rpt;
        }
        transport_timer.stop();
        double tt = transport_timer.secsSinceStart();
        std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
        ttime += tt;

	Opm::computeAverageSat(porevol, state.saturation(), aver_sats);
	std::cout << "\nAverage saturations are    " << aver_sats[0] << "    " << aver_sats[1] << std::endl;

        current_time += stepsize;
    }
    total_timer.stop();

    std::cout << "\n\n================    End of simulation     ===============\n"
            << "Total time taken: " << total_timer.secsSinceStart()
            << "\n  Pressure time:  " << ptime
            << "\n  Transport time: " << ttime << std::endl;

    if (output) {
        outputState(grid->c_grid(), state, num_psteps, output_dir);
    }

    destroy_transport_source(tsrc);
}
