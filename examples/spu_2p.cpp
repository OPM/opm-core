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
#include "config.h"

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

#include <opm/core/linalg/sparse_sys.h>

#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/utility/cart_grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/SimpleFluid2p.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

#include <opm/core/transport/reorder/twophasetransport.hpp>




class ReservoirState {
public:
    ReservoirState(const UnstructuredGrid* g, const int num_phases = 2)
        : press_ (g->number_of_cells),
          fpress_(g->number_of_faces),
          flux_  (g->number_of_faces),
          sat_   (num_phases * g->number_of_cells)
    {}

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




class PressureSolver {
public:
    PressureSolver(UnstructuredGrid* g,
		   const Opm::IncompPropertiesInterface& props)
        : htrans_(g->cell_facepos[ g->number_of_cells ]),
          trans_ (g->number_of_faces),
          gpress_(g->cell_facepos[ g->number_of_cells ])
    {
        tpfa_htrans_compute(g, props.permeability(), &htrans_[0]);

        h_ = ifs_tpfa_construct(g);
    }

    ~PressureSolver() {
        ifs_tpfa_destroy(h_);
    }

    template <class State>
    void
    solve(UnstructuredGrid*            g     ,
          const ::std::vector<double>& totmob,
          const ::std::vector<double>& src   ,
          State&                       state ) {

        tpfa_eff_trans_compute(g, &totmob[0], &htrans_[0], &trans_[0]);

        // No gravity
        ::std::fill(gpress_.begin(), gpress_.end(), double(0.0));

        ifs_tpfa_assemble(g, &trans_[0], &src[0], &gpress_[0], h_);

        using Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;

        CSRMatrixUmfpackSolver linsolve;
        linsolve.solve(h_->A, h_->b, h_->x);

        ifs_tpfa_press_flux(g, &trans_[0], h_,
                            &state.pressure()[0],
                            &state.faceflux()[0]);
    }

private:
    ::std::vector<double> htrans_;
    ::std::vector<double> trans_ ;
    ::std::vector<double> gpress_;

    struct ifs_tpfa_data* h_;
};




typedef Opm::SimpleFluid2p<2>                         TwophaseFluid;
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




static void
compute_porevolume(const UnstructuredGrid* g,
                   const Opm::IncompPropertiesInterface& props,
                   std::vector<double>& porevol)
{
    int num_cells = g->number_of_cells;
    porevol.resize(num_cells);
    const double* poro = props.porosity();
    ::std::transform(poro, poro + num_cells,
                     g->cell_volumes,
                     porevol.begin(),
                     ::std::multiplies<double>());
}




template <class State>
void outputState(const std::tr1::array<int, 3>& grid_dims,
		 const std::tr1::array<double, 3>& cell_size,
		 const State& state,
		 const int step)
{
    std::ostringstream vtkfilename;
    vtkfilename << "output-" << std::setw(3) << std::setfill('0') << step << ".vtk";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
	THROW("Failed to open " << vtkfilename.str());
    }
    writeVtkDataAllCartesian(grid_dims, cell_size, state, vtkfile);
}




template <class State>
void writeVtkDataAllCartesian(const std::tr1::array<int, 3>& dims,
			      const std::tr1::array<double, 3>& cell_size,
			      const State& state,
			      std::ostream& vtk_file)
{
    // Dimension is hardcoded in the prototype and the next two lines,
    // but the rest is flexible (allows dimension == 2 or 3).
    int dimension = 3;
    int num_cells = dims[0]*dims[1]*dims[2];

    ASSERT(dimension == 2 || dimension == 3);
    ASSERT(num_cells = dims[0]*dims[1]* (dimension == 2 ? 1 : dims[2]));

    vtk_file << "# vtk DataFile Version 2.0\n";
    vtk_file << "Structured Grid\n \n";
    vtk_file << "ASCII \n";
    vtk_file << "DATASET STRUCTURED_POINTS\n";

    vtk_file << "DIMENSIONS "
	     << dims[0] + 1 << " "
	     << dims[1] + 1 << " ";
    if (dimension == 3) {
	vtk_file << dims[2] + 1;
    } else {
	vtk_file << 1;
    }
    vtk_file << "\n";
	
    vtk_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";

    vtk_file << "SPACING " << cell_size[0] << " " << cell_size[1];
    if (dimension == 3) {
	vtk_file << " " << cell_size[2];
    } else {
	vtk_file << " " << 0.0;
    }
    vtk_file << "\n";

    vtk_file << "CELL_DATA " << num_cells << '\n';
    vtk_file << "SCALARS pressure float" << '\n';
    vtk_file << "LOOKUP_TABLE pressure_table " << '\n';
    for (int i = 0; i < num_cells; ++i) {
	vtk_file << state.pressure()[i] << '\n';
    }

    ASSERT(state.numPhases() == 2);
    vtk_file << "SCALARS saturation float" << '\n';
    vtk_file << "LOOKUP_TABLE saturation_table " << '\n';
    for (int i = 0; i < num_cells; ++i) {
	double s = state.saturation()[2*i];
    	if (s > 1e-10) {
    	    vtk_file << s << '\n';
    	} else {
    	    vtk_file << 0.0 << '\n';
    	}
    }
}




static void toWaterSat(const std::vector<double>& sboth, std::vector<double>& sw)
{
    int num = sboth.size()/2;
    sw.resize(num);
    for (int i = 0; i < num; ++i) {
	sw[i] = sboth[2*i];
    }
}

static void toBothSat(const std::vector<double>& sw, std::vector<double>& sboth)
{
    int num = sw.size();
    sboth.resize(2*num);
    for (int i = 0; i < num; ++i) {
	sboth[2*i] = sw[i];
	sboth[2*i + 1] = 1.0 - sw[i];
    }
}




// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    Opm::parameter::ParameterGroup param(argc, argv, false);
    const int nx = param.getDefault("nx", 100);
    const int ny = param.getDefault("ny", 100);
    const int nz = param.getDefault("nz", 1);
    const int num_psteps = param.getDefault("num_psteps", 1);
    const double stepsize_days = param.getDefault("stepsize_days", 1.0);
    const double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);
    const bool guess_old_solution = param.getDefault("guess_old_solution", false);
    const bool use_reorder = param.getDefault("use_reorder", false);
    const bool output = param.getDefault("output", true);

    // Grid init.
    std::tr1::array<int, 3> grid_dims = {{ nx, ny, nz }};
    std::tr1::array<double, 3> cell_size = {{ 1.0, 1.0, 1.0 }};
    UnstructuredGrid* grid = create_cart_grid(nx, ny, nz);

    // Rock init.
    Opm::IncompPropertiesBasic props(param, grid->dimensions, grid->number_of_cells);
    std::vector<double> porevol;
    compute_porevolume(grid, props, porevol);
    double tot_porevol = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Fluid init.
    std::tr1::array<double, 2> mu  = {{ 0.001, 0.003 }};
    std::tr1::array<double, 2> rho = {{ 0.0, 0.0 }};
    TwophaseFluid fluid(mu, rho);

    // Solvers init.
    PressureSolver psolver(grid, props);
    TransportModel  model  (fluid, *grid, porevol, 0, guess_old_solution);
    TransportSolver tsolver(model);

    // State-related and source-related variables init.
    std::vector<double> totmob(grid->number_of_cells, 1.0);
    ReservoirState state(grid);
    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(grid->number_of_cells);
    double flow_per_sec = 0.1*tot_porevol/Opm::unit::day;
    std::vector<double> src   (grid->number_of_cells, 0.0);
    src[0]                         =  flow_per_sec;
    src[grid->number_of_cells - 1] = -flow_per_sec;
    TransportSource* tsrc = create_transport_source(2, 2);
    double ssrc[]   = { 1.0, 0.0 };
    double ssink[]  = { 0.0, 1.0 };
    double zdummy[] = { 0.0, 0.0 };
    append_transport_source(0, 2, 0, src[0], ssrc, zdummy, tsrc);
    append_transport_source(grid->number_of_cells - 1, 2, 0,
			    src.back(), ssink, zdummy, tsrc);
    std::vector<double> reorder_src = src;

    // Control init.
    Opm::ImplicitTransportDetails::NRReport  rpt;
    Opm::ImplicitTransportDetails::NRControl ctrl;
    double current_time = 0.0;
    double total_time = stepsize*num_psteps;
    ctrl.max_it = param.getDefault("max_it", 20);
    ctrl.verbosity = param.getDefault("verbosity", 0);
    ctrl.max_it_ls = param.getDefault("max_it_ls", 5);

    // Linear solver init.
    using Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;
    CSRMatrixUmfpackSolver linsolve;

    // Main simulation loop.
    for (int pstep = 0; pstep < num_psteps; ++pstep) {
        std::cout << "\n\n================    Simulation step number " << pstep
                  << "    ==============="
                  << "\n      Current time (days)     " << Opm::unit::convert::to(current_time, Opm::unit::day)
                  << "\n      Current stepsize (days) " << Opm::unit::convert::to(stepsize, Opm::unit::day)
                  << "\n      Total time (days)       " << Opm::unit::convert::to(total_time, Opm::unit::day)
                  << "\n" << std::endl;

	if (output) {
	    outputState(grid_dims, cell_size, state, pstep);
	}

	psolver.solve(grid, totmob, src, state);

	if (use_reorder) {
	    toWaterSat(state.saturation(), reorder_sat);
	    // We must treat reorder_src here,
	    // if we are to handle anything but simple water
	    // injection, since it is expected to be
	    // equal to total outflow (if negative)
	    // and water inflow (if positive).
	    // Also, for anything but noflow boundaries,
	    // boundary flows must be accumulated into
	    // source term following the same convention.
	    twophasetransport(&porevol[0],
			      &reorder_src[0],
			      stepsize,
			      grid,
			      &state.faceflux()[0],
			      NULL,
			      &reorder_sat[0]);
	    toBothSat(reorder_sat, state.saturation());
	} else {
	    tsolver.solve(*grid, tsrc, stepsize, ctrl, state, linsolve, rpt);
	    std::cout << rpt;
	}

	current_time += stepsize;
    }

    if (output) {
	outputState(grid_dims, cell_size, state, num_psteps);
    }

    destroy_transport_source(tsrc);
    destroy_cart_grid(grid);
}
