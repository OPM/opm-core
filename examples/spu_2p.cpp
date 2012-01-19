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

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>

#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>

#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>

#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

// template <class Ostream, class Collection>
// Ostream&
// operator<<(Ostream& os, const Collection& c)
// {
//     typedef typename Collection::value_type VT;

//     os << "[ ";
//     std::copy(c.begin(), c.end(), ::std::ostream_iterator<VT>(os, " "));
//     os << "]";

//     return os;
// }

class Rock {
public:
    Rock(::std::size_t nc, ::std::size_t dim)
        : dim_ (dim           ),
          perm_(nc * dim * dim),
          poro_(nc            ) {}

    const ::std::vector<double>& perm() const { return perm_; }
    const ::std::vector<double>& poro() const { return poro_; }

    void
    perm_homogeneous(double k) {
        setVector(0.0, perm_);

        const ::std::size_t d2 = dim_ * dim_;

        for (::std::size_t c = 0, nc = poro_.size(); c < nc; ++c) {
            for (::std::size_t i = 0; i < dim_; ++i) {
                perm_[c*d2 + i*(dim_ + 1)] = k;
            }
        }
    }

    void
    poro_homogeneous(double phi) {
        setVector(phi, poro_);
    }

private:
    void
    setVector(double x, ::std::vector<double>& v) {
        ::std::fill(v.begin(), v.end(), x);
    }

    ::std::size_t         dim_ ;
    ::std::vector<double> perm_;
    ::std::vector<double> poro_;
};

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
    PressureSolver(UnstructuredGrid* g, const Rock& rock)
        : htrans_(g->cell_facepos[ g->number_of_cells ]),
          trans_ (g->number_of_faces),
          gpress_(g->cell_facepos[ g->number_of_cells ])
    {
        tpfa_htrans_compute(g, &rock.perm()[0], &htrans_[0]);

        h_ = ifs_tpfa_construct(g);
    }

    ~PressureSolver() {
        ifs_tpfa_destroy(h_);
    }

    template <class State>
    void
    solve(UnstructuredGrid*                      g     ,
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


typedef Opm::SimpleFluid2p<>                          TwophaseFluid;
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
compute_porevolume(const UnstructuredGrid*        g,
                   const Rock&          rock,
                   std::vector<double>& porevol)
{
    const ::std::vector<double>& poro = rock.poro();

#if 0
    assert (poro.size() == static_cast<::std::size_t>(g->number_of_cells));
#endif

    porevol.resize(rock.poro().size());

    ::std::transform(poro.begin(), poro.end(),
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



// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    Opm::parameter::ParameterGroup param(argc, argv, false);
    const int nx = param.getDefault("nx", 100);
    const int ny = param.getDefault("ny", 100);
    const int nz = param.getDefault("nz", 1);
    const int num_psteps = param.getDefault("num_psteps", 1);
    double stepsize_days = param.getDefault("stepsize_days", 1.0);
    double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);

    std::tr1::array<int, 3> grid_dims = {{ nx, ny, nz }};
    std::tr1::array<double, 3> cell_size = {{ 1.0, 1.0, 1.0 }};

    UnstructuredGrid* grid = create_cart_grid(nx, ny, nz);
    Rock rock(grid->number_of_cells, grid->dimensions);
    rock.perm_homogeneous(1);
    rock.poro_homogeneous(1);

    PressureSolver psolver(grid, rock);

    std::vector<double> totmob(grid->number_of_cells, 1.0);
    std::vector<double> src   (grid->number_of_cells, 0.0);

    src[0]                         =  1.0;
    src[grid->number_of_cells - 1] = -1.0;

    ReservoirState state(grid);

    TransportSource* tsrc = create_transport_source(2, 2);
    double ssrc[]   = { 1.0, 0.0 };
    double ssink[]  = { 0.0, 1.0 };
    double zdummy[] = { 0.0, 0.0 };
    append_transport_source(0, 2, 0, src[0], ssrc, zdummy, tsrc);
    append_transport_source(grid->number_of_cells - 1, 2, 0,
			    src.back(), ssink, zdummy, tsrc);

    std::tr1::array<double, 2> mu  = {{ 1.0, 1.0 }};
    std::tr1::array<double, 2> rho = {{ 0.0, 0.0 }};
    TwophaseFluid fluid(mu, rho);

    std::vector<double> porevol;
    compute_porevolume(grid, rock, porevol);

    TransportModel  model  (fluid, *grid, porevol, 0, false);
    TransportSolver tsolver(model);

    Opm::ImplicitTransportDetails::NRReport  rpt;
    Opm::ImplicitTransportDetails::NRControl ctrl;
    double current_time = 0.0;
    double total_time = stepsize*num_psteps;
    ctrl.max_it = 20;

    using Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;
    CSRMatrixUmfpackSolver linsolve;

    for (int pstep = 0; pstep < num_psteps; ++pstep) {
        std::cout << "\n\n================    Simulation step number " << pstep
                  << "    ==============="
                  << "\n      Current time (days)     " << Opm::unit::convert::to(current_time, Opm::unit::day)
                  << "\n      Current stepsize (days) " << Opm::unit::convert::to(stepsize, Opm::unit::day)
                  << "\n      Total time (days)       " << Opm::unit::convert::to(total_time, Opm::unit::day)
                  << "\n" << std::endl;

	outputState(grid_dims, cell_size, state, pstep);

	psolver.solve(grid, totmob, src, state);

	tsolver.solve(*grid, tsrc, stepsize, ctrl, state, linsolve, rpt);

	std::cout << rpt;

	current_time += stepsize;
    }

    outputState(grid_dims, cell_size, state, num_psteps);

    destroy_transport_source(tsrc);
    destroy_cart_grid(grid);
}
