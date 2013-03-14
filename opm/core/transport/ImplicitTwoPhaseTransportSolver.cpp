/*===========================================================================
//
// File: ImpliciteTwoPhaseTransportSolver.cpp
//
// Author: hnil <hnil@sintef.no>
//
// Created: 9 Nov 2012
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


#include <opm/core/transport/ImplicitTwoPhaseTransportSolver.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/utility/miscUtilities.hpp>
namespace Opm{

    ImplicitTwoPhaseTransportSolver::ImplicitTwoPhaseTransportSolver(
            const Opm::WellsManager& wells,
            const Opm::RockCompressibility& rock_comp,
            const ImplicitTransportDetails::NRControl& ctrl,
            SinglePointUpwindTwoPhase<Opm::SimpleFluid2pWrappingProps>& model,
            const UnstructuredGrid& grid,
            const Opm::IncompPropertiesInterface& props,
            const parameter::ParameterGroup& param)
        : tsolver_(model),
          grid_(grid),
          ctrl_(ctrl),
          props_(props),
          rock_comp_(rock_comp),
          wells_(wells)
    {
        tsrc_ = create_transport_source(2, 2);
        //linsolver_(param),
        //src_(num_cells, 0.0);
    }

    void ImplicitTwoPhaseTransportSolver::solve(const double* porevolume,
                                                 const double* source,
                                                 const double dt,
                                                 TwophaseState& state,
                                                 const WellState& well_state)
    {
        std::vector<double> porevol;
        if (rock_comp_.isActive()) {
            computePorevolume(grid_, props_.porosity(), rock_comp_, state.pressure(), porevol);
        }
        std::vector<double> src(grid_.number_of_cells, 0.0);
        //Opm::wellsToSrc(*wells->c_wells(), num_cells, src);
        Opm::computeTransportSource(grid_, src, state.faceflux(), 1.0,
                                    wells_.c_wells(), well_state.perfRates(), src);
        double ssrc[]   = { 1.0, 0.0 };
        double ssink[]  = { 0.0, 1.0 };
        double zdummy[] = { 0.0, 0.0 };
        for (int cell = 0; cell < grid_.number_of_cells; ++cell) {
            clear_transport_source(tsrc_);
            if (src[cell] > 0.0) {
                append_transport_source(cell, 2, 0, src[cell], ssrc, zdummy, tsrc_);
            } else if (src[cell] < 0.0) {
                append_transport_source(cell, 2, 0, src[cell], ssink, zdummy, tsrc_);
            }
        }
        // Boundary conditions.
        Opm::ImplicitTransportDetails::NRReport  rpt;
        tsolver_.solve(grid_, tsrc_, dt, ctrl_, state, linsolver_, rpt);
        std::cout << rpt;

    }
}
