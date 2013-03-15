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


#include <opm/core/transport/TransportSolverTwophaseImplicit.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/utility/miscUtilities.hpp>

namespace Opm
{

    TransportSolverTwophaseImplicit::TransportSolverTwophaseImplicit(
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
    }

    TransportSolverTwophaseImplicit::~TransportSolverTwophaseImplicit()
    {
        destroy_transport_source(tsrc_);
    }

    /// Solve for saturation at next timestep.
    /// \param[in]      porevolume   Array of pore volumes.
    /// \param[in]      source       Transport source term. For interpretation see Opm::computeTransportSource().
    /// \param[in]      dt           Time step.
    /// \param[in, out] state        Reservoir state. Calling solve() will read state.faceflux() and
    ///                              read and write state.saturation().
    void TransportSolverTwophaseImplicit::solve(const double* porevolume,
                                                const double* source,
                                                const double dt,
                                                TwophaseState& state)
    {
        std::vector<double> porevol;
        if (rock_comp_.isActive()) {
            computePorevolume(grid_, props_.porosity(), rock_comp_, state.pressure(), porevol);
        }
        double ssrc[] = { 1.0, 0.0 };
        double dummy[] = { 0.0, 0.0 };
        clear_transport_source(tsrc_);
        for (int cell = 0; cell < grid_.number_of_cells; ++cell) {
            if (source[cell] > 0.0) {
                append_transport_source(cell, 2, 0, source[cell], ssrc, dummy, tsrc_);
            } else if (source[cell] < 0.0) {
                append_transport_source(cell, 2, 0, source[cell], dummy, dummy, tsrc_);
            }
        }
        Opm::ImplicitTransportDetails::NRReport  rpt;
        tsolver_.solve(grid_, tsrc_, dt, ctrl_, state, linsolver_, rpt);
        std::cout << rpt;
    }

} // namespace Opm
