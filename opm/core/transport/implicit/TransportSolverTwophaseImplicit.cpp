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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/transport/implicit/TransportSolverTwophaseImplicit.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/utility/miscUtilities.hpp>

#include <iostream>

namespace Opm
{

    TransportSolverTwophaseImplicit::TransportSolverTwophaseImplicit(
            const UnstructuredGrid& grid,
            const Opm::IncompPropertiesInterface& props,
            const std::vector<double>& porevol,
            const double* gravity,
            const std::vector<double>& half_trans,
            const parameter::ParameterGroup& param)
        : fluid_(props),
          model_(fluid_, grid, porevol, gravity, param.getDefault("guess_old_solution", false)),
          tsolver_(model_),
          grid_(grid),
          props_(props)
    {
        ctrl_.max_it = param.getDefault("max_it", 20);
        ctrl_.verbosity = param.getDefault("verbosity", 0);
        ctrl_.max_it_ls = param.getDefault("max_it_ls", 5);
        model_.initGravityTrans(grid_, half_trans);
        tsrc_ = create_transport_source(2, 2);
        initial_porevolume_cell0_ = porevol[0];
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
        // A very crude check for constant porosity (i.e. no rock-compressibility).
        if (porevolume[0] != initial_porevolume_cell0_) {
            THROW("Detected changed pore volumes, but solver cannot handle rock compressibility.");
        }
        double ssrc[] = { 1.0, 0.0 };
        double dummy[] = { 0.0, 0.0 };
        clear_transport_source(tsrc_);
        const int num_phases = 2;
        for (int cell = 0; cell < grid_.number_of_cells; ++cell) {
            int success = 1;
            if (source[cell] > 0.0) {
                success = append_transport_source(cell, num_phases, state.pressure()[cell], source[cell], ssrc, dummy, tsrc_);
            } else if (source[cell] < 0.0) {
                success = append_transport_source(cell, num_phases, state.pressure()[cell], source[cell], dummy, dummy, tsrc_);
            }
            if (!success) {
                THROW("Failed building TransportSource struct.");
            }
        }
        Opm::ImplicitTransportDetails::NRReport  rpt;
        tsolver_.solve(grid_, tsrc_, dt, ctrl_, state, linsolver_, rpt);
        std::cout << rpt;
    }

} // namespace Opm
