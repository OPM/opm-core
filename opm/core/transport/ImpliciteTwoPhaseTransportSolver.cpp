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


#include "ImpliciteTwoPhaseTransportSolver.hpp"

ImpliciteTwoPhaseTransportSolver::ImpliciteTwoPhaseTransportSolver()
{
    ImpliciteTwoPhaseTransportSolver::ImpliciteTwoPhaseTransportSolver(TwoPhaseTransprotModel& model,
                const Opm::IncompPropertiesInterface& param)
        : tsolver_(model),
          grid_(grid),
          ctrl_(ctrl),
          props_(props),
          rock_comp_(rock_comp),
          wells_(wells),
    {
        tsrc_ = create_transport_source(2, 2);
        //src_(num_cells, 0.0);
    }

    void ImpliciteTwoPhaseTransportSolver::solve(const double* darcyflux,
                                       const double* porevolume,
                                       const double* source,
                                       const double dt,
                                       std::vector<double>& saturation)
    {
        if (rock_comp->isActive()) {
            computePorevolume(grid_->c_grid(), props->porosity(), *rock_comp, state.pressure(), porevol);
        }
        std::vector<double> src(num_cells, 0.0);
        //Opm::wellsToSrc(*wells->c_wells(), num_cells, src);
        Opm::computeTransportSource(*grid->c_grid(), src, state.faceflux(), 1.0,
                                            wells->c_wells(), well_state.perfRates(), src);
        double ssrc[]   = { 1.0, 0.0 };
        double ssink[]  = { 0.0, 1.0 };
        double zdummy[] = { 0.0, 0.0 };
        for (int cell = 0; cell < num_cells; ++cell) {
            clear_transport_source(tsrc);
            if (src[cell] > 0.0) {
                append_transport_source(cell, 2, 0, src[cell], ssrc, zdummy, tsrc);
            } else if (src[cell] < 0.0) {
                append_transport_source(cell, 2, 0, src[cell], ssink, zdummy, tsrc);
            }
        }
        // Boundary conditions.
        Opm::ImplicitTransportDetails::NRReport  rpt;
        tsolver.solve(grid_->c_grid(), tsrc, dt, ctrl_, state, linsolvet_, rpt);
        std::cout << rpt;

    }
}
