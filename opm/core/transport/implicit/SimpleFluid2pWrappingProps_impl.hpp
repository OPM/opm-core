/*===========================================================================
//
// File: SimpleFluid2pWrappingProps.cpp
//
// Author: hnil <hnil@sintef.no>
//
// Created: 15 Nov 2012
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
#ifndef SIMPLEFLUID2PWRAPPINGPROPS_IMPL_HPP
#define SIMPLEFLUID2PWRAPPINGPROPS_IMPL_HPP

#include <opm/core/transport/implicit/SimpleFluid2pWrappingProps.hpp>
#include <cassert>
#include <opm/core/utility/ErrorMacros.hpp>
namespace Opm{

    inline SimpleFluid2pWrappingProps::SimpleFluid2pWrappingProps(const Opm::IncompPropertiesInterface& props)
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

    inline double SimpleFluid2pWrappingProps::density(int phase) const
    {
        return props_.density()[phase];
    }

    template <class Sat,
              class Mob,
              class DMob>
    void SimpleFluid2pWrappingProps::mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const
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
    void SimpleFluid2pWrappingProps::pc(int c, const Sat& s, Pcap& pcap, DPcap& dpcap) const
    {
        double pcow[2];
        double dpcow[4];
        props_.capPress(1, &s[0], &c, pcow, dpcow);
        pcap = pcow[0];
        ASSERT(pcow[1] == 0.0);
        dpcap = dpcow[0];
        ASSERT(dpcow[1] == 0.0);
        ASSERT(dpcow[2] == 0.0);
        ASSERT(dpcow[3] == 0.0);
    }

    inline double SimpleFluid2pWrappingProps::s_min(int c) const
    {
        return smin_[2*c + 0];
    }

    inline double SimpleFluid2pWrappingProps::s_max(int c) const
    {
        return smax_[2*c + 0];
    }

}
#endif // SIMPLEFLUID2PWRAPPINGPROPS_IMPL_HPP
