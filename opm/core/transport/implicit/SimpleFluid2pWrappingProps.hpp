/*===========================================================================
//
// File: SimpleFluid2pWrappingProps.hpp
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


#ifndef SIMPLEFLUID2PWRAPPINGPROPS_HPP
#define SIMPLEFLUID2PWRAPPINGPROPS_HPP

#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <vector>

namespace Opm{
    class SimpleFluid2pWrappingProps
    {
    public:
        SimpleFluid2pWrappingProps(const Opm::IncompPropertiesInterface& props);

        double density(int phase) const;


        template <class Sat,
                  class Mob,
                  class DMob>
        void mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const;


        template <class Sat,
                  class Pcap,
                  class DPcap>
        void pc(int c, const Sat& s, Pcap& pcap, DPcap& dpcap) const;

        double s_min(int c) const;


        double s_max(int c) const;


    private:
        const Opm::IncompPropertiesInterface& props_;
        std::vector<double> smin_;
        std::vector<double> smax_;
    };
}

#include <opm/core/transport/implicit/SimpleFluid2pWrappingProps_impl.hpp>
#endif // SIMPLEFLUID2PWRAPPINGPROPS_HPP
