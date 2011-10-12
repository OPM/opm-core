/*===========================================================================
//
// File: SimpleFluid2pWrapper.hpp
//
// Created: 2011-09-30 11:38:28+0200
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

#ifndef OPM_SimpleFluid2pWrapper_HPP_HEADER
#define OPM_SimpleFluid2pWrapper_HPP_HEADER

#include <array>
#include <cmath>

namespace Opm {
    template <class ReservoirProperties>
    class SimpleFluid2pWrapper {
    public:
        SimpleFluid2pWrapper	(const ReservoirProperties& resprop)
        {
           resprop_ = resprop;
        }

        double density(int p) const { return 0; }

        template <class Sat ,
                  class Mob ,
                  class DMob>
        void
        mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const {
            const double s1 = s[0];
            const double s2 = 1 - s1;
            mob[0]  = resprop_.mobilityFirstPhase(c, s[0]) ;
            mob[1]  = resprop_.mobilitySecondPhase(c, s[0]) ;

            dmob[0] =  resprop_.dmobilityFirstPhase(c, s[0]) ;
            dmob[3] = -resprop_.dmobilitySecondPhase(c, s[0]) ;
            dmob[1] = dmob[2] = 0.0;
        }

    private:
        ReservoirProperties& resprop_;
        std::array<double, 2> mu_ ;
        std::array<double, 2> rho_;
    };
}

#endif  /* OPM_SimpleFluid2pWrapper_HPP_HEADER */
