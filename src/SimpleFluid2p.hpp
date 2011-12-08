/*===========================================================================
//
// File: SimpleFluid2p.hpp
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

#ifndef OPM_SIMPLEFLUID2P_HPP_HEADER
#define OPM_SIMPLEFLUID2P_HPP_HEADER

#include <array>
#include <cmath>

namespace Opm {
    template <int n = 2>
    class SimpleFluid2p {
    public:
        SimpleFluid2p(const std::array<double, 2>& mu,
                      const std::array<double, 2>& rho)
            : mu_(mu), rho_(rho)
        {
        }

        double density(int p) const { return rho_[p]; }

        template <class Sat ,
                  class Mob ,
                  class DMob>
        void
        mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const {
            (void) c;           // Unused

            const double s1 = s[0];
            const double s2 = 1 - s1;

            if (n == 1) {
                mob [      0] = s1;   mob [      1] = s2;

                dmob[0*2 + 0] = 1 ;   dmob[1*2 + 1] = 1 ;
            } else if (n == 2) {
                mob [      0] = s1 * s1;   mob [      1] = s2 * s2;

                dmob[0*2 + 0] = 2 * s1 ;   dmob[1*2 + 1] = 2 * s2 ;
            } else {
                mob [      0] = pow(s1, double(n));
                mob [      1] = pow(s2, double(n));

                dmob[0*2 + 0] = double(n) * pow(s1, double(n) - 1);
                dmob[1*2 + 1] = double(n) * pow(s2, double(n) - 1);
            }

            mob[0] /= mu_[0];  dmob[0*2 + 0] /= mu_[0];
            mob[1] /= mu_[1];  dmob[1*2 + 1] /= mu_[1];

            dmob[0*2 + 1] = dmob[1*2 + 0] = 0;
        }

        template <class Sat  ,
                  class Pcap ,
                  class DPcap>
        void
        pc(int c, const Sat& s, Pcap& pcap, DPcap& dpcap) const {
            (void) c;           // Unused
            (void) s;           // Unused

            pcap  = 0.0;
            dpcap = 0.0;
        }

        double s_min(int c) const { (void) c; return 0.0; }
        double s_max(int c) const { (void) c; return 1.0; }

    private:
        std::array<double, 2> mu_ ;
        std::array<double, 2> rho_;
    };
}

#endif  /* OPM_SIMPLEFLUID2P_HPP_HEADER */
