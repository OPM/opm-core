/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

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


#include "BlackoilPVT.hpp"
#include <dune/common/EclipseGridParser.hpp>
#include "MiscibilityDead.hpp"
#include "MiscibilityLiveOil.hpp"
#include "MiscibilityLiveGas.hpp"
#include "MiscibilityWater.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/linInt.hpp>

using namespace Dune;

namespace Opm
{


    void BlackoilPVT::init(const Dune::EclipseGridParser& parser)
    {
        typedef std::vector<std::vector<std::vector<double> > > table_t;
	region_number_ = 0;

	// Surface densities. Accounting for different orders in eclipse and our code.
	if (parser.hasField("DENSITY")) {
	    const int region_number = 0;
	    enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
	    const std::vector<double>& d = parser.getDENSITY().densities_[region_number];
            const double du = parser.units().density;
            using namespace Dune::unit;
	    densities_[Aqua] = convert::from(d[ECL_water], du);
	    densities_[Vapour] = convert::from(d[ECL_gas], du);
	    densities_[Liquid] = convert::from(d[ECL_oil], du);
	} else {
	    THROW("Input is missing DENSITY\n");
	}

        // Water PVT
        if (parser.hasField("PVTW")) {
            water_props_.reset(new MiscibilityWater(parser.getPVTW().pvtw_, parser.units()));
        } else {
            water_props_.reset(new MiscibilityWater(0.5*Dune::prefix::centi*Dune::unit::Poise)); // Eclipse 100 default 
        }

        // Oil PVT
        if (parser.hasField("PVDO")) {
            oil_props_.reset(new MiscibilityDead(parser.getPVDO().pvdo_, parser.units()));
        } else if (parser.hasField("PVTO")) {
            oil_props_.reset(new MiscibilityLiveOil(parser.getPVTO().pvto_, parser.units()));
        } else {
            THROW("Input is missing PVDO and PVTO\n");
        }

	// Gas PVT
        if (parser.hasField("PVDG")) {
            gas_props_.reset(new MiscibilityDead(parser.getPVDG().pvdg_, parser.units()));
        } else if (parser.hasField("PVTG")) {
            gas_props_.reset(new MiscibilityLiveGas(parser.getPVTG().pvtg_, parser.units()));
        } else {
            THROW("Input is missing PVDG and PVTG\n");
        }
    }

    BlackoilPVT::CompVec BlackoilPVT::surfaceDensities() const
    {
        return densities_;
    }

    double BlackoilPVT::getViscosity(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).getViscosity(region_number_, press, surfvol);
    }

    double BlackoilPVT::B(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).B(region_number_, press, surfvol);
    }

    double BlackoilPVT::dBdp(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).dBdp(region_number_, press, surfvol);
    }

    double BlackoilPVT::R(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).R(region_number_, press, surfvol);
    }

    double BlackoilPVT::dRdp(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).dRdp(region_number_, press, surfvol);
    }

    const MiscibilityProps& BlackoilPVT::propsForPhase(PhaseIndex phase) const
    {
        switch (phase) {
        case Aqua:
            return *water_props_;
        case Liquid:
            return *oil_props_;
        case Vapour:
            return *gas_props_;
        default:
            THROW("Unknown phase accessed: " << phase);
        }
    }

} // namespace Opm
