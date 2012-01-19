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


#include <opm/core/fluid/blackoil/BlackoilPVT.hpp>
#include <opm/core/fluid/blackoil/MiscibilityDead.hpp>
#include <opm/core/fluid/blackoil/MiscibilityLiveOil.hpp>
#include <opm/core/fluid/blackoil/MiscibilityLiveGas.hpp>
#include <opm/core/fluid/blackoil/MiscibilityWater.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linInt.hpp>

using namespace Opm;

namespace Opm
{


    void BlackoilPVT::init(const EclipseGridParser& parser)
    {
        typedef std::vector<std::vector<std::vector<double> > > table_t;
	region_number_ = 0;

	// Surface densities. Accounting for different orders in eclipse and our code.
	if (parser.hasField("DENSITY")) {
	    const int region_number = 0;
	    enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
	    const std::vector<double>& d = parser.getDENSITY().densities_[region_number];
	    densities_[Aqua]   = d[ECL_water];
	    densities_[Vapour] = d[ECL_gas];
	    densities_[Liquid] = d[ECL_oil];
	} else {
	    THROW("Input is missing DENSITY\n");
	}

        // Water PVT
        if (parser.hasField("PVTW")) {
            water_props_.reset(new MiscibilityWater(parser.getPVTW().pvtw_));
        } else {
            water_props_.reset(new MiscibilityWater(0.5*Opm::prefix::centi*Opm::unit::Poise)); // Eclipse 100 default 
        }

        // Oil PVT
        if (parser.hasField("PVDO")) {
            oil_props_.reset(new MiscibilityDead(parser.getPVDO().pvdo_));
        } else if (parser.hasField("PVTO")) {
            oil_props_.reset(new MiscibilityLiveOil(parser.getPVTO().pvto_));
        } else if (parser.hasField("PVCDO")) {
            oil_props_.reset(new MiscibilityWater(parser.getPVCDO().pvcdo_));
        } else {
            THROW("Input is missing PVDO and PVTO\n");
        }

	// Gas PVT
        if (parser.hasField("PVDG")) {
            gas_props_.reset(new MiscibilityDead(parser.getPVDG().pvdg_));
        } else if (parser.hasField("PVTG")) {
            gas_props_.reset(new MiscibilityLiveGas(parser.getPVTG().pvtg_));
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

    void BlackoilPVT::getViscosity(const std::vector<PhaseVec>& pressures,
                                   const std::vector<CompVec>& surfvol,
                                   std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).getViscosity(pressures, surfvol, phase, data1_);
            for (int i = 0; i < num; ++i) {
                output[i][phase] = data1_[i];
            }
        }
    }

    void BlackoilPVT::B(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).B(pressures, surfvol, phase, data1_);
            for (int i = 0; i < num; ++i) {
                output[i][phase] = data1_[i];
            }
        }
    }

    void BlackoilPVT::dBdp(const std::vector<PhaseVec>& pressures,
                           const std::vector<CompVec>& surfvol,
                           std::vector<PhaseVec>& output_B,
                           std::vector<PhaseVec>& output_dBdp) const
    {
        int num = pressures.size();
        output_B.resize(num);
        output_dBdp.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).dBdp(pressures, surfvol, phase, data1_, data2_);
            for (int i = 0; i < num; ++i) {
                output_B[i][phase] = data1_[i];
                output_dBdp[i][phase] = data2_[i];
            }
        }
    }

    void BlackoilPVT::R(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).R(pressures, surfvol, phase, data1_);
            for (int i = 0; i < num; ++i) {
                output[i][phase] = data1_[i];
            }
        }
    }

    void BlackoilPVT::dRdp(const std::vector<PhaseVec>& pressures,
                           const std::vector<CompVec>& surfvol,
                           std::vector<PhaseVec>& output_R,
                           std::vector<PhaseVec>& output_dRdp) const
    {
        int num = pressures.size();
        output_R.resize(num);
        output_dRdp.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).dRdp(pressures, surfvol, phase, data1_, data2_);
            for (int i = 0; i < num; ++i) {
                output_R[i][phase] = data1_[i];
                output_dRdp[i][phase] = data2_[i];
            }
        }
    }

} // namespace Opm
