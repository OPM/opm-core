//===========================================================================
//                                                                           
// File: FluidMiscibilityThreePhase.cpp                                      
//                                                                           
// Created: Wed Feb 10 09:25:57 2010                                         
//                                                                           
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
//===========================================================================

#include "FluidMiscibilityThreePhase.hpp"
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


    void FluidMiscibilityThreePhase::init(const std::string& pvt_filename)
    {
        typedef std::vector<std::vector<std::vector<double> > > table_t;
	region_number_ = 0;
        Dune::EclipseGridParser eclipse_props(pvt_filename);

	// Surface densities. Accounting for different orders in eclipse and our code.
	if (eclipse_props.hasField("DENSITY")) {
	    const int region_number = 0;
	    enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
	    const std::vector<std::vector<double> > d_tmp = 
		eclipse_props.getDENSITY().densities_;
	    densities_[Aqua] = d_tmp[region_number][ECL_water];
	    densities_[Liquid] = d_tmp[region_number][ECL_oil];
	    densities_[Vapour] = d_tmp[region_number][ECL_gas];
	} else {
	    THROW("DENSITY not defined");
	}

        // Water PVT
        if (eclipse_props.hasField("PVTW")) {
            water_props_.reset(new MiscibilityWater(eclipse_props.getPVTW().pvtw_));
        } else {
            water_props_.reset(new MiscibilityWater(3e-4)); // Default is 0.3 cP.
        }

        // Oil PVT
        if (eclipse_props.hasField("PVDO")) {
            oil_props_.reset(new MiscibilityDead(eclipse_props.getPVDO().pvdo_));
        } else if (eclipse_props.hasField("PVTO")) {
            oil_props_.reset(new MiscibilityLiveOil(eclipse_props.getPVTO().pvto_));
        } else {
            THROW("File " << pvt_filename << " is missing PVDO and PVTO\n");
        }

	// Gas PVT
        if (eclipse_props.hasField("PVDG")) {
            gas_props_.reset(new MiscibilityDead(eclipse_props.getPVDG().pvdg_));
        } else if (eclipse_props.hasField("PVTG")) {
            gas_props_.reset(new MiscibilityLiveGas(eclipse_props.getPVTG().pvtg_));
        } else {
            THROW("File " << pvt_filename << " is missing PVDG and PVTG\n");
        }

	//SWOF
	if (eclipse_props.hasField("SWOF")) {
	    swof_ = eclipse_props.getSWOF().swof_;
	}

	//SGOF
	if (eclipse_props.hasField("SGOF")) {
	    sgof_ = eclipse_props.getSGOF().sgof_;
	}

    }
    FluidMiscibilityThreePhase::surfvol_t FluidMiscibilityThreePhase::surfaceDensities() const
    {
        return densities_;
    }

    double FluidMiscibilityThreePhase::getViscosity(double press, const surfvol_t& surfvol, PhaseNames phase) const
    {
        return propsForPhase(phase).getViscosity(region_number_, press, surfvol);
    }

    FluidMiscibilityThreePhase::surfvol_t 
    FluidMiscibilityThreePhase::getMobilities(double press, const surfvol_t& sat,
					      const surfvol_t& surfvol) const
    {
        if (swof_.empty() || sgof_.empty()) {
            THROW("The SWOF and SGOF keywords were not given, cannot compute mobilities. Try tracer_flow=true.");
        }
	surfvol_t mobilities;
	double sw = sat[Aqua];
	double sg = sat[Vapour];
	mobilities[Aqua]   = krw(sw)/getViscosity(press, surfvol, Aqua);
	mobilities[Liquid] = krow(sw)*krog(sg)/getViscosity(press, surfvol, Liquid);
	mobilities[Vapour] = krg(sg)/getViscosity(press, surfvol, Vapour);
	return mobilities;
    }

    double FluidMiscibilityThreePhase::B(double press, const surfvol_t& surfvol, PhaseNames phase) const
    {
        return propsForPhase(phase).B(region_number_, press, surfvol);
    }

    double FluidMiscibilityThreePhase::dBdp(double press, const surfvol_t& surfvol, PhaseNames phase) const
    {
        return propsForPhase(phase).dBdp(region_number_, press, surfvol);
    }

    double FluidMiscibilityThreePhase::R(double press, const surfvol_t& surfvol, PhaseNames phase) const
    {
        return propsForPhase(phase).R(region_number_, press, surfvol);
    }

    double FluidMiscibilityThreePhase::dRdp(double press, const surfvol_t& surfvol, PhaseNames phase) const
    {
        return propsForPhase(phase).dRdp(region_number_, press, surfvol);
    }

    const MiscibilityProps& FluidMiscibilityThreePhase::propsForPhase(PhaseNames phase) const
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

    // Stone's first model(Modified)
    double FluidMiscibilityThreePhase::kro(double sw, double sg) const
    {
	double so = 1.0-sw-sg;
	double fw = krow(sw);
	double fg = krog(sg);
	return so*fw*fg;
    }

    // kro partial derivative Sw
    double FluidMiscibilityThreePhase::dkro_dsw(double sw, double sg) const
    {
	double so = 1.0-sw-sg;
	double fw = krow(sw);
	double fg = krog(sg);
	double fwder = dkrow_dsw(sw);
	return fg*(-fw + so*fwder);
    }

    // kro partial derivative Sg
    double FluidMiscibilityThreePhase::dkro_dsg(double sw, double sg) const
    {
	double so = 1.0-sw-sg;
	double fw = krow(sw);
	double fg = krog(sg);
	double fgder = dkrog_dsg(sg);
	return fw*(-fg + so*fgder);
    }

    // kro partial derivatives
    void FluidMiscibilityThreePhase::dkro(double sw, double sg, double& dkro_dsw,
			     double& dkro_dsg) const
    {
	double so = 1.0-sw-sg;
	double fw = krow(sw);
	double fg = krog(sg);
	double fwder = dkrow_dsw(sw);
	double fgder = dkrog_dsg(sg);
	dkro_dsw = fg*(-fw + so*fwder);
	dkro_dsg = fw*(-fg + so*fgder);
    }

    // Water relative permeability
    double FluidMiscibilityThreePhase::krw  (double sw) const
    {
	return linearInterpolation(swof_[region_number_][0],
				   swof_[region_number_][1], sw);
    }

    // Oil relative permeability
    double FluidMiscibilityThreePhase::krow  (double sw) const
    {
	return linearInterpolation(swof_[region_number_][0],
				   swof_[region_number_][2], sw);
    }

    // krow derivative
    double FluidMiscibilityThreePhase::dkrow_dsw(double sw) const
    {
	return linearInterpolDerivative(swof_[region_number_][0],
					swof_[region_number_][2], sw);
    }

    // Water-oil capillary pressure
    double  FluidMiscibilityThreePhase::Pcow  (double sw) const
    {
	return linearInterpolation(swof_[region_number_][0],
				   swof_[region_number_][3], sw);
    }

    // Gas relative permeability
    double FluidMiscibilityThreePhase::krg  (double sg) const
    {
	return linearInterpolation(sgof_[region_number_][0],
				   sgof_[region_number_][1], sg);
    }

    // Oil relative permeability
    double FluidMiscibilityThreePhase::krog  (double sg) const
    {
	return linearInterpolation(sgof_[region_number_][0],
				   sgof_[region_number_][2], sg);
    }

    // krog derivative
    double FluidMiscibilityThreePhase::dkrog_dsg(double sg) const
    {
	return linearInterpolDerivative(sgof_[region_number_][0],
					sgof_[region_number_][2], sg);
    }

    // Oil-gas capillary pressure
    double FluidMiscibilityThreePhase::Pcog(double sg) const
    {
	return linearInterpolation(sgof_[region_number_][0],
				   sgof_[region_number_][3], sg);
    }


} // namespace Opm
