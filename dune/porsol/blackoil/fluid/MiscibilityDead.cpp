//===========================================================================
//                                                                           
// File: MiscibilityDead.cpp                                                  
//                                                                           
// Created: Wed Feb 10 09:06:04 2010                                         
//                                                                           
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
//===========================================================================

#include <algorithm>
#include "MiscibilityDead.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/linInt.hpp>

using namespace std;
using namespace Dune;

namespace samcode
{

    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    MiscibilityDead::MiscibilityDead(const table_t& pvd_table)
	: pvdx_(pvd_table)
    {
	const int region_number = 0;
	if (pvd_table.size() != 1) {
	    THROW("More than one PVT-region");
	}
	// Convert units
	const double bar = 1e5;
	const double VISCOSITY_UNIT = 1e-3;
	const int sz =  pvdx_[region_number][0].size();
	for (int i=0; i<sz; ++i) {
	    pvdx_[region_number][0][i] *= bar;  // Pressure
	    pvdx_[region_number][2][i] *= VISCOSITY_UNIT;
	}

	// Interpolate 1/B 
	for (int i=0; i<sz; ++i) {
	    pvdx_[region_number][1][i] = 1.0/pvdx_[region_number][1][i];
	}
    }

    // Destructor
    MiscibilityDead::~MiscibilityDead()
    {
    }

    double MiscibilityDead::getViscosity(int region, double press, const surfvol_t& /*surfvol*/) const
    {
	return linearInterpolationExtrap(pvdx_[region][0],
					 pvdx_[region][2], press);
    }


    double MiscibilityDead::B(int region, double press, const surfvol_t& /*surfvol*/) const
    {
	// Interpolate 1/B 
	return 1.0/linearInterpolationExtrap(pvdx_[region][0],
					     pvdx_[region][1], press);

	//return linearInterpolationExtrap(pvdx_[region][0],
	//pvdx_[region_number_][1], press);
    }

    double MiscibilityDead::dBdp(int region, double press, const surfvol_t& /*surfvol*/) const
    {
	// Interpolate 1/B
	surfvol_t dummy_surfvol;
	double Bg = B(region, press, dummy_surfvol);
	return -Bg*Bg*
	    linearInterpolDerivative(pvdx_[region][0],
				     pvdx_[region][1], press);

	//return linearInterpolDerivative(pvdx_[region][0],
	//			pvdx_[region][1], press);
    }

    double MiscibilityDead::R(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
    {
        return 0.0;
    }

    double MiscibilityDead::dRdp(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
    {
        return 0.0;
    }

}
