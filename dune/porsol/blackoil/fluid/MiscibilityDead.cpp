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

#include <algorithm>
#include "MiscibilityDead.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/linInt.hpp>
#include <dune/common/Units.hpp>

using namespace std;
using namespace Dune;

namespace Opm
{

    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    MiscibilityDead::MiscibilityDead(const table_t& pvd_table, const Dune::EclipseUnits& units)
	: pvdx_(pvd_table)
    {
	const int region_number = 0;
	if (pvd_table.size() != 1) {
	    THROW("More than one PVT-region");
	}
	// Convert units
	const int sz =  pvdx_[region_number][0].size();
        using namespace Dune::unit;
	for (int i=0; i<sz; ++i) {
	    pvdx_[region_number][0][i] = convert::from(pvdx_[region_number][0][i], units.pressure);
	    pvdx_[region_number][2][i] = convert::from(pvdx_[region_number][2][i], units.viscosity);
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
