//===========================================================================
//
// File: MiscibilityWater.hpp
//
// Created: Tue May 18 10:26:13 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bjørn Spjelkavik    <bsp@sintef.no>
//
// $Date$
//
// $Revision$
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

#ifndef OPENRS_MISCIBILITYWATER_HEADER
#define OPENRS_MISCIBILITYWATER_HEADER

#include "MiscibilityProps.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/Units.hpp>

// Forward declaration.
class PVTW;

namespace Opm
{
    class MiscibilityWater : public MiscibilityProps
    {
    public:
        typedef std::vector<std::vector<double> > table_t;
	MiscibilityWater(const table_t& pvtw, const Dune::EclipseUnits& units)
        {
	    const int region_number = 0;
	    if (pvtw.size() != 1) {
		THROW("More than one PVD-region");
	    }
            using namespace Dune::unit;
            ref_press_ = convert::from(pvtw[region_number][0], units.pressure);
            ref_B_ = convert::from(pvtw[region_number][1], units.liqvol_r/units.liqvol_s);
            comp_ = convert::from(pvtw[region_number][2], units.compressibility);
            viscosity_ = convert::from(pvtw[region_number][3], units.viscosity);
            if (pvtw[region_number].size() > 4 && pvtw[region_number][4] != 0.0) {
                THROW("MiscibilityWater does not support 'viscosibility'.");
            }
        }
	MiscibilityWater(double visc)
            : ref_press_(0.0),
              ref_B_(1.0),
              comp_(0.0),
              viscosity_(visc)
        {
        }
	virtual ~MiscibilityWater()
        {
        }

        virtual double getViscosity(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return viscosity_;
        }
        virtual double B(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 1.0;
        }
	virtual double dBdp(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }
	virtual double R(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }
	virtual double dRdp(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }
    private:
        double ref_press_;
        double ref_B_;
        double comp_;
        double viscosity_;
    };

}

#endif // OPENRS_MISCIBILITYWATER_HEADER
