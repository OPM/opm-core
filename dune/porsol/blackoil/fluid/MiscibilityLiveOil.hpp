//===========================================================================
//                                                                           
// File: MiscibilityLiveOil.hpp                                               
//                                                                           
// Created: Wed Feb 10 09:08:09 2010                                         
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

#ifndef SINTEF_MISCIBILITYLIVEOIL_HEADER
#define SINTEF_MISCIBILITYLIVEOIL_HEADER

    /** Class for miscible live oil.
     *  Detailed description.
     */

#include "MiscibilityProps.hpp"
#include <dune/common/EclipseUnits.hpp>

namespace Opm
{
    class MiscibilityLiveOil : public MiscibilityProps
    {
    public:
	typedef std::vector<std::vector<std::vector<double> > > table_t;

	MiscibilityLiveOil(const table_t& pvto, const Dune::EclipseUnits& units);
	virtual ~MiscibilityLiveOil();

        virtual double getViscosity(int region, double press, const surfvol_t& surfvol) const;
	virtual double R   (int region, double press, const surfvol_t& surfvol) const;
	virtual double dRdp(int region, double press, const surfvol_t& surfvol) const;
        virtual double B   (int region, double press, const surfvol_t& surfvol) const;
	virtual double dBdp(int region, double press, const surfvol_t& surfvol) const;

    protected:
	// item:  1=B  2=mu;
	double miscible_oil(double press, const surfvol_t& surfvol, int item,
			    bool deriv = false) const;

	// PVT properties of live oil (with dissolved gas)
	std::vector<std::vector<double> > saturated_oil_table_;	
	std::vector<std::vector<std::vector<double> > > undersat_oil_tables_;
    };

}

#endif // SINTEF_MISCIBILITYLIVEOIL_HEADER

