//===========================================================================
//                                                                           
// File: MiscibilityLiveGas.hpp                                               
//                                                                           
// Created: Wed Feb 10 09:21:26 2010                                         
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

#ifndef SINTEF_MISCIBILITYLIVEGAS_HEADER
#define SINTEF_MISCIBILITYLIVEGAS_HEADER

    /** Class for miscible wet gas.
     *  Detailed description.
     */

#include "MiscibilityProps.hpp"

namespace Opm
{
    class MiscibilityLiveGas : public MiscibilityProps
    {
    public:
	typedef std::vector<std::vector<std::vector<double> > > table_t;

	MiscibilityLiveGas(const table_t& pvto);
	virtual ~MiscibilityLiveGas();

        virtual double getViscosity(int region, double press, const surfvol_t& surfvol) const;
	virtual double R(int region, double press, const surfvol_t& surfvol) const;
	virtual double dRdp(int region, double press, const surfvol_t& surfvol) const;
        virtual double B(int region, double press, const surfvol_t& surfvol) const;
	virtual double dBdp(int region, double press, const surfvol_t& surfvol) const;

    protected:
	// item:  1=B  2=mu;
	double miscible_gas(double press, const surfvol_t& surfvol, int item,
			    bool deriv = false) const;
	// PVT properties of wet gas (with vaporised oil)
	std::vector<std::vector<double> > saturated_gas_table_;	
	std::vector<std::vector<std::vector<double> > > undersat_gas_tables_;

    };

}

#endif // SINTEF_MISCIBILITYLIVEGAS_HEADER

