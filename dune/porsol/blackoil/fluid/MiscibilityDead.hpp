//===========================================================================
//                                                                           
// File: MiscibilityDead.hpp                                                  
//                                                                           
// Created: Wed Feb 10 09:05:47 2010                                         
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

#ifndef SINTEF_MISCIBILITYDEAD_HEADER
#define SINTEF_MISCIBILITYDEAD_HEADER

/** Class for immiscible dead oil and dry gas.
 *  Detailed description.
 */

#include "MiscibilityProps.hpp"
#include <dune/porsol/common/UniformTableLinear.hpp>

namespace Opm
{
    class MiscibilityDead : public MiscibilityProps
    {
    public:
	typedef std::vector<std::vector<std::vector<double> > > table_t;

	MiscibilityDead(const table_t& pvd_table);
	virtual ~MiscibilityDead();

        virtual double getViscosity(int region, double press, const surfvol_t& surfvol) const;
        virtual double B(int region, double press, const surfvol_t& surfvol) const;
	virtual double dBdp(int region, double press, const surfvol_t& surfvol) const;
	virtual double R(int region, double press, const surfvol_t& surfvol) const;
	virtual double dRdp(int region, double press, const surfvol_t& surfvol) const;
        

    private:
	// PVT properties of dry gas or dead oil
        Dune::utils::UniformTableLinear<double> one_over_B_;
        Dune::utils::UniformTableLinear<double> viscosity_;
    };

}

#endif // SINTEF_MISCIBILITYDEAD_HEADER

