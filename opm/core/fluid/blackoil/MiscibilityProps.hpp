//===========================================================================
//                                                                           
// File: MiscibilityProps.hpp                                                 
//                                                                           
// Created: Wed Feb 10 09:04:35 2010                                         
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

#ifndef SINTEF_MISCIBILITYPROPS_HEADER
#define SINTEF_MISCIBILITYPROPS_HEADER

    /** Base class for properties of fluids and rocks.
     *  Detailed description.
     */

#include "BlackoilDefs.hpp"
#include <vector>
#include <tr1/array>

namespace Opm
{


    class MiscibilityProps : public BlackoilDefs
    {
    public:
	typedef CompVec surfvol_t;

	MiscibilityProps();
	virtual ~MiscibilityProps();

        virtual double getViscosity(int region, double press, const surfvol_t& surfvol) const = 0;
        virtual double B   (int region, double press, const surfvol_t& surfvol) const = 0;
	virtual double dBdp(int region, double press, const surfvol_t& surfvol) const = 0;
	virtual double R   (int region, double press, const surfvol_t& surfvol) const = 0;
	virtual double dRdp(int region, double press, const surfvol_t& surfvol) const = 0;

        virtual void getViscosity(const std::vector<PhaseVec>& pressures,
                                  const std::vector<CompVec>& surfvol,
                                  int phase,
                                  std::vector<double>& output) const = 0;
        virtual void B(const std::vector<PhaseVec>& pressures,
                       const std::vector<CompVec>& surfvol,
                       int phase,
                       std::vector<double>& output) const = 0;
        virtual void dBdp(const std::vector<PhaseVec>& pressures,
                          const std::vector<CompVec>& surfvol,
                          int phase,
                          std::vector<double>& output_B,
                          std::vector<double>& output_dBdp) const = 0;
        virtual void R(const std::vector<PhaseVec>& pressures,
                       const std::vector<CompVec>& surfvol,
                       int phase,
                       std::vector<double>& output) const = 0;
        virtual void dRdp(const std::vector<PhaseVec>& pressures,
                          const std::vector<CompVec>& surfvol,
                          int phase,
                          std::vector<double>& output_R,
                          std::vector<double>& output_dRdp) const = 0;
    };

} // namespace Opm

#endif // SINTEF_MISCIBILITYPROPS_HEADER

