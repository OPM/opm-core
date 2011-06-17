//===========================================================================
//                                                                           
// File: MiscibiltyLiveOil.cpp                                               
//                                                                           
// Created: Wed Feb 10 09:08:25 2010                                         
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
#include "MiscibilityLiveOil.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/linInt.hpp>

using namespace std;
using namespace Dune;

namespace Opm
{


    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    MiscibilityLiveOil::MiscibilityLiveOil(const table_t& pvto)
    {
	// OIL, PVTO
	const int region_number = 0;
	if (pvto.size() != 1) {
	    THROW("More than one PVD-region");
	}
	saturated_oil_table_.resize(4);
	const int sz =  pvto[region_number].size();
	for (int k=0; k<4; ++k) {
	    saturated_oil_table_[k].resize(sz);
	}
	for (int i=0; i<sz; ++i) {
	    saturated_oil_table_[0][i] = pvto[region_number][i][1]; // p
	    saturated_oil_table_[1][i] = 1.0/pvto[region_number][i][2]; // 1/Bo
	    saturated_oil_table_[2][i] = pvto[region_number][i][3];   // mu_o
	    saturated_oil_table_[3][i] = pvto[region_number][i][0];     // Rs
	}
	
	undersat_oil_tables_.resize(sz);
	for (int i=0; i<sz; ++i) {
	    undersat_oil_tables_[i].resize(3);
	    int tsize = (pvto[region_number][i].size() - 1)/3;
	    undersat_oil_tables_[i][0].resize(tsize);
	    undersat_oil_tables_[i][1].resize(tsize);
	    undersat_oil_tables_[i][2].resize(tsize);
	    for (int j=0, k=0; j<tsize; ++j) {
		undersat_oil_tables_[i][0][j] = pvto[region_number][i][++k];  // p
		undersat_oil_tables_[i][1][j] = 1.0/pvto[region_number][i][++k];  // 1/Bo
		undersat_oil_tables_[i][2][j] = pvto[region_number][i][++k];  // mu_o
	    }
	}
	
	
	// Fill in additional entries in undersaturated tables by interpolating/extrapolating 1/Bo and mu_o ...
	int iPrev = -1;
	int iNext = 1;
	while (undersat_oil_tables_[iNext][0].size() < 2) {
		++iNext;
	}
	assert(iNext < sz); 
	for (int i=0; i<sz; ++i) {
		if (undersat_oil_tables_[i][0].size() > 1) {
			iPrev = i;
			continue;
		}
		
		bool flagPrev = (iPrev >= 0); 
		bool flagNext = true;			
		if (iNext < i) {
			iPrev = iNext;
			flagPrev = true;
			iNext = i+1;
			while (undersat_oil_tables_[iNext][0].size() < 2) {
				++iNext;
			}
		}				
		double slopePrevBinv = 0.0;
		double slopePrevVisc = 0.0;
		double slopeNextBinv = 0.0;
		double slopeNextVisc = 0.0;
		while (flagPrev || flagNext) {
			double pressure0 = undersat_oil_tables_[i][0].back();
			double pressure = 1.0e47;
			if (flagPrev) {
				std::vector<double>::iterator itPrev = upper_bound(undersat_oil_tables_[iPrev][0].begin(),
				                                                   undersat_oil_tables_[iPrev][0].end(),pressure0+1.);
				if (itPrev == undersat_oil_tables_[iPrev][0].end()) {
					--itPrev; // Extrapolation ...
				} else if (itPrev == undersat_oil_tables_[iPrev][0].begin()) {
					++itPrev;
				}
				if (itPrev == undersat_oil_tables_[iPrev][0].end()-1) {
					flagPrev = false; // Last data set for "prev" ...
				}
				double dPPrev = *itPrev - *(itPrev-1);
				pressure = *itPrev;
				int index = int(itPrev - undersat_oil_tables_[iPrev][0].begin());
				slopePrevBinv = (undersat_oil_tables_[iPrev][1][index] - undersat_oil_tables_[iPrev][1][index-1])/dPPrev;
				slopePrevVisc = (undersat_oil_tables_[iPrev][2][index] - undersat_oil_tables_[iPrev][2][index-1])/dPPrev;
			}
			if (flagNext) {
				std::vector<double>::iterator itNext = upper_bound(undersat_oil_tables_[iNext][0].begin(),
				                                                   undersat_oil_tables_[iNext][0].end(),pressure0+1.);
				if (itNext == undersat_oil_tables_[iNext][0].end()) {
					--itNext; // Extrapolation ...
				} else if (itNext == undersat_oil_tables_[iNext][0].begin()) {
					++itNext;
				}
				if (itNext == undersat_oil_tables_[iNext][0].end()-1) {
					flagNext = false; // Last data set for "next" ...
				}
				double dPNext = *itNext - *(itNext-1);
				if (flagPrev) {
					pressure = std::min(pressure,*itNext);
				} else {
					pressure = *itNext;
				}			
				int index = int(itNext - undersat_oil_tables_[iNext][0].begin());
				slopeNextBinv = (undersat_oil_tables_[iNext][1][index] - undersat_oil_tables_[iNext][1][index-1])/dPNext;
				slopeNextVisc = (undersat_oil_tables_[iNext][2][index] - undersat_oil_tables_[iNext][2][index-1])/dPNext;
			}	
			double dP = pressure - pressure0;
			if (iPrev >= 0) {
				double w = (saturated_oil_table_[3][i] - saturated_oil_table_[3][iPrev]) /
				           (saturated_oil_table_[3][iNext] - saturated_oil_table_[3][iPrev]);
				undersat_oil_tables_[i][0].push_back(pressure0+dP);
				undersat_oil_tables_[i][1].push_back(undersat_oil_tables_[i][1].back() +
				                                     dP*(slopePrevBinv+w*(slopeNextBinv-slopePrevBinv)));
				undersat_oil_tables_[i][2].push_back(undersat_oil_tables_[i][2].back() +
				                                     dP*(slopePrevVisc+w*(slopeNextVisc-slopePrevVisc)));
			} else {
				undersat_oil_tables_[i][0].push_back(pressure0+dP);
				undersat_oil_tables_[i][1].push_back(undersat_oil_tables_[i][1].back()+dP*slopeNextBinv);
				undersat_oil_tables_[i][2].push_back(undersat_oil_tables_[i][2].back()+dP*slopeNextVisc);
			} 
		}
	}
    }
    // Destructor
     MiscibilityLiveOil::~MiscibilityLiveOil()
    {
    }

    double MiscibilityLiveOil::getViscosity(int /*region*/, double press, const surfvol_t& surfvol) const
    {
	return miscible_oil(press, surfvol, 2, false);
    }

    // Dissolved gas-oil ratio   
    double MiscibilityLiveOil::R(int /*region*/, double press, const surfvol_t& surfvol) const
    {
        if (surfvol[Vapour] == 0.0) {
            return 0.0;
        }	
	double R = linearInterpolationExtrap(saturated_oil_table_[0],
					     saturated_oil_table_[3], press);
	double maxR = surfvol[Vapour]/surfvol[Liquid];
	if (R < maxR ) {  // Saturated case
	    return R;
	} else {
	    return maxR;  // Undersaturated case
	}
    }

    //  Dissolved gas-oil ratio derivative
    double MiscibilityLiveOil::dRdp(int /*region*/, double press, const surfvol_t& surfvol) const
    {
	double R = linearInterpolationExtrap(saturated_oil_table_[0],
					     saturated_oil_table_[3], press);
	double maxR = surfvol[Vapour]/surfvol[Liquid];
	if (R < maxR ) {  // Saturated case
	    return linearInterpolDerivative(saturated_oil_table_[0],
					    saturated_oil_table_[3],
					    press);
	} else {
	    return 0.0;  // Undersaturated case
	}	
    }

    double MiscibilityLiveOil::B(int /*region*/, double press, const surfvol_t& surfvol) const
    {
        // if (surfvol[Liquid] == 0.0) return 1.0; // To handle no-oil case.
	return 1.0/miscible_oil(press, surfvol, 1, false);
    }

    double MiscibilityLiveOil::dBdp(int region, double press, const surfvol_t& surfvol) const
    {	
        // if (surfvol[Liquid] == 0.0) return 0.0; // To handle no-oil case.
	double Bo = B(region, press, surfvol);
	return -Bo*Bo*miscible_oil(press, surfvol, 1, true);
    }

    double MiscibilityLiveOil::miscible_oil(double press, const surfvol_t& surfvol,
					    int item, bool deriv) const
    {
	int section;
	double R = linearInterpolationExtrap(saturated_oil_table_[0],
					     saturated_oil_table_[3],
					     press, section);
	double maxR = (surfvol[Liquid] == 0.0) ? 0.0 : surfvol[Vapour]/surfvol[Liquid];
	if (deriv) {
	    if (R < maxR ) {  // Saturated case
		return linearInterpolDerivative(saturated_oil_table_[0],
						saturated_oil_table_[item],
						press);
	    } else {  // Undersaturated case
		int is = tableIndex(saturated_oil_table_[3], maxR);
		double w = (maxR - saturated_oil_table_[3][is]) /
		    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                ASSERT(undersat_oil_tables_[is][0].size() >= 2);
                ASSERT(undersat_oil_tables_[is+1][0].size() >= 2);
		double val1 =
		    linearInterpolDerivative(undersat_oil_tables_[is][0],
					     undersat_oil_tables_[is][item],
					     press);
		double val2 = 
		    linearInterpolDerivative(undersat_oil_tables_[is+1][0],
					     undersat_oil_tables_[is+1][item],
					     press);
		double val = val1 + w*(val2 - val1);
		return val;
	    }
	} else {
	    if (R < maxR ) {  // Saturated case
		return linearInterpolationExtrap(saturated_oil_table_[0],
						 saturated_oil_table_[item],
						 press);
	    } else {  // Undersaturated case
		// Interpolate between table sections
                int is = tableIndex(saturated_oil_table_[3], maxR);
		double w = (maxR - saturated_oil_table_[3][is]) /
		    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                ASSERT(undersat_oil_tables_[is][0].size() >= 2);
                ASSERT(undersat_oil_tables_[is+1][0].size() >= 2);
		double val1 =
		    linearInterpolationExtrap(undersat_oil_tables_[is][0],
					      undersat_oil_tables_[is][item],
					      press);
		double val2 = 
		    linearInterpolationExtrap(undersat_oil_tables_[is+1][0],
					      undersat_oil_tables_[is+1][item],
					      press);
		double val = val1 + w*(val2 - val1);
		return val;
	    }
	}
    }

} // namespace Opm
