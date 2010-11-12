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
	const double bar = 1e5;
	const double VISCOSITY_UNIT = 1e-3;
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
	    saturated_oil_table_[0][i] = pvto[region_number][i][1]*bar; // p
	    saturated_oil_table_[1][i] = 1.0/pvto[region_number][i][2]; // 1/Bo
	    saturated_oil_table_[2][i] = pvto[region_number][i][3]*VISCOSITY_UNIT;   // mu_o
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
		undersat_oil_tables_[i][0][j] = pvto[region_number][i][++k]*bar;  // p
		undersat_oil_tables_[i][1][j] = 1.0/pvto[region_number][i][++k];  // 1/Bo
		undersat_oil_tables_[i][2][j] = pvto[region_number][i][++k]*VISCOSITY_UNIT;  // mu_o
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
        if (surfvol[Liquid] == 0.0) return 1.0; // To handle no-oil case.
	return 1.0/miscible_oil(press, surfvol, 1, false);
    }

    double MiscibilityLiveOil::dBdp(int region, double press, const surfvol_t& surfvol) const
    {	
        if (surfvol[Liquid] == 0.0) return 0.0; // To handle no-oil case.
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
	double maxR = surfvol[Vapour]/surfvol[Liquid];
	if (deriv) {
	    if (R < maxR ) {  // Saturated case
		return linearInterpolDerivative(saturated_oil_table_[0],
						saturated_oil_table_[item],
						press);
	    } else {  // Undersaturated case
		int is = tableIndex(saturated_oil_table_[3], maxR);
		if (undersat_oil_tables_[is][0].size() < 2) {
		    double val = (saturated_oil_table_[item][is+1]
				  - saturated_oil_table_[item][is]) /
			(saturated_oil_table_[0][is+1] -
			 saturated_oil_table_[0][is]);

		    return val;
		}
		double w = (maxR - saturated_oil_table_[3][is]) /
		    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);

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
		int is = tableIndex(saturated_oil_table_[3], maxR);

		// Extrapolate from first table section
		if (is == 0 && press < saturated_oil_table_[0][0]) {
		    return linearInterpolationExtrap(undersat_oil_tables_[0][0],
						     undersat_oil_tables_[0][item],
						     maxR);
		}

		// Extrapolate from last table section
		int ltp = saturated_oil_table_[0].size() - 1;
		if (is+1 == ltp && press > saturated_oil_table_[0][ltp]) {
		    return linearInterpolationExtrap(undersat_oil_tables_[ltp][0],
						     undersat_oil_tables_[ltp][item],
						     maxR);
		}

		// Interpolate between table sections
		double w = (maxR - saturated_oil_table_[3][is]) /
		    (saturated_oil_table_[3][is+1] - 
		     saturated_oil_table_[3][is]);
		if (undersat_oil_tables_[is][0].size() < 2) {
		    double val = saturated_oil_table_[item][is] +
			w*(saturated_oil_table_[item][is+1] -
			   saturated_oil_table_[item][is]);
		    return val;
		}
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
