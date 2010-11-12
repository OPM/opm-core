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

#ifndef SINTEF_MISCIBILITYLIVEOIL_HEADER
#define SINTEF_MISCIBILITYLIVEOIL_HEADER

    /** Class for miscible live oil.
     *  Detailed description.
     */

#include "MiscibilityProps.hpp"

namespace Opm
{
    class MiscibilityLiveOil : public MiscibilityProps
    {
    public:
	typedef std::vector<std::vector<std::vector<double> > > table_t;

	MiscibilityLiveOil(const table_t& pvto);
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

