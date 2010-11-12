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

#ifndef SINTEF_MISCIBILITYDEAD_HEADER
#define SINTEF_MISCIBILITYDEAD_HEADER

/** Class for immiscible dead oil and dry gas.
 *  Detailed description.
 */

#include "MiscibilityProps.hpp"

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
	table_t pvdx_;
    };

}

#endif // SINTEF_MISCIBILITYDEAD_HEADER

