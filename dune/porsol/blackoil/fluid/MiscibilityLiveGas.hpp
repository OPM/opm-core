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

