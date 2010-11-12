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

#ifndef SINTEF_MISCIBILITYPROPS_HEADER
#define SINTEF_MISCIBILITYPROPS_HEADER

    /** Base class for properties of fluids and rocks.
     *  Detailed description.
     */

#include <vector>
#include <tr1/array>

namespace Opm
{


    class MiscibilityProps
    {
    public:
	typedef std::tr1::array<double, 3> surfvol_t;
        enum PhaseNames { Aqua = 0, Liquid = 1, Vapour = 2 };

	MiscibilityProps();
	virtual ~MiscibilityProps();

        virtual double getViscosity(int region, double press, const surfvol_t& surfvol) const = 0;
        virtual double B   (int region, double press, const surfvol_t& surfvol) const = 0;
	virtual double dBdp(int region, double press, const surfvol_t& surfvol) const = 0;
	virtual double R   (int region, double press, const surfvol_t& surfvol) const = 0;
	virtual double dRdp(int region, double press, const surfvol_t& surfvol) const = 0;
    };

} // namespace Opm

#endif // SINTEF_MISCIBILITYPROPS_HEADER

