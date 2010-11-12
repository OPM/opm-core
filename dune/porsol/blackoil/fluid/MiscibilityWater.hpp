//===========================================================================
//
// File: MiscibilityWater.hpp
//
// Created: Tue May 18 10:26:13 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bjørn Spjelkavik    <bsp@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

#ifndef OPENRS_MISCIBILITYWATER_HEADER
#define OPENRS_MISCIBILITYWATER_HEADER

#include "MiscibilityProps.hpp"
#include <dune/common/ErrorMacros.hpp>

// Forward declaration.
class PVTW;

namespace Opm
{
    class MiscibilityWater : public MiscibilityProps
    {
    public:
        typedef std::vector<std::vector<double> > table_t;
	MiscibilityWater(const table_t& pvtw)
        {
	    const int region_number = 0;
	    if (pvtw.size() != 1) {
		THROW("More than one PVD-region");
	    }
            double b = pvtw[region_number][1];
            double comp = pvtw[region_number][2];
            double visc = pvtw[region_number][3];
            const double VISCOSITY_UNIT = 1e-3;
            if (b == 1.0 && comp == 0) {
                viscosity_ = visc*VISCOSITY_UNIT;
            } else {
                THROW("Compressible water not implemented.");
            }
        }
	MiscibilityWater(double visc)
            : viscosity_(visc)
        {
        }
	virtual ~MiscibilityWater()
        {
        }

        virtual double getViscosity(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return viscosity_;
        }
        virtual double B(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 1.0;
        }
	virtual double dBdp(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }
	virtual double R(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }
	virtual double dRdp(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }
    private:
        double viscosity_;
    };

}

#endif // OPENRS_MISCIBILITYWATER_HEADER
