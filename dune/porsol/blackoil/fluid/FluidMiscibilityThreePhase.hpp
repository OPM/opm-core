//===========================================================================
//                                                                           
// File: FluidMiscibilityThreePhase.hpp                                      
//                                                                           
// Created: Thu Feb 11 10:35:05 2010                                         
//                                                                           
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
//===========================================================================

#ifndef SINTEF_FLUIDMISCIBILITYTHREEPHASE_HEADER
#define SINTEF_FLUIDMISCIBILITYTHREEPHASE_HEADER

    /** Temporary class for testing Miscibility* classes
     *  Detailed description.
     */
#include <string>
#include <boost/scoped_ptr.hpp>
#include "MiscibilityProps.hpp"

namespace Opm
{
    class FluidMiscibilityThreePhase
    {
    public:
        typedef MiscibilityProps::surfvol_t surfvol_t;

        enum PhaseNames { Aqua = 0, Liquid = 1, Vapour = 2 };

	FluidMiscibilityThreePhase(){}
	~FluidMiscibilityThreePhase(){}
	void init(const std::string& pvt_filename);

        double getViscosity(double press, const surfvol_t& surfvol,
			    PhaseNames phase) const;
        surfvol_t getMobilities(double press, const surfvol_t& sat, const surfvol_t& surfvol) const;
	surfvol_t surfaceDensities() const;
        double B   (double press, const surfvol_t& surfvol,
		    PhaseNames phase) const;
        double dBdp(double press, const surfvol_t& surfvol,
		    PhaseNames phase) const;
        double R   (double press, const surfvol_t& surfvol,
		    PhaseNames phase) const;
        double dRdp(double press, const surfvol_t& surfvol,
		    PhaseNames phase) const;


    private:
	int region_number_;
        const MiscibilityProps& propsForPhase(PhaseNames phase) const;
	double kro(double sw, double sg) const; // Stone's first model(Modified)
	double dkro_dsw(double sw, double sg) const; // kro partial derivative dSw
	double dkro_dsg(double sw, double sg) const; // kro partial derivative dSg
	void   dkro (double sw, double sg,   // kro partial derivatives. dSw and dSg
		     double& dkro_dsw, double& dkro_dsg) const;
        double krw  (double sw) const;  // Water relative permeability
        double krow (double sw) const;  // Oil relative permeability
        double dkrow_dsw (double sw) const; // krow derivative
        double Pcow (double sw) const;  // Water-oil capillary pressure
        double krg  (double sg) const;  // Gas relative permeability
        double krog (double sg) const;  // Oil relative permeability
        double dkrog_dsg(double sg) const; // krog derivative
        double Pcog (double sg) const;  // Oil-gas capillary pressure

	boost::scoped_ptr<MiscibilityProps> water_props_;
	boost::scoped_ptr<MiscibilityProps> oil_props_;
	boost::scoped_ptr<MiscibilityProps> gas_props_;
	surfvol_t densities_;
	std::vector<std::vector<std::vector<double> > > sgof_;  // Gas/Oil saturation
	std::vector<std::vector<std::vector<double> > > swof_;  // Water/Oil saturation
    };

}

#endif // SINTEF_FLUIDMISCIBILITYTHREEPHASE_HEADER

