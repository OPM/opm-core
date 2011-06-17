//===========================================================================
//                                                                           
// File: MiscibilityDead.cpp                                                  
//                                                                           
// Created: Wed Feb 10 09:06:04 2010                                         
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
#include "MiscibilityDead.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/linInt.hpp>
#include <dune/porsol/common/buildUniformMonotoneTable.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

using namespace std;
using namespace Dune;

namespace Opm
{

    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    MiscibilityDead::MiscibilityDead(const table_t& pvd_table)
    {
	const int region_number = 0;
	if (pvd_table.size() != 1) {
	    THROW("More than one PVT-region");
	}

	// Copy data
	const int sz = pvd_table[region_number][0].size();
        std::vector<double> press(sz);
        std::vector<double> B_inv(sz);
        std::vector<double> visc(sz);
	for (int i = 0; i < sz; ++i) {
            press[i] = pvd_table[region_number][0][i];
            B_inv[i] = 1.0 / pvd_table[region_number][1][i];
            visc[i]  = pvd_table[region_number][2][i];
	}
        int samples = 1025;
        buildUniformMonotoneTable(press, B_inv, samples, one_over_B_);
        buildUniformMonotoneTable(press, visc, samples, viscosity_);

        // Dumping the created tables.
//         static int count = 0;
//         std::ofstream os((std::string("dump-") + boost::lexical_cast<std::string>(count++)).c_str());
//         os.precision(15);
//         os << "1/B\n\n" << one_over_B_
//            << "\n\nvisc\n\n" << viscosity_ << std::endl;
    }

    // Destructor
    MiscibilityDead::~MiscibilityDead()
    {
    }

    double MiscibilityDead::getViscosity(int region, double press, const surfvol_t& /*surfvol*/) const
    {
	return viscosity_(press);
    }

    void MiscibilityDead::getViscosity(const std::vector<PhaseVec>& pressures,
                                       const std::vector<CompVec>&,
                                       int phase,
                                       std::vector<double>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i] = viscosity_(pressures[i][phase]);
        }
    }

    double MiscibilityDead::B(int region, double press, const surfvol_t& /*surfvol*/) const
    {
	// Interpolate 1/B 
	return 1.0/one_over_B_(press);
    }

    void MiscibilityDead::B(const std::vector<PhaseVec>& pressures,
                            const std::vector<CompVec>&,
                            int phase,
                            std::vector<double>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i] = 1.0/one_over_B_(pressures[i][phase]);
        }
    }

    double MiscibilityDead::dBdp(int region, double press, const surfvol_t& /*surfvol*/) const
    {
	// Interpolate 1/B
	surfvol_t dummy_surfvol;
	double Bg = B(region, press, dummy_surfvol);
	return -Bg*Bg*one_over_B_.derivative(press);
    }

    void MiscibilityDead::dBdp(const std::vector<PhaseVec>& pressures,
                               const std::vector<CompVec>& surfvols,
                               int phase,
                               std::vector<double>& output_B,
                               std::vector<double>& output_dBdp) const
    {
        B(pressures, surfvols, phase, output_B);
        int num = pressures.size();
        output_dBdp.resize(num);
        for (int i = 0; i < num; ++i) {
            double Bg = output_B[i];
            output_dBdp[i] = -Bg*Bg*one_over_B_.derivative(pressures[i][phase]);
        }
    }

    double MiscibilityDead::R(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
    {
        return 0.0;
    }

    void MiscibilityDead::R(const std::vector<PhaseVec>& pressures,
                            const std::vector<CompVec>&,
                            int,
                            std::vector<double>& output) const
    {
        int num = pressures.size();
        output.clear();
        output.resize(num, 0.0);
    }

    double MiscibilityDead::dRdp(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
    {
        return 0.0;
    }

    void MiscibilityDead::dRdp(const std::vector<PhaseVec>& pressures,
                               const std::vector<CompVec>&,
                               int,
                               std::vector<double>& output_R,
                               std::vector<double>& output_dRdp) const
    {
        int num = pressures.size();
        output_R.clear();
        output_R.resize(num, 0.0);
        output_dRdp.clear();
        output_dRdp.resize(num, 0.0);
    }

}
