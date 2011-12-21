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

#ifndef OPENRS_MISCIBILITYWATER_HEADER
#define OPENRS_MISCIBILITYWATER_HEADER

#include <opm/core/fluid/blackoil/MiscibilityProps.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

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
            ref_press_ = pvtw[region_number][0];
            ref_B_     = pvtw[region_number][1];
            comp_      = pvtw[region_number][2];
            viscosity_ = pvtw[region_number][3];
            if (pvtw[region_number].size() > 4 && pvtw[region_number][4] != 0.0) {
                THROW("MiscibilityWater does not support 'viscosibility'.");
            }
        }

	MiscibilityWater(double visc)
            : ref_press_(0.0),
              ref_B_(1.0),
              comp_(0.0),
              viscosity_(visc)
        {
        }

	virtual ~MiscibilityWater()
        {
        }

        virtual double getViscosity(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return viscosity_;
        }

        virtual void getViscosity(const std::vector<PhaseVec>& pressures,
                                  const std::vector<CompVec>&,
                                  int,
                                  std::vector<double>& output) const
        {
            int num = pressures.size();
            output.clear();
            output.resize(num, viscosity_);
        }

        virtual double B(int /*region*/, double press, const surfvol_t& /*surfvol*/) const
        {
            if (comp_) {
                // Computing a polynomial approximation to the exponential.
                double x = comp_*(press - ref_press_);
                return ref_B_/(1.0 + x + 0.5*x*x);
            } else {
                return ref_B_;
            }
        }

        virtual void B(const std::vector<PhaseVec>& pressures,
                       const std::vector<CompVec>&,
                       int phase,
                       std::vector<double>& output) const
        {
            int num = pressures.size();
            if (comp_) {
                output.resize(num);
#pragma omp parallel for
                for (int i = 0; i < num; ++i) {
                    // Computing a polynomial approximation to the exponential.
                    double x = comp_*(pressures[i][phase] - ref_press_);
                    output[i] = ref_B_/(1.0 + x + 0.5*x*x);
                }
            } else {
                output.clear();
                output.resize(num, ref_B_);
            }
        }

	virtual double dBdp(int region, double press, const surfvol_t& surfvol) const
        {
            if (comp_) {
                return -comp_*B(region, press, surfvol);
            } else {
                return 0.0;
            }
        }

        virtual void dBdp(const std::vector<PhaseVec>& pressures,
                          const std::vector<CompVec>& surfvols,
                          int phase,
                          std::vector<double>& output_B,
                          std::vector<double>& output_dBdp) const
        {
            B(pressures, surfvols, phase, output_B);
            int num = pressures.size();
            if (comp_) {
                output_dBdp.resize(num);
#pragma omp parallel for
                for (int i = 0; i < num; ++i) {
                    output_dBdp[i] = -comp_*output_B[i];
                }
            } else {
                output_dBdp.clear();
                output_dBdp.resize(num, 0.0);
            }
        }

	virtual double R(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }

        virtual void R(const std::vector<PhaseVec>& pressures,
                       const std::vector<CompVec>&,
                       int,
                       std::vector<double>& output) const
        {
            int num = pressures.size();
            output.clear();
            output.resize(num, 0.0);
        }

	virtual double dRdp(int /*region*/, double /*press*/, const surfvol_t& /*surfvol*/) const
        {
            return 0.0;
        }

        virtual void dRdp(const std::vector<PhaseVec>& pressures,
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

    private:
        double ref_press_;
        double ref_B_;
        double comp_;
        double viscosity_;
    };

}

#endif // OPENRS_MISCIBILITYWATER_HEADER
