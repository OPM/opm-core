/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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



#include <opm/core/fluid/blackoil/BlackoilPvtProperties.hpp>
#include <opm/core/fluid/blackoil/SinglePvtDead.hpp>
#include <opm/core/fluid/blackoil/SinglePvtLiveOil.hpp>
#include <opm/core/fluid/blackoil/SinglePvtLiveGas.hpp>
#include <opm/core/fluid/blackoil/SinglePvtConstCompr.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linInt.hpp>


namespace Opm
{

    void BlackoilPvtProperties::init(const Dune::EclipseGridParser& deck)
    {
        typedef std::vector<std::vector<std::vector<double> > > table_t;
        // If we need multiple regions, this class and the SinglePvt* classes must change.
	region_number_ = 0;

        // Discover phase usage.
        std::fill(phase_used_, phase_used_ + MaxNumPhases, 0);
        if (deck.hasField("WATER")) {
            phase_used_[Aqua] = 1;
        }
        if (deck.hasField("OIL")) {
            phase_used_[Liquid] = 1;
        }
        if (deck.hasField("GAS")) {
            phase_used_[Vapour] = 1;
        }
        num_phases_ = 0;
        for (int i = 0; i < MaxNumPhases; ++i) {
            phase_pos_[i] = num_phases_;
            num_phases_ += phase_used_[i];
        }

        // Only 2 or 3 phase systems handled.
        if (num_phases_ < 2 || num_phases_ > 3) {
            THROW("Cannot handle cases with " << num_phases_ << " phases.");
        }

        // We need oil systems, since we do not support the keywords needed for
        // water-gas systems.
        if (!phase_used_[Liquid]) {
            THROW("Cannot handle cases with no OIL, i.e. water-gas systems.");
        }

	// Surface densities. Accounting for different orders in eclipse and our code.
	if (deck.hasField("DENSITY")) {
	    const std::vector<double>& d = deck.getDENSITY().densities_[region_number_];
	    enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
	    densities_[Aqua]   = d[ECL_water];
	    densities_[Vapour] = d[ECL_gas];
	    densities_[Liquid] = d[ECL_oil];
	} else {
	    THROW("Input is missing DENSITY\n");
	}

        // Set the properties.
        props_.resize(num_phases_);
        // Water PVT
        if (phase_used_[Aqua]) {
            if (deck.hasField("PVTW")) {
                props_[phase_pos_[Aqua]].reset(new SinglePvtConstCompr(deck.getPVTW().pvtw_));
            } else {
                // Eclipse 100 default.
                props_[phase_pos_[Aqua]].reset(new SinglePvtConstCompr(0.5*Dune::prefix::centi*Dune::unit::Poise));
            }
        }
        // Oil PVT
        if (phase_used_[Liquid]) {
            if (deck.hasField("PVDO")) {
                props_[phase_pos_[Liquid]].reset(new SinglePvtDead(deck.getPVDO().pvdo_));
            } else if (deck.hasField("PVTO")) {
                props_[phase_pos_[Liquid]].reset(new SinglePvtLiveOil(deck.getPVTO().pvto_));
            } else if (deck.hasField("PVCDO")) {
                props_[phase_pos_[Liquid]].reset(new SinglePvtConstCompr(deck.getPVCDO().pvcdo_));
            } else {
                THROW("Input is missing PVDO or PVTO\n");
            }
        }
	// Gas PVT
        if (phase_used_[Vapour]) {
            if (deck.hasField("PVDG")) {
                props_[phase_pos_[Vapour]].reset(new SinglePvtDead(deck.getPVDG().pvdg_));
            } else if (deck.hasField("PVTG")) {
                props_[phase_pos_[Vapour]].reset(new SinglePvtLiveGas(deck.getPVTG().pvtg_));
            } else {
                THROW("Input is missing PVDG or PVTG\n");
            }
        }

        // Must inform pvt property objects of phase structure.
        for (int i = 0; i < num_phases_; ++i) {
            props_[i]->setPhaseConfiguration(num_phases_, phase_pos_);
        }
    }

    const double* BlackoilPvtProperties::surfaceDensities() const
    {
        return densities_;
    }


    void BlackoilPvtProperties::mu(const int n,
                                   const double* p,
                                   const double* z,
                                   double* output_mu) const
    {
        data1_.resize(n);
        for (int phase = 0; phase < num_phases_; ++phase) {
            props_[phase]->mu(n, p, z, &data1_[0]);
#pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_mu[num_phases_*i + phase] = data1_[i];
            }
        }
    }

    void BlackoilPvtProperties::B(const int n,
                                  const double* p,
                                  const double* z,
                                  double* output_B) const
    {
        data1_.resize(n);
        for (int phase = 0; phase < num_phases_; ++phase) {
            props_[phase]->B(n, p, z, &data1_[0]);
#pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_B[num_phases_*i + phase] = data1_[i];
            }
        }
    }

    void BlackoilPvtProperties::dBdp(const int n,
                                     const double* p,
                                     const double* z,
                                     double* output_B,
                                     double* output_dBdp) const
    {
        data1_.resize(n);
        data2_.resize(n);
        for (int phase = 0; phase < num_phases_; ++phase) {
            props_[phase]->dBdp(n, p, z, &data1_[0], &data2_[0]);
#pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_B[num_phases_*i + phase] = data1_[i];
                output_dBdp[num_phases_*i + phase] = data2_[i];
            }
        }
    }


    void BlackoilPvtProperties::R(const int n,
                                  const double* p,
                                  const double* z,
                                  double* output_R) const
    {
        data1_.resize(n);
        for (int phase = 0; phase < num_phases_; ++phase) {
            props_[phase]->R(n, p, z, &data1_[0]);
#pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_R[num_phases_*i + phase] = data1_[i];
            }
        }
    }

    void BlackoilPvtProperties::dRdp(const int n,
                                     const double* p,
                                     const double* z,
                                     double* output_R,
                                     double* output_dRdp) const
    {
        data1_.resize(n);
        data2_.resize(n);
        for (int phase = 0; phase < num_phases_; ++phase) {
            props_[phase]->dRdp(n, p, z, &data1_[0], &data2_[0]);
#pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_R[num_phases_*i + phase] = data1_[i];
                output_dRdp[num_phases_*i + phase] = data2_[i];
            }
        }
    }

} // namespace Opm
