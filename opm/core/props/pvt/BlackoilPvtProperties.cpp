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



#include "config.h"
#include <opm/core/props/pvt/BlackoilPvtProperties.hpp>
#include <opm/core/props/pvt/SinglePvtDead.hpp>
#include <opm/core/props/pvt/SinglePvtDeadSpline.hpp>
#include <opm/core/props/pvt/SinglePvtLiveOil.hpp>
#include <opm/core/props/pvt/SinglePvtLiveGas.hpp>
#include <opm/core/props/pvt/SinglePvtConstCompr.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>


namespace Opm
{

    BlackoilPvtProperties::BlackoilPvtProperties()
    {
    }


    void BlackoilPvtProperties::init(const EclipseGridParser& deck, const int samples)
    {
        typedef std::vector<std::vector<std::vector<double> > > table_t;
        // If we need multiple regions, this class and the SinglePvt* classes must change.
        region_number_ = 0;

        phase_usage_ = phaseUsageFromDeck(deck);

        // Surface densities. Accounting for different orders in eclipse and our code.
        if (deck.hasField("DENSITY")) {
            const std::vector<double>& d = deck.getDENSITY().densities_[region_number_];
            enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
            if (phase_usage_.phase_used[Aqua]) {
                densities_[phase_usage_.phase_pos[Aqua]]   = d[ECL_water];
            }
            if (phase_usage_.phase_used[Vapour]) {
                densities_[phase_usage_.phase_pos[Vapour]] = d[ECL_gas];
            }
            if (phase_usage_.phase_used[Liquid]) {
                densities_[phase_usage_.phase_pos[Liquid]] = d[ECL_oil];
            }
        } else {
            THROW("Input is missing DENSITY\n");
        }

        // Set the properties.
        props_.resize(phase_usage_.num_phases);
        // Water PVT
        if (phase_usage_.phase_used[Aqua]) {
            if (deck.hasField("PVTW")) {
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(deck.getPVTW().pvtw_));
            } else {
                // Eclipse 100 default.
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(0.5*Opm::prefix::centi*Opm::unit::Poise));
            }
        }
        // Oil PVT
        if (phase_usage_.phase_used[Liquid]) {
            if (deck.hasField("PVDO")) {
                if (samples > 0) {
                    props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDeadSpline(deck.getPVDO().pvdo_, samples));
                } else {
                    props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDead(deck.getPVDO().pvdo_));
                }
            } else if (deck.hasField("PVTO")) {
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtLiveOil(deck.getPVTO().pvto_));
            } else if (deck.hasField("PVCDO")) {
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtConstCompr(deck.getPVCDO().pvcdo_));
            } else {
                THROW("Input is missing PVDO or PVTO\n");
            }
        }
        // Gas PVT
        if (phase_usage_.phase_used[Vapour]) {
            if (deck.hasField("PVDG")) {
                if (samples > 0) {
                    props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDeadSpline(deck.getPVDG().pvdg_, samples));
                } else {
                    props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDead(deck.getPVDG().pvdg_));
                }
            } else if (deck.hasField("PVTG")) {
                props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtLiveGas(deck.getPVTG().pvtg_));
            } else {
                THROW("Input is missing PVDG or PVTG\n");
            }
        }

        // Must inform pvt property objects of phase structure.
        for (int i = 0; i < phase_usage_.num_phases; ++i) {
            props_[i]->setPhaseConfiguration(phase_usage_.num_phases, phase_usage_.phase_pos);
        }
    }

    const double* BlackoilPvtProperties::surfaceDensities() const
    {
        return densities_;
    }


    int BlackoilPvtProperties::numPhases() const
    {
        return phase_usage_.num_phases;
    }

    const int* BlackoilPvtProperties::phaseUsed() const
    {
        return phase_usage_.phase_used;
    }

    const int* BlackoilPvtProperties::phasePosition() const
    {
        return phase_usage_.phase_pos;
    }


    void BlackoilPvtProperties::mu(const int n,
                                   const double* p,
                                   const double* z,
                                   double* output_mu) const
    {
        data1_.resize(n);
        for (int phase = 0; phase < phase_usage_.num_phases; ++phase) {
            props_[phase]->mu(n, p, z, &data1_[0]);
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_mu[phase_usage_.num_phases*i + phase] = data1_[i];
            }
        }
    }

    void BlackoilPvtProperties::B(const int n,
                                  const double* p,
                                  const double* z,
                                  double* output_B) const
    {
        data1_.resize(n);
        for (int phase = 0; phase < phase_usage_.num_phases; ++phase) {
            props_[phase]->B(n, p, z, &data1_[0]);
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_B[phase_usage_.num_phases*i + phase] = data1_[i];
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
        for (int phase = 0; phase < phase_usage_.num_phases; ++phase) {
            props_[phase]->dBdp(n, p, z, &data1_[0], &data2_[0]);
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_B[phase_usage_.num_phases*i + phase] = data1_[i];
                output_dBdp[phase_usage_.num_phases*i + phase] = data2_[i];
            }
        }
    }


    void BlackoilPvtProperties::R(const int n,
                                  const double* p,
                                  const double* z,
                                  double* output_R) const
    {
        data1_.resize(n);
        for (int phase = 0; phase < phase_usage_.num_phases; ++phase) {
            props_[phase]->R(n, p, z, &data1_[0]);
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_R[phase_usage_.num_phases*i + phase] = data1_[i];
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
        for (int phase = 0; phase < phase_usage_.num_phases; ++phase) {
            props_[phase]->dRdp(n, p, z, &data1_[0], &data2_[0]);
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_R[phase_usage_.num_phases*i + phase] = data1_[i];
                output_dRdp[phase_usage_.num_phases*i + phase] = data2_[i];
            }
        }
    }

} // namespace Opm
