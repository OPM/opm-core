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
#include <opm/core/props/pvt/PvtDead.hpp>
#include <opm/core/props/pvt/PvtDeadSpline.hpp>
#include <opm/core/props/pvt/PvtLiveOil.hpp>
#include <opm/core/props/pvt/PvtLiveGas.hpp>
#include <opm/core/props/pvt/PvtConstCompr.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>

#include <opm/parser/eclipse/Utility/PvtwTable.hpp>
#include <opm/parser/eclipse/Utility/PvcdoTable.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

namespace Opm
{

    BlackoilPvtProperties::BlackoilPvtProperties()
    {
    }

    void BlackoilPvtProperties::init(Opm::DeckConstPtr deck,
                                     int numSamples,
                                     int numCompressedCells,
                                     const int *compressedToCartesianCellIdx)
    {
        phase_usage_ = phaseUsageFromDeck(deck);

        // Surface densities. Accounting for different orders in eclipse and our code.
        Opm::DeckKeywordConstPtr densityKeyword = deck->getKeyword("DENSITY");
        int numRegions = densityKeyword->size();

        densities_.resize(numRegions);
        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            if (phase_usage_.phase_used[Liquid]) {
                densities_[regionIdx][phase_usage_.phase_pos[Liquid]]
                    = densityKeyword->getRecord(regionIdx)->getItem("OIL")->getSIDouble(0);
            }
            if (phase_usage_.phase_used[Aqua]) {
                densities_[regionIdx][phase_usage_.phase_pos[Aqua]]
                    = densityKeyword->getRecord(regionIdx)->getItem("WATER")->getSIDouble(0);
            }
            if (phase_usage_.phase_used[Vapour]) {
                densities_[regionIdx][phase_usage_.phase_pos[Vapour]]
                    = densityKeyword->getRecord(regionIdx)->getItem("GAS")->getSIDouble(0);
            }
        }

        // first, calculate the PVT table index for each compressed
        // cell. This array is required to construct the PVT classes
        // below.
        Opm::extractPvtTableIndex(pvtTableIdx_,
                                  deck,
                                  numCompressedCells,
                                  compressedToCartesianCellIdx);

        // Resize the property objects container
        props_.resize(phase_usage_.num_phases);

        // Water PVT
        if (phase_usage_.phase_used[Aqua]) {
            // if water is used, we require the presence of the "PVTW"
            // keyword for now...
            std::shared_ptr<PvtConstCompr> pvtw(new PvtConstCompr);
            pvtw->initFromWater(deck->getKeyword("PVTW"), pvtTableIdx_);

            props_[phase_usage_.phase_pos[Aqua]] = pvtw;
        }
        // Oil PVT
        if (phase_usage_.phase_used[Liquid]) {
            // for oil, we support the "PVDO", "PVTO" and "PVCDO"
            // keywords...
            if (deck->hasKeyword("PVDO")) {
                Opm::DeckKeywordConstPtr pvdoKeyword = deck->getKeyword("PVDO");
                if (numSamples > 0) {
                    auto splinePvt = std::shared_ptr<PvtDeadSpline>(new PvtDeadSpline);
                    splinePvt->initFromOil(pvdoKeyword, pvtTableIdx_, numSamples);
                    props_[phase_usage_.phase_pos[Liquid]] = splinePvt;
                } else {
                    auto deadPvt = std::shared_ptr<PvtDead>(new PvtDead);
                    deadPvt->initFromOil(pvdoKeyword, pvtTableIdx_);
                    props_[phase_usage_.phase_pos[Liquid]] = deadPvt;
                }
            } else if (deck->hasKeyword("PVTO")) {
                props_[phase_usage_.phase_pos[Liquid]].reset(new PvtLiveOil(deck->getKeyword("PVTO"), pvtTableIdx_));
            } else if (deck->hasKeyword("PVCDO")) {
                std::shared_ptr<PvtConstCompr> pvcdo(new PvtConstCompr);
                pvcdo->initFromOil(deck->getKeyword("PVCDO"), pvtTableIdx_);

                props_[phase_usage_.phase_pos[Liquid]] = pvcdo;
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDO, PVCDO or PVTO\n");
            }
        }
        // Gas PVT
        if (phase_usage_.phase_used[Vapour]) {
            // gas can be specified using the "PVDG" or "PVTG" keywords...
            if (deck->hasKeyword("PVDG")) {
                Opm::DeckKeywordConstPtr pvdgKeyword = deck->getKeyword("PVDG");

                if (numSamples > 0) {
                    std::shared_ptr<PvtDeadSpline> splinePvt(new PvtDeadSpline);
                    splinePvt->initFromGas(pvdgKeyword, pvtTableIdx_, numSamples);

                    props_[phase_usage_.phase_pos[Vapour]] = splinePvt;
                } else {
                    std::shared_ptr<PvtDead> deadPvt(new PvtDead);
                    deadPvt->initFromGas(pvdgKeyword, pvtTableIdx_);

                    props_[phase_usage_.phase_pos[Vapour]] = deadPvt;
                }
            } else if (deck->hasKeyword("PVTG")) {
                props_[phase_usage_.phase_pos[Vapour]].reset(new PvtLiveGas(deck->getKeyword("PVTG"), pvtTableIdx_));
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDG or PVTG\n");
            }
        }
    }

    const double* BlackoilPvtProperties::surfaceDensities(int regionIdx) const
    {
        return &densities_[regionIdx][0];
    }


    PhaseUsage BlackoilPvtProperties::phaseUsage() const
    {
        return phase_usage_;
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
