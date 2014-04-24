/*
  Copyright 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_PVTCONSTCOMPR_HEADER_INCLUDED
#define OPM_PVTCONSTCOMPR_HEADER_INCLUDED

#include <opm/core/props/pvt/PvtInterface.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/parser/eclipse/Utility/PvtwTable.hpp>
#include <opm/parser/eclipse/Utility/PvcdoTable.hpp>

#include <vector>
#include <algorithm>


namespace Opm
{

    /// Class for constant compressible phases (PVTW or PVCDO).
    /// The PVT properties can either be given as a function of pressure (p) and surface volume (z)
    /// or pressure (p) and gas resolution factor (r).
    /// For all the virtual methods, the following apply: p, r and z
    /// are expected to be of size n, size n and n*num_phases, respectively.
    /// Output arrays shall be of size n, and must be valid before
    /// calling the method.
    class PvtConstCompr : public PvtInterface
    {
    public:
        PvtConstCompr()
        {}

        void initFromWater(Opm::DeckKeywordConstPtr pvtwKeyword, const std::vector<int> &pvtTableIdx)
        {
            pvtTableIdx_ = pvtTableIdx;

            int numRegions = pvtwKeyword->size();

            ref_press_.resize(numRegions);
            ref_B_.resize(numRegions);
            comp_.resize(numRegions);
            viscosity_.resize(numRegions);
            visc_comp_.resize(numRegions);

            for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                Opm::DeckRecordConstPtr pvtwRecord = pvtwKeyword->getRecord(regionIdx);

                ref_press_[regionIdx] = pvtwRecord->getItem("P_REF")->getSIDouble(0);
                ref_B_[regionIdx]     = pvtwRecord->getItem("WATER_VOL_FACTOR")->getSIDouble(0);
                comp_[regionIdx]      = pvtwRecord->getItem("WATER_COMPRESSIBILITY")->getSIDouble(0);
                viscosity_[regionIdx] = pvtwRecord->getItem("WATER_VISCOSITY")->getSIDouble(0);
                visc_comp_[regionIdx] = pvtwRecord->getItem("WATER_VISCOSIBILITY")->getSIDouble(0);
            }
        }

        void initFromOil(Opm::DeckKeywordConstPtr pvcdoKeyword, const std::vector<int> &pvtTableIdx)
        {
            pvtTableIdx_ = pvtTableIdx;

            int numRegions = pvcdoKeyword->size();

            ref_press_.resize(numRegions);
            ref_B_.resize(numRegions);
            comp_.resize(numRegions);
            viscosity_.resize(numRegions);
            visc_comp_.resize(numRegions);

            for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                Opm::DeckRecordConstPtr pvcdoRecord = pvcdoKeyword->getRecord(regionIdx);

                ref_press_[regionIdx] = pvcdoRecord->getItem("P_REF")->getSIDouble(0);
                ref_B_[regionIdx]     = pvcdoRecord->getItem("BO_REF")->getSIDouble(0);
                comp_[regionIdx]      = pvcdoRecord->getItem("C_OIL_REF")->getSIDouble(0);
                viscosity_[regionIdx] = pvcdoRecord->getItem("MUO_REF")->getSIDouble(0);
                visc_comp_[regionIdx] = pvcdoRecord->getItem("OIL_VISCOSIBILITY")->getSIDouble(0);
            }
        }

        PvtConstCompr(double visc)
            : ref_press_(1, 0.0),
              ref_B_(1, 1.0),
              comp_(1, 0.0),
              viscosity_(1, visc),
              visc_comp_(1, 0.0)
        {
        }

        virtual ~PvtConstCompr()
        {
        }

        virtual void mu(const int n,
                        const double* p,
                        const double* /*z*/,
                        double* output_mu) const
        {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                // Computing a polynomial approximation to the exponential.
                int tableIdx = getTableIndex_(i);
                double x = -visc_comp_[tableIdx]*(p[i] - ref_press_[tableIdx]);
                output_mu[i] = viscosity_[tableIdx]/(1.0 + x + 0.5*x*x);
            }
        }

        virtual void mu(const int n,
                        const double* p,
                        const double* /*r*/,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const
        {
            // #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                // Computing a polynomial approximation to the exponential.
                int tableIdx = getTableIndex_(i);
                double x = -visc_comp_[tableIdx]*(p[i] - ref_press_[tableIdx]);
                double d = (1.0 + x + 0.5*x*x);
                output_mu[i] = viscosity_[tableIdx]/d;
                output_dmudp[i] = (viscosity_[tableIdx]/(d*d))*(1+x) * visc_comp_[tableIdx];
            }
            std::fill(output_dmudr, output_dmudr + n, 0.0);
        }

        virtual void mu(const int n,
                        const double* p,
                        const double* /*r*/,
                        const PhasePresence* /*cond*/,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const
        {
            // #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                // Computing a polynomial approximation to the exponential.
                int tableIdx = getTableIndex_(i);
                double x = -visc_comp_[tableIdx]*(p[i] - ref_press_[tableIdx]);
                double d = (1.0 + x + 0.5*x*x);
                output_mu[i] = viscosity_[tableIdx]/d;
                output_dmudp[i] = (viscosity_[tableIdx]/(d*d))*(1+x) * visc_comp_[tableIdx];
            }
            std::fill(output_dmudr, output_dmudr + n, 0.0);
        }

        virtual void B(const int n,
                       const double* p,
                       const double* /*z*/,
                       double* output_B) const
        {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                // Computing a polynomial approximation to the exponential.
                int tableIdx = getTableIndex_(i);
                double x = comp_[tableIdx]*(p[i] - ref_press_[tableIdx]);
                output_B[i] = ref_B_[tableIdx]/(1.0 + x + 0.5*x*x);
            }
        }

        virtual void dBdp(const int n,
                          const double* p,
                          const double* /*z*/,
                          double* output_B,
                          double* output_dBdp) const
        {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                int tableIdx = getTableIndex_(i);
                double x = comp_[tableIdx]*(p[i] - ref_press_[tableIdx]);
                double d = (1.0 + x + 0.5*x*x);
                output_B[i] = ref_B_[tableIdx]/d;
                output_dBdp[i] = (-ref_B_[tableIdx]/(d*d))*(1 + x) * comp_[tableIdx];
            }
        }

        virtual void b(const int n,
                       const double* p,
                       const double* /*r*/,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const
        {
        // #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                // Computing a polynomial approximation to the exponential.
                int tableIdx = getTableIndex_(i);
                double x = comp_[tableIdx]*(p[i] - ref_press_[tableIdx]);
                double d = (1.0 + x + 0.5*x*x);

                // b = 1/B = d/ref_B_B;
                output_b[i] = d/ref_B_[tableIdx];
                output_dbdp[i] = (1 + x) * comp_[tableIdx]/ref_B_[tableIdx];
            }

            std::fill(output_dbdr, output_dbdr + n, 0.0);
        }

        virtual void b(const int n,
                       const double* p,
                       const double* /*r*/,
                       const PhasePresence* /*cond*/,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const
        {
        // #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                // Computing a polynomial approximation to the exponential.
                int tableIdx = getTableIndex_(i);
                double x = comp_[tableIdx]*(p[i] - ref_press_[tableIdx]);
                double d = (1.0 + x + 0.5*x*x);

                // b = 1/B = d/ref_B_[tableIdx]B;
                output_b[i] = d/ref_B_[tableIdx];
                output_dbdp[i] = (1 + x) * comp_[tableIdx]/ref_B_[tableIdx];
            }
            std::fill(output_dbdr, output_dbdr + n, 0.0);
        }

        virtual void rsSat(const int n,
                           const double* /*p*/,
                           double* output_rsSat,
                           double* output_drsSatdp) const
        {
            std::fill(output_rsSat, output_rsSat + n, 0.0);
            std::fill(output_drsSatdp, output_drsSatdp + n, 0.0);
        }

        virtual void rvSat(const int n,
                           const double* /*p*/,
                           double* output_rvSat,
                           double* output_drvSatdp) const
        {
            std::fill(output_rvSat, output_rvSat + n, 0.0);
            std::fill(output_drvSatdp, output_drvSatdp + n, 0.0);
        }

        virtual void R(const int n,
                       const double* /*p*/,
                       const double* /*z*/,
                       double* output_R) const
        {
            std::fill(output_R, output_R + n, 0.0);
        }

        virtual void dRdp(const int n,
                          const double* /*p*/,
                          const double* /*z*/,
                          double* output_R,
                          double* output_dRdp) const
        {
            std::fill(output_R, output_R + n, 0.0);
            std::fill(output_dRdp, output_dRdp + n, 0.0);
        }

    private:
        int getTableIndex_(int cellIdx) const
        {
            if (pvtTableIdx_.empty())
                return 0;
            return pvtTableIdx_[cellIdx];
        }

        std::vector<int> pvtTableIdx_;
        std::vector<double> ref_press_;
        std::vector<double> ref_B_;
        std::vector<double> comp_;
        std::vector<double> viscosity_;
        std::vector<double> visc_comp_;
    };

}


#endif // OPM_PVTCONSTCOMPR_HEADER_INCLUDED

