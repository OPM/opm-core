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

#ifndef OPM_PVTLIVEGAS_HEADER_INCLUDED
#define OPM_PVTLIVEGAS_HEADER_INCLUDED

#include <opm/core/props/pvt/PvtInterface.hpp>

#include <opm/parser/eclipse/Utility/PvtgTable.hpp>

#include <vector>

namespace Opm
{
    /// Class for miscible wet gas (with vaporized oil in vapour phase).
    /// The PVT properties can either be given as a function of pressure (p) and surface volume (z)
    /// or pressure (p) and gas resolution factor (r).
    /// For all the virtual methods, the following apply: p, r and z
    /// are expected to be of size n, size n and n*num_phases, respectively.
    /// Output arrays shall be of size n, and must be valid before
    /// calling the method.
    class PvtLiveGas : public PvtInterface
    {
    public:
        PvtLiveGas(Opm::DeckKeywordConstPtr pvtgKeyword);
        virtual ~PvtLiveGas();

        /// Viscosity as a function of p and z.
        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* z,
                        double* output_mu) const;

        /// Viscosity and its derivatives as a function of p and r.
        /// The fluid is considered saturated if r >= rsSat(p).
        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* r,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const;

        /// Viscosity as a function of p and r.
        /// State condition determined by 'cond'.
        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* r,
                        const PhasePresence* cond,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const;

        /// Formation volume factor as a function of p and z.
        virtual void B(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* z,
                       double* output_B) const;

        /// Formation volume factor and p-derivative as functions of p and z.
        virtual void dBdp(const int n,
                          const int* pvtRegionIdx,
                          const double* p,
                          const double* z,
                          double* output_B,
                          double* output_dBdp) const;

        /// The inverse of the formation volume factor b = 1 / B, and its derivatives as a function of p and r.
        /// The fluid is considered saturated if r >= rsSat(p).
        virtual void b(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* r,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const;

        /// The inverse of the formation volume factor b = 1 / B, and its derivatives as a function of p and r.
        /// State condition determined by 'cond'.
        virtual void b(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* r,
                       const PhasePresence* cond,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const;



        /// Solution gas/oil ratio and its derivatives at saturated conditions as a function of p.
        virtual void rsSat(const int n,
                           const int* pvtRegionIdx,
                          const double* p,
                          double* output_rsSat,
                          double* output_drsSatdp) const;

        /// Vapor oil/gas ratio and its derivatives at saturated conditions as a function of p.
        virtual void rvSat(const int n,
                           const int* pvtRegionIdx,
                          const double* p,
                          double* output_rvSat,
                          double* output_drvSatdp) const;

        /// Solution factor as a function of p and z.
        virtual void R(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* z,
                       double* output_R) const;

        /// Solution factor and p-derivative as functions of p and z.
        virtual void dRdp(const int n,
                          const int* pvtRegionIdx,
                          const double* p,
                          const double* z,
                          double* output_R,
                          double* output_dRdp) const;

    protected:
        int getTableIndex_(const int* pvtTableIdx, int cellIdx) const
        {
            if (!pvtTableIdx)
                return 0;
            return pvtTableIdx[cellIdx];
        }

        double evalB(double press, const double* surfvol, int pvtTableIdx) const;
        void evalBDeriv(double press, const double* surfvol, int pvtTableIdx, double& B, double& dBdp) const;
        double evalR(double press, const double* surfvol, int pvtTableIdx) const;
        void evalRDeriv(double press, const double* surfvol, int pvtTableIdx, double& R, double& dRdp) const;

        // item:  1=>B  2=>mu;
        double miscible_gas(const double press,
                            const double* surfvol,
                            const int pvtTableIdx,
                            const int item,
                            const bool deriv = false) const;
        double miscible_gas(const double press,
                            const double r,
                            const PhasePresence& cond,
                            const int pvtTableIdx,
                            const int item,
                            const int deriv = 0) const;
        // PVT properties of wet gas (with vaporised oil). We need to
        // store one table per PVT region.
        std::vector< std::vector<std::vector<double> > > saturated_gas_table_;
        std::vector< std::vector<std::vector<std::vector<double> > > > undersat_gas_tables_;
    };

}

#endif // OPM_PVTLIVEGAS_HEADER_INCLUDED

