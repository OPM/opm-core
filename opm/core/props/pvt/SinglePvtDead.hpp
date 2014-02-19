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

#ifndef OPM_SINGLEPVTDEAD_HEADER_INCLUDED
#define OPM_SINGLEPVTDEAD_HEADER_INCLUDED


#include <opm/core/props/pvt/SinglePvtInterface.hpp>
#include <opm/core/utility/NonuniformTableLinear.hpp>

#include <opm/parser/eclipse/Utility/PvdoTable.hpp>

#include <vector>

namespace Opm
{

    /// Class for immiscible dead oil and dry gas.
    /// The PVT properties can either be given as a function of pressure (p) and surface volume (z)
    /// or pressure (p) and gas resolution factor (r).
    /// For all the virtual methods, the following apply: p, r and z
    /// are expected to be of size n, size n and n*num_phases, respectively.
    /// Output arrays shall be of size n, and must be valid before
    /// calling the method.
    class SinglePvtDead : public SinglePvtInterface
    {
    public:
        typedef std::vector<std::vector<std::vector<double> > > table_t;
        SinglePvtDead(const table_t& pvd_table);
        SinglePvtDead(const Opm::PvdoTable &pvdoTable);
        virtual ~SinglePvtDead();

        /// Viscosity as a function of p and z.
        virtual void mu(const int n,
                        const double* p,
                        const double* z,
                        double* output_mu) const;

        /// Viscosity and its derivatives as a function of p and r.
        /// The fluid is considered saturated if r >= rsSat(p).
        virtual void mu(const int n,
                        const double* p,
                        const double* r,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const;

        /// Viscosity as a function of p and r.
        /// State condition determined by 'cond'.
        virtual void mu(const int n,
                        const double* p,
                        const double* r,
                        const PhasePresence* cond,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const;

        /// Formation volume factor as a function of p and z.
        virtual void B(const int n,
                       const double* p,
                       const double* z,
                       double* output_B) const;

        /// Formation volume factor and p-derivative as functions of p and z.
        virtual void dBdp(const int n,
                          const double* p,
                          const double* z,
                          double* output_B,
                          double* output_dBdp) const;

        /// The inverse of the formation volume factor b = 1 / B, and its derivatives as a function of p and r.
        /// The fluid is considered saturated if r >= rsSat(p).
        virtual void b(const int n,
                       const double* p,
                       const double* r,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const;

        /// The inverse of the formation volume factor b = 1 / B, and its derivatives as a function of p and r.
        /// State condition determined by 'cond'.
        virtual void b(const int n,
                       const double* p,
                       const double* r,
                       const PhasePresence* cond,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const;


        /// Solution gas/oil ratio and its derivatives at saturated conditions as a function of p.
        virtual void rsSat(const int n,
                          const double* p,
                          double* output_rsSat,
                          double* output_drsSatdp) const;

        /// Vapor oil/gas ratio and its derivatives at saturated conditions as a function of p.
        virtual void rvSat(const int n,
                          const double* p,
                          double* output_rvSat,
                          double* output_drvSatdp) const;

        /// Solution factor as a function of p and z.
        virtual void R(const int n,
                       const double* p,
                       const double* z,
                       double* output_R) const;

        /// Solution factor and p-derivative as functions of p and z.
        virtual void dRdp(const int n,
                          const double* p,
                          const double* z,
                          double* output_R,
                          double* output_dRdp) const;
    private:
        // PVT properties of dry gas or dead oil
        NonuniformTableLinear<double> b_;
        NonuniformTableLinear<double> viscosity_;
    };

}


#endif // OPM_SINGLEPVTDEAD_HEADER_INCLUDED
