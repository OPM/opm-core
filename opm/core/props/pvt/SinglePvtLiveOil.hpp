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

#ifndef OPM_SINGLEPVTLIVEOIL_HEADER_INCLUDED
#define OPM_SINGLEPVTLIVEOIL_HEADER_INCLUDED


#include <opm/core/props/pvt/SinglePvtInterface.hpp>
#include <vector>

namespace Opm
{
    /// Class for miscible live oil (with dissolved gas in liquid phase).
    /// The PVT properties can either be given as a function of pressure (p) and surface volume (z)
    /// or pressure (p) and gas resolution factor (r).
    /// For all the virtual methods, the following apply: p, r and z
    /// are expected to be of size n, size n and n*num_phases, respectively.
    /// Output arrays shall be of size n, and must be valid before
    /// calling the method.
    class SinglePvtLiveOil : public SinglePvtInterface
    {
    public:
        typedef std::vector<std::vector<std::vector<double> > > table_t;

        SinglePvtLiveOil(const table_t& pvto);
        virtual ~SinglePvtLiveOil();

        /// Viscosity as a function of p and z.
        virtual void mu(const int n,
                        const double* p,
                        const double* z,
                        double* output_mu) const;

        /// Viscosity and its derivatives as a function of p and r.
        /// The fluid is considered saturated if r >= rbub(p).
        virtual void mu(const int n,
                        const double* p,
                        const double* r,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const;

        /// Viscosity as a function of p and r.
        /// Whether the fluid is saturated or not is given explicitly by isSat.
        virtual void mu(const int n,
                        const double* p,
                        const double* r,
                        const bool* isSat,
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
        /// The fluid is considered saturated if r >= rbub(p).
        virtual void b(const int n,
                       const double* p,
                       const double* r,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const;

        /// The inverse of the formation volume factor b = 1 / B, and its derivatives as a function of p and r.
        /// Whether the fluid is saturated or not is given explicitly by isSat.
        virtual void b(const int n,
                       const double* p,
                       const double* r,
                       const bool* isSat,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const;

        /// Gas resolution and its derivatives at bublepoint as a function of p.
        virtual void rbub(const int n,
                          const double* p,
                          double* output_rbub,
                          double* output_drbubdp) const;

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
        double evalB(double press, const double* surfvol) const;
        void evalBDeriv(double press, const double* surfvol, double& B, double& dBdp) const;
        double evalR(double press, const double* surfvol) const;
        void evalRDeriv(double press, const double* surfvol, double& R, double& dRdp) const;

        // item:  1=>1/B  2=>mu;
        double miscible_oil(const double press,
                            const double* surfvol,
                            const int item,
                            const bool deriv = false) const;

        double miscible_oil(const double press,
                            const double r,
                            const int item,
                            const int deriv = 0) const;

        double miscible_oil(const double press,
                            const double r,
                            const bool isSat,
                            const int item,
                            const int deriv = 0) const;

        // PVT properties of live oil (with dissolved gas)
        std::vector<std::vector<double> > saturated_oil_table_;
        std::vector<std::vector<std::vector<double> > > undersat_oil_tables_;
    };

}

#endif // OPM_SINGLEPVTLIVEOIL_HEADER_INCLUDED

