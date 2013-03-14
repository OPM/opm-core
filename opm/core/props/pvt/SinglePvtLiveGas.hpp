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

#ifndef OPM_SINGLEPVTLIVEGAS_HEADER_INCLUDED
#define OPM_SINGLEPVTLIVEGAS_HEADER_INCLUDED

#include <opm/core/fluid/blackoil/SinglePvtInterface.hpp>
#include <vector>

namespace Opm
{
    /// Class for miscible wet gas (with vaporized oil in vapour phase).
    /// For all the virtual methods, the following apply: p and z
    /// are expected to be of size n and n*num_phases, respectively.
    /// Output arrays shall be of size n, and must be valid before
    /// calling the method.
    class SinglePvtLiveGas : public SinglePvtInterface
    {
    public:
        typedef std::vector<std::vector<std::vector<double> > > table_t;

        SinglePvtLiveGas(const table_t& pvto);
        virtual ~SinglePvtLiveGas();

        /// Viscosity as a function of p and z.
        virtual void mu(const int n,
                        const double* p,
                        const double* z,
                        double* output_mu) const;

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

    protected:
        double evalB(double press, const double* surfvol) const;
        void evalBDeriv(double press, const double* surfvol, double& B, double& dBdp) const;
        double evalR(double press, const double* surfvol) const;
        void evalRDeriv(double press, const double* surfvol, double& R, double& dRdp) const;

        // item:  1=>B  2=>mu;
        double miscible_gas(const double press,
                            const double* surfvol,
                            const int item,
                            const bool deriv = false) const;
        // PVT properties of wet gas (with vaporised oil)
        std::vector<std::vector<double> > saturated_gas_table_;
        std::vector<std::vector<std::vector<double> > > undersat_gas_tables_;

    };

}

#endif // OPM_SINGLEPVTLIVEGAS_HEADER_INCLUDED

