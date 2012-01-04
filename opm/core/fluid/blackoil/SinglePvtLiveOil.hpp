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


#include <opm/core/fluid/blackoil/SinglePvtInterface.hpp>
#include <vector>

namespace Opm
{
    /// Class for miscible live oil.
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

	// PVT properties of live oil (with dissolved gas)
	std::vector<std::vector<double> > saturated_oil_table_;	
	std::vector<std::vector<std::vector<double> > > undersat_oil_tables_;
    };

}

#endif // OPM_SINGLEPVTLIVEOIL_HEADER_INCLUDED

