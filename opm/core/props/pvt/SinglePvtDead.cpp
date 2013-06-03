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
#include <opm/core/props/pvt/SinglePvtDead.hpp>
#include <algorithm>

// Extra includes for debug dumping of tables.
// #include <boost/lexical_cast.hpp>
// #include <string>
// #include <fstream>

namespace Opm
{

    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------
    /// Constructor
    SinglePvtDead::SinglePvtDead(const table_t& pvd_table)
    {
        const int region_number = 0;
        if (pvd_table.size() != 1) {
            THROW("More than one PVT-region");
        }

        // Copy data
        const int sz = pvd_table[region_number][0].size();
        std::vector<double> press(sz);
        std::vector<double> b(sz);
        std::vector<double> visc(sz);
        for (int i = 0; i < sz; ++i) {
            press[i] = pvd_table[region_number][0][i];
            b[i] = 1.0 / pvd_table[region_number][1][i];
            visc[i]  = pvd_table[region_number][2][i];
        }
        b_ = NonuniformTableLinear<double>(press, b);
        viscosity_ = NonuniformTableLinear<double>(press, visc);

        // Dumping the created tables.
//         static int count = 0;
//         std::ofstream os((std::string("dump-") + boost::lexical_cast<std::string>(count++)).c_str());
//         os.precision(15);
//         os << "1/B\n\n" << one_over_B_
//            << "\n\nvisc\n\n" << viscosity_ << std::endl;
    }

    // Destructor
    SinglePvtDead::~SinglePvtDead()
    {
    }



    void SinglePvtDead::mu(const int n,
                           const double* p,
                           const double* /*z*/,
                           double* output_mu) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            output_mu[i] = viscosity_(p[i]);
        }
    }

    void SinglePvtDead::mu(const int n,
                               const double* p,
                               const double* /*r*/,
                               double* output_mu,
                               double* output_dmudp,
                               double* output_dmudr) const
        {
    // #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_mu[i] = viscosity_(p[i]);
                output_dmudp[i] = viscosity_.derivative(p[i]);
            }
            std::fill(output_dmudr, output_dmudr + n, 0.0);

        }

    void SinglePvtDead::B(const int n,
                          const double* p,
                          const double* /*z*/,
                          double* output_B) const
    {
// #pragma omp parallel for
        // B = 1/b
        for (int i = 0; i < n; ++i) {
            output_B[i] = 1.0/b_(p[i]);
        }
    }

    void SinglePvtDead::dBdp(const int n,
                             const double* p,
                             const double* /*z*/,
                             double* output_B,
                             double* output_dBdp) const
    {
        B(n, p, 0, output_B);
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double Bg = output_B[i];
            output_dBdp[i] = -Bg*Bg*b_.derivative(p[i]);
        }
    }

    void SinglePvtDead::b(const int n,
                              const double* p,
                              const double* /*r*/,
                              double* output_b,
                              double* output_dbdp,
                              double* output_dbdr) const

        {
    // #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                output_b[i] = b_(p[i]);
                output_dbdp[i] = b_.derivative(p[i]);

            }
            std::fill(output_dbdr, output_dbdr + n, 0.0);

        }

    void SinglePvtDead::rbub(const int n,
                             const double* /*p*/,
                             double* output_rbub,
                             double* output_drbubdp) const
    {
        std::fill(output_rbub, output_rbub + n, 0.0);
        std::fill(output_drbubdp, output_drbubdp + n, 0.0);
    }

    void SinglePvtDead::R(const int n,
                          const double* /*p*/,
                          const double* /*z*/,
                          double* output_R) const
    {
        std::fill(output_R, output_R + n, 0.0);
    }

    void SinglePvtDead::dRdp(const int n,
                             const double* /*p*/,
                             const double* /*z*/,
                             double* output_R,
                             double* output_dRdp) const
    {
        std::fill(output_R, output_R + n, 0.0);
        std::fill(output_dRdp, output_dRdp + n, 0.0);
    }

}
