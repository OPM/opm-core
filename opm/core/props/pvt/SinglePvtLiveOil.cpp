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

#include "config.h"
#include <opm/core/props/pvt/SinglePvtLiveOil.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <algorithm>


namespace Opm
{

    using Opm::linearInterpolation;
    using Opm::linearInterpolationDerivative;
    using Opm::tableIndex;


    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    SinglePvtLiveOil::SinglePvtLiveOil(const table_t& pvto)
    {
        // OIL, PVTO
        const int region_number = 0;
        if (pvto.size() != 1) {
            OPM_THROW(std::runtime_error, "More than one PVD-region");
        }
        saturated_oil_table_.resize(4);
        const int sz =  pvto[region_number].size();
        for (int k=0; k<4; ++k) {
            saturated_oil_table_[k].resize(sz);
        }
        for (int i=0; i<sz; ++i) {
            saturated_oil_table_[0][i] = pvto[region_number][i][1]; // p
            saturated_oil_table_[1][i] = 1.0/pvto[region_number][i][2]; // 1/Bo
            saturated_oil_table_[2][i] = pvto[region_number][i][3];   // mu_o
            saturated_oil_table_[3][i] = pvto[region_number][i][0];     // Rs
        }

        undersat_oil_tables_.resize(sz);
        for (int i=0; i<sz; ++i) {
            undersat_oil_tables_[i].resize(3);
            int tsize = (pvto[region_number][i].size() - 1)/3;
            undersat_oil_tables_[i][0].resize(tsize);
            undersat_oil_tables_[i][1].resize(tsize);
            undersat_oil_tables_[i][2].resize(tsize);
            for (int j=0, k=0; j<tsize; ++j) {
                undersat_oil_tables_[i][0][j] = pvto[region_number][i][++k];  // p
                undersat_oil_tables_[i][1][j] = 1.0/pvto[region_number][i][++k];  // 1/Bo
                undersat_oil_tables_[i][2][j] = pvto[region_number][i][++k];  // mu_o
            }
        }

        // Complete undersaturated tables by extrapolating from existing data
        // as is done in Eclipse and Mrst
        int iNext = -1;
        for (int i=0; i<sz; ++i) {
            // Skip records already containing undersaturated data
            if (undersat_oil_tables_[i][0].size() > 1) {
                continue;
            }
            // Look ahead for next record containing undersaturated data
            if (iNext < i) {
                iNext = i+1;
                while (iNext<sz && undersat_oil_tables_[iNext][0].size() < 2) {
                    ++iNext;
                }
                if (iNext == sz) OPM_THROW(std::runtime_error,"Unable to complete undersaturated table.");
            }
            // Add undersaturated data to current record while maintaining compressibility and viscosibility
            typedef std::vector<std::vector<std::vector<double> > >::size_type sz_t;
            for (sz_t j=1; j<undersat_oil_tables_[iNext][0].size(); ++j) {
                double diffPressure = undersat_oil_tables_[iNext][0][j]-undersat_oil_tables_[iNext][0][j-1];
                double pressure = undersat_oil_tables_[i][0].back()+diffPressure;
                undersat_oil_tables_[i][0].push_back(pressure);
                double compr = (1.0/undersat_oil_tables_[iNext][1][j]-1.0/undersat_oil_tables_[iNext][1][j-1])
                        / (0.5*(1.0/undersat_oil_tables_[iNext][1][j]+1.0/undersat_oil_tables_[iNext][1][j-1]));
                double B = (1.0/undersat_oil_tables_[i][1].back())*(1.0+0.5*compr)/(1.0-0.5*compr);
                undersat_oil_tables_[i][1].push_back(1.0/B);
                double visc = (undersat_oil_tables_[iNext][2][j]-undersat_oil_tables_[iNext][2][j-1])
                        / (0.5*(undersat_oil_tables_[iNext][2][j]+undersat_oil_tables_[iNext][2][j-1]));
                double mu = (undersat_oil_tables_[i][2].back())*(1.0+0.5*visc)/(1.0-0.5*visc);
                undersat_oil_tables_[i][2].push_back(mu);
            }
        }

    }

    /// Destructor.
    SinglePvtLiveOil::~SinglePvtLiveOil()
    {
    }


    /// Viscosity as a function of p and z.
    void SinglePvtLiveOil::mu(const int n,
                              const double* p,
                              const double* z,
                              double* output_mu) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            output_mu[i] = miscible_oil(p[i], z + num_phases_*i, 2, false);
        }
    }

    /// Viscosity and its derivatives as a function of p and r.
    void SinglePvtLiveOil::mu(const int n,
                               const double* p,
                               const double* r,
                               double* output_mu,
                               double* output_dmudp,
                               double* output_dmudr) const
    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    output_mu[i] = miscible_oil(p[i], r[i], 2, 0);
                    output_dmudp[i] = miscible_oil(p[i], r[i], 2, 1);
                    output_dmudr[i] = miscible_oil(p[i], r[i], 2, 2);

                }
    }

    /// Viscosity and its derivatives as a function of p and r.
    void SinglePvtLiveOil::mu(const int n,
                               const double* p,
                               const double* r,
                               const bool* isSat,
                               double* output_mu,
                               double* output_dmudp,
                               double* output_dmudr) const
    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    output_mu[i] = miscible_oil(p[i], r[i],isSat[i], 2, 0);
                    output_dmudp[i] = miscible_oil(p[i], r[i],isSat[i], 2, 1);
                    output_dmudr[i] = miscible_oil(p[i], r[i],isSat[i], 2, 2);

                }
    }


    /// Formation volume factor as a function of p and z.
    void SinglePvtLiveOil::B(const int n,
                             const double* p,
                             const double* z,
                             double* output_B) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            output_B[i] = evalB(p[i], z + num_phases_*i);
        }

    }


    /// Formation volume factor and p-derivative as functions of p and z.
    void SinglePvtLiveOil::dBdp(const int n,
                                const double* p,
                                const double* z,
                                double* output_B,
                                double* output_dBdp) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            evalBDeriv(p[i], z + num_phases_*i, output_B[i], output_dBdp[i]);
        }
    }

    void SinglePvtLiveOil::b(const int n,
                          const double* p,
                          const double* r,
                          double* output_b,
                          double* output_dbdp,
                          double* output_dbdr) const

    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    output_b[i] = miscible_oil(p[i], r[i], 1, 0);
                    output_dbdp[i] = miscible_oil(p[i], r[i], 1, 1);
                    output_dbdr[i] = miscible_oil(p[i], r[i], 1, 2);

                }
    }

    void SinglePvtLiveOil::b(const int n,
                          const double* p,
                          const double* r,
                          const bool* isSat,
                          double* output_b,
                          double* output_dbdp,
                          double* output_dbdr) const

    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    output_b[i] = miscible_oil(p[i], r[i], isSat[i],1, 0);
                    output_dbdp[i] = miscible_oil(p[i], r[i], isSat[i], 1, 1);
                    output_dbdr[i] = miscible_oil(p[i], r[i], isSat[i],1, 2);

                }
    }

    void SinglePvtLiveOil::rbub(const int n,
                             const double* p,
                             double* output_rbub,
                             double* output_drbubdp) const
    {

        for (int i = 0; i < n; ++i) {
            output_rbub[i] = linearInterpolation(saturated_oil_table_[0],
                    saturated_oil_table_[3],p[i]);
            output_drbubdp[i] = linearInterpolationDerivative(saturated_oil_table_[0],
                    saturated_oil_table_[3],p[i]);

        }
    }

    /// Solution factor as a function of p and z.
    void SinglePvtLiveOil::R(const int n,
                             const double* p,
                             const double* z,
                             double* output_R) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            output_R[i] = evalR(p[i], z + num_phases_*i);
        }

    }


    /// Solution factor and p-derivative as functions of p and z.
    void SinglePvtLiveOil::dRdp(const int n,
                                const double* p,
                                const double* z,
                                double* output_R,
                                double* output_dRdp) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            evalRDeriv(p[i], z + num_phases_*i, output_R[i], output_dRdp[i]);
        }
    }




    // ---- Private methods ----

    double SinglePvtLiveOil::evalB(double press, const double* surfvol) const
    {
        // if (surfvol[phase_pos_[Liquid]] == 0.0) return 1.0; // To handle no-oil case.
        return 1.0/miscible_oil(press, surfvol, 1, false);
    }


    void SinglePvtLiveOil::evalBDeriv(const double press, const double* surfvol,
                                      double& Bval, double& dBdpval) const
    {
        Bval = evalB(press, surfvol);
        dBdpval = -Bval*Bval*miscible_oil(press, surfvol, 1, true);
    }

    double SinglePvtLiveOil::evalR(double press, const double* surfvol) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            return 0.0;
        }
        double Rval = linearInterpolation(saturated_oil_table_[0],
                                             saturated_oil_table_[3], press);
        double maxR = surfvol[phase_pos_[Vapour]]/surfvol[phase_pos_[Liquid]];
        if (Rval < maxR ) {  // Saturated case
            return Rval;
        } else {
            return maxR;  // Undersaturated case
        }
    }

    void SinglePvtLiveOil::evalRDeriv(const double press, const double* surfvol,
                                      double& Rval, double& dRdpval) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            Rval = 0.0;
            dRdpval = 0.0;
            return;
        }
        Rval = linearInterpolation(saturated_oil_table_[0],
                                      saturated_oil_table_[3], press);
        double maxR = surfvol[phase_pos_[Vapour]]/surfvol[phase_pos_[Liquid]];
        if (Rval < maxR ) {
            // Saturated case
            dRdpval = linearInterpolationDerivative(saturated_oil_table_[0],
                                            saturated_oil_table_[3],
                                            press);
        } else {
            // Undersaturated case
            Rval = maxR;
            dRdpval = 0.0;
        }
    }


    double SinglePvtLiveOil::miscible_oil(const double press,
                                          const double* surfvol,
                                          const int item,
                                          const bool deriv) const
    {
        int section;
        double Rval = linearInterpolation(saturated_oil_table_[0],
                                                saturated_oil_table_[3],
                                                press, section);
        double maxR = (surfvol[phase_pos_[Liquid]] == 0.0) ? 0.0 : surfvol[phase_pos_[Vapour]]/surfvol[phase_pos_[Liquid]];
        if (deriv) {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolationDerivative(saturated_oil_table_[0],
                                                saturated_oil_table_[item],
                                                press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[3], maxR);
                double w = (maxR - saturated_oil_table_[3][is]) /
                    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersat_oil_tables_[is][0],
                                             undersat_oil_tables_[is][item],
                                             press);
                double val2 =
                    linearInterpolationDerivative(undersat_oil_tables_[is+1][0],
                                             undersat_oil_tables_[is+1][item],
                                             press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        } else {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolation(saturated_oil_table_[0],
                                                 saturated_oil_table_[item],
                                                 press);
            } else {  // Undersaturated case
                // Interpolate between table sections
                int is = tableIndex(saturated_oil_table_[3], maxR);
                double w = (maxR - saturated_oil_table_[3][is]) /
                    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[is][0],
                                              undersat_oil_tables_[is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[is+1][0],
                                              undersat_oil_tables_[is+1][item],
                                              press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }

    double SinglePvtLiveOil::miscible_oil(const double press,
                                          const double r,
                                          const int item,
                                          const int deriv) const
    {
        int section;
        double Rval = linearInterpolation(saturated_oil_table_[0],
                                                saturated_oil_table_[3],
                                                press, section);
        // derivative with respect to frist component (pressure)
        if (deriv == 1) {
            if (Rval <= r ) {  // Saturated case
                return linearInterpolationDerivative(saturated_oil_table_[0],
                                                saturated_oil_table_[item],
                                                press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[3], r);
                double w = (r - saturated_oil_table_[3][is]) /
                    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersat_oil_tables_[is][0],
                                             undersat_oil_tables_[is][item],
                                             press);
                double val2 =
                    linearInterpolationDerivative(undersat_oil_tables_[is+1][0],
                                             undersat_oil_tables_[is+1][item],
                                             press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
            // derivative with respect to second component (r)
        } else if (deriv == 2)  {
            if (Rval <= r ) {  // Saturated case
                return 0;
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[3], r);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[is][0],
                                              undersat_oil_tables_[is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[is+1][0],
                                              undersat_oil_tables_[is+1][item],
                                              press);

                double val = (val2 - val1)/(saturated_oil_table_[3][is+1]-saturated_oil_table_[3][is]);
                return val;
            }


            }    else {
            if (Rval <= r ) {  // Saturated case
                return linearInterpolation(saturated_oil_table_[0],
                                                 saturated_oil_table_[item],
                                                 press);
            } else {  // Undersaturated case
                // Interpolate between table sections
                int is = tableIndex(saturated_oil_table_[3], r);
                double w = (r - saturated_oil_table_[3][is]) /
                    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[is][0],
                                              undersat_oil_tables_[is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[is+1][0],
                                              undersat_oil_tables_[is+1][item],
                                              press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }

    double SinglePvtLiveOil::miscible_oil(const double press,
                                          const double r,
                                          const bool isSat,
                                          const int item,
                                          const int deriv) const
    {
        // derivative with respect to frist component (pressure)
        if (deriv == 1) {
            if (isSat) {  // Saturated case
                return linearInterpolationDerivative(saturated_oil_table_[0],
                                                saturated_oil_table_[item],
                                                press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[3], r);
                double w = (r - saturated_oil_table_[3][is]) /
                    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersat_oil_tables_[is][0],
                                             undersat_oil_tables_[is][item],
                                             press);
                double val2 =
                    linearInterpolationDerivative(undersat_oil_tables_[is+1][0],
                                             undersat_oil_tables_[is+1][item],
                                             press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
            // derivative with respect to second component (r)
        } else if (deriv == 2)  {
            if (isSat) {  // Saturated case
                return 0;
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[3], r);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[is][0],
                                              undersat_oil_tables_[is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[is+1][0],
                                              undersat_oil_tables_[is+1][item],
                                              press);

                double val = (val2 - val1)/(saturated_oil_table_[3][is+1]-saturated_oil_table_[3][is]);
                return val;
            }


            }    else {
            if (isSat) {  // Saturated case
                return linearInterpolation(saturated_oil_table_[0],
                                                 saturated_oil_table_[item],
                                                 press);
            } else {  // Undersaturated case
                // Interpolate between table sections
                int is = tableIndex(saturated_oil_table_[3], r);
                double w = (r - saturated_oil_table_[3][is]) /
                    (saturated_oil_table_[3][is+1] - saturated_oil_table_[3][is]);
                assert(undersat_oil_tables_[is][0].size() >= 2);
                assert(undersat_oil_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[is][0],
                                              undersat_oil_tables_[is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[is+1][0],
                                              undersat_oil_tables_[is+1][item],
                                              press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }

} // namespace Opm
