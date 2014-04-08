//===========================================================================
//
// File: MiscibilityLiveGas.cpp
//
// Created: Wed Feb 10 09:21:53 2010
//
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//
// Revision: $Id$
//
//===========================================================================
/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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
#include <opm/core/props/pvt/SinglePvtLiveGas.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>

#include <algorithm>
#include <iostream>

namespace Opm
{

    using Opm::linearInterpolation;
    using Opm::linearInterpolationDerivative;


    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    SinglePvtLiveGas::SinglePvtLiveGas(const table_t& pvtg)
    {
        // GAS, PVTG
        const int region_number = 0;
        if (pvtg.size() != 1) {
            OPM_THROW(std::runtime_error, "More than one PVD-region");
        }
        saturated_gas_table_.resize(4);
        const int sz = pvtg[region_number].size();
        for (int k=0; k<4; ++k) {
            saturated_gas_table_[k].resize(sz);
        }

        for (int i=0; i<sz; ++i) {
            saturated_gas_table_[0][i] = pvtg[region_number][i][0];  // p
            saturated_gas_table_[1][i] = 1.0/pvtg[region_number][i][2];  // 1/Bg
            saturated_gas_table_[2][i] = pvtg[region_number][i][3];  // mu_g
            saturated_gas_table_[3][i] = pvtg[region_number][i][1]; // Rv
        }

        undersat_gas_tables_.resize(sz);
        for (int i=0; i<sz; ++i) {
            undersat_gas_tables_[i].resize(3);
            int tsize = (pvtg[region_number][i].size() - 1)/3;
            undersat_gas_tables_[i][0].resize(tsize);
            undersat_gas_tables_[i][1].resize(tsize);
            undersat_gas_tables_[i][2].resize(tsize);
            for (int j=0, k=0; j<tsize; ++j) {
                undersat_gas_tables_[i][0][j] = pvtg[region_number][i][++k]; // Rv
                undersat_gas_tables_[i][1][j] = 1.0/pvtg[region_number][i][++k]; // 1/Bg
                undersat_gas_tables_[i][2][j] = pvtg[region_number][i][++k]; // mu_g
            }
        }
    }

    SinglePvtLiveGas::SinglePvtLiveGas(const Opm::PvtgTable& pvtgTable)
    {
        // GAS, PVTG
        saturated_gas_table_.resize(4);
        saturated_gas_table_[0] = pvtgTable.getOuterTable()->getPressureColumn();
        saturated_gas_table_[1] = pvtgTable.getOuterTable()->getGasFormationFactorColumn();
        saturated_gas_table_[2] = pvtgTable.getOuterTable()->getGasViscosityColumn();
        saturated_gas_table_[3] = pvtgTable.getOuterTable()->getOilSolubilityColumn();

        int sz = pvtgTable.getOuterTable()->numRows();
        undersat_gas_tables_.resize(sz);
        for (int i=0; i<sz; ++i) {
            const auto &undersatTable = *pvtgTable.getInnerTable(i);

            undersat_gas_tables_[i].resize(3);
            undersat_gas_tables_[i][0] = undersatTable.getOilSolubilityColumn();
            undersat_gas_tables_[i][1] = undersatTable.getGasFormationFactorColumn();
            undersat_gas_tables_[i][2] = undersatTable.getGasViscosityColumn();
        }
        
        // Bg -> 1/Bg      
        for (int i=0; i<sz; ++i) {
            saturated_gas_table_[1][i] = 1.0/saturated_gas_table_[1][i];
            for (int j=0; j<undersat_gas_tables_[i][1].size(); ++j) {
                undersat_gas_tables_[i][1][j] = 1.0/undersat_gas_tables_[i][1][j];
            }
        }
    }

    // Destructor
     SinglePvtLiveGas::~SinglePvtLiveGas()
    {
    }


    void SinglePvtLiveGas::mu(const int n,
                              const double* p,
                              const double* z,
                              double* output_mu) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            output_mu[i] = miscible_gas(p[i], z + num_phases_*i, 2, false);
        }
    }

    /// Viscosity and its derivatives as a function of p and r.
    void SinglePvtLiveGas::mu(const int /*n*/,
                              const double* /*p*/,
                              const double* /*r*/,
                              double* /*output_mu*/,
                              double* /*output_dmudp*/,
                              double* /*output_dmudr*/) const
    {
        OPM_THROW(std::runtime_error, "The new fluid interface not yet implemented");
    }

    /// Viscosity and its derivatives as a function of p and r.
    void SinglePvtLiveGas::mu(const int n,
                               const double* p,
                               const double* r,
                               const PhasePresence* cond,
                               double* output_mu,
                               double* output_dmudp,
                               double* output_dmudr) const
    {
        for (int i = 0; i < n; ++i) {
            const PhasePresence& cnd = cond[i];
            output_mu[i] = miscible_gas(p[i], r[i], cnd,2, 0);
            output_dmudp[i] = miscible_gas(p[i], r[i], cnd, 2, 1);
            output_dmudr[i] = miscible_gas(p[i], r[i], cnd, 2, 2);
        }

    }


    /// Formation volume factor as a function of p and z.
    void SinglePvtLiveGas::B(const int n,
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
    void SinglePvtLiveGas::dBdp(const int n,
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

    /// The inverse of the formation volume factor b = 1 / B, and its derivatives as a function of p and r.
    void SinglePvtLiveGas::b(const int /*n*/,
                             const double* /*p*/,
                             const double* /*r*/,
                             double* /*output_b*/,
                             double* /*output_dbdp*/,
                             double* /*output_dbdr*/) const

    {
        OPM_THROW(std::runtime_error, "The new fluid interface not yet implemented");
    }

    /// The inverse of the formation volume factor b = 1 / B, and its derivatives as a function of p and r.
    void SinglePvtLiveGas::b(const int n,
                          const double* p,
                          const double* r,
                          const PhasePresence* cond,
                          double* output_b,
                          double* output_dbdp,
                          double* output_dbdr) const

    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    const PhasePresence& cnd = cond[i];

                    output_b[i] = miscible_gas(p[i], r[i], cnd, 1, 0);
                    output_dbdp[i] = miscible_gas(p[i], r[i], cnd, 1, 1);
                    output_dbdr[i] = miscible_gas(p[i], r[i], cnd, 1, 2);

                }
    }

    /// Gas resolution and its derivatives at bublepoint as a function of p.
    void SinglePvtLiveGas::rvSat(const int n,
                             const double* p,
                             double* output_rvSat,
                             double* output_drvSatdp) const
    {
        for (int i = 0; i < n; ++i) {
            output_rvSat[i] = linearInterpolation(saturated_gas_table_[0],
                    saturated_gas_table_[3],p[i]);
            output_drvSatdp[i] = linearInterpolationDerivative(saturated_gas_table_[0],
                    saturated_gas_table_[3],p[i]);

        }
    }

    void SinglePvtLiveGas::rsSat(const int n,
                             const double* /*p*/,
                             double* output_rsSat,
                             double* output_drsSatdp) const
    {
        std::fill(output_rsSat, output_rsSat + n, 0.0);
        std::fill(output_drsSatdp, output_drsSatdp + n, 0.0);
    }

    /// Solution factor as a function of p and z.
    void SinglePvtLiveGas::R(const int n,
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
    void SinglePvtLiveGas::dRdp(const int n,
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

    double SinglePvtLiveGas::evalB(const double press, const double* surfvol) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            // To handle no-gas case.
            return 1.0;
        }
        return 1.0/miscible_gas(press, surfvol, 1, false);
    }

    void SinglePvtLiveGas::evalBDeriv(const double press, const double* surfvol,
                                      double& Bval, double& dBdpval) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            // To handle no-gas case.
            Bval = 1.0;
            dBdpval = 0.0;
            return;
        }
        Bval = evalB(press, surfvol);
        dBdpval =  -Bval*Bval*miscible_gas(press, surfvol, 1, true);
    }

    double SinglePvtLiveGas::evalR(const double press, const double* surfvol) const
    {
        if (surfvol[phase_pos_[Liquid]] == 0.0) {
            // To handle no-gas case.
            return 0.0;
        }
        double satR = linearInterpolation(saturated_gas_table_[0],
                                             saturated_gas_table_[3], press);
        double maxR = surfvol[phase_pos_[Liquid]]/surfvol[phase_pos_[Vapour]];
        if (satR < maxR ) {
            // Saturated case
            return satR;
        } else {
            // Undersaturated case
            return maxR;
        }
    }

    void SinglePvtLiveGas::evalRDeriv(const double press, const double* surfvol,
                                      double& Rval, double& dRdpval) const
    {
        if (surfvol[phase_pos_[Liquid]] == 0.0) {
            // To handle no-gas case.
            Rval = 0.0;
            dRdpval = 0.0;
            return;
        }
        double satR = linearInterpolation(saturated_gas_table_[0],
                                             saturated_gas_table_[3], press);
        double maxR = surfvol[phase_pos_[Liquid]]/surfvol[phase_pos_[Vapour]];
        if (satR < maxR ) {
            // Saturated case
            Rval = satR;
            dRdpval = linearInterpolationDerivative(saturated_gas_table_[0],
                                            saturated_gas_table_[3],
                                            press);
        } else {
            // Undersaturated case
            Rval = maxR;
            dRdpval = 0.0;
        }
    }

    double SinglePvtLiveGas::miscible_gas(const double press,
                                          const double* surfvol,
                                          const int item,
                                          const bool deriv) const
    {
        int section;
        double Rval = linearInterpolation(saturated_gas_table_[0],
                                                saturated_gas_table_[3], press,
                                                section);
        double maxR = surfvol[phase_pos_[Liquid]]/surfvol[phase_pos_[Vapour]];
        if (deriv) {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolationDerivative(saturated_gas_table_[0],
                                                saturated_gas_table_[item],
                                                press);
            } else {  // Undersaturated case
                int is = section;
                if (undersat_gas_tables_[is][0].size() < 2) {
                    double val = (saturated_gas_table_[item][is+1]
                                  - saturated_gas_table_[item][is]) /
                        (saturated_gas_table_[0][is+1] -
                         saturated_gas_table_[0][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersat_gas_tables_[is][0],
                                              undersat_gas_tables_[is][item],
                                              maxR);
                double val2 =
                    linearInterpolation(undersat_gas_tables_[is+1][0],
                                              undersat_gas_tables_[is+1][item],
                                              maxR);
                double val = (val2 - val1)/
                    (saturated_gas_table_[0][is+1] - saturated_gas_table_[0][is]);
                return val;
            }
        } else {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolation(saturated_gas_table_[0],
                                                 saturated_gas_table_[item],
                                                 press);
            } else {  // Undersaturated case
                int is = section;
                // Extrapolate from first table section
                if (is == 0 && press < saturated_gas_table_[0][0]) {
                    return linearInterpolation(undersat_gas_tables_[0][0],
                                                     undersat_gas_tables_[0][item],
                                                     maxR);
                }

                // Extrapolate from last table section
                int ltp = saturated_gas_table_[0].size() - 1;
                if (is+1 == ltp && press > saturated_gas_table_[0][ltp]) {
                    return linearInterpolation(undersat_gas_tables_[ltp][0],
                                                     undersat_gas_tables_[ltp][item],
                                                     maxR);
                }

                // Interpolate between table sections
                double w = (press - saturated_gas_table_[0][is]) /
                    (saturated_gas_table_[0][is+1] -
                     saturated_gas_table_[0][is]);
                if (undersat_gas_tables_[is][0].size() < 2) {
                    double val = saturated_gas_table_[item][is] +
                        w*(saturated_gas_table_[item][is+1] -
                           saturated_gas_table_[item][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersat_gas_tables_[is][0],
                                              undersat_gas_tables_[is][item],
                                              maxR);
                double val2 =
                    linearInterpolation(undersat_gas_tables_[is+1][0],
                                              undersat_gas_tables_[is+1][item],
                                              maxR);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }

    double SinglePvtLiveGas::miscible_gas(const double press,
                                          const double r,
                                          const PhasePresence& cond,
                                          const int item,
                                          const int deriv) const
{
    const bool isSat = cond.hasFreeOil();

        // Derivative w.r.t p
        if (deriv == 1) {
            if (isSat) {  // Saturated case
                return linearInterpolationDerivative(saturated_gas_table_[0],
                                                saturated_gas_table_[item],
                                                press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_gas_table_[0], press);
                if (undersat_gas_tables_[is][0].size() < 2) {
                    double val = (saturated_gas_table_[item][is+1]
                                  - saturated_gas_table_[item][is]) /
                        (saturated_gas_table_[0][is+1] -
                         saturated_gas_table_[0][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersat_gas_tables_[is][0],
                                              undersat_gas_tables_[is][item],
                                              r);
                double val2 =
                    linearInterpolation(undersat_gas_tables_[is+1][0],
                                              undersat_gas_tables_[is+1][item],
                                              r);
                double val = (val2 - val1)/
                    (saturated_gas_table_[0][is+1] - saturated_gas_table_[0][is]);
                return val;
            }
        } else if (deriv == 2){
            if (isSat) {
                return 0;
            } else {
                int is = tableIndex(saturated_gas_table_[0], press);
                double w = (press - saturated_gas_table_[0][is]) /
                    (saturated_gas_table_[0][is+1] - saturated_gas_table_[0][is]);
                assert(undersat_gas_tables_[is][0].size() >= 2);
                assert(undersat_gas_tables_[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersat_gas_tables_[is][0],
                                              undersat_gas_tables_[is][item],
                                              r);
                double val2 =
                    linearInterpolationDerivative(undersat_gas_tables_[is+1][0],
                                              undersat_gas_tables_[is+1][item],
                                              r);

                double val = val1 + w * (val2 - val1);
                return val;

            }
        } else {
            if (isSat) {  // Saturated case
                return linearInterpolation(saturated_gas_table_[0],
                                                 saturated_gas_table_[item],
                                                 press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_gas_table_[0], press);
                // Extrapolate from first table section
                if (is == 0 && press < saturated_gas_table_[0][0]) {
                    return linearInterpolation(undersat_gas_tables_[0][0],
                                                     undersat_gas_tables_[0][item],
                                                     r);
                }

                // Extrapolate from last table section
                //int ltp = saturated_gas_table_[0].size() - 1;
                //if (is+1 == ltp && press > saturated_gas_table_[0][ltp]) {
                //    return linearInterpolation(undersat_gas_tables_[ltp][0],
                //                                    undersat_gas_tables_[ltp][item],
                //                                    r);
                //}

                // Interpolate between table sections
                double w = (press - saturated_gas_table_[0][is]) /
                    (saturated_gas_table_[0][is+1] -
                     saturated_gas_table_[0][is]);
                if (undersat_gas_tables_[is][0].size() < 2) {
                    double val = saturated_gas_table_[item][is] +
                        w*(saturated_gas_table_[item][is+1] -
                           saturated_gas_table_[item][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersat_gas_tables_[is][0],
                                              undersat_gas_tables_[is][item],
                                              r);
                double val2 =
                    linearInterpolation(undersat_gas_tables_[is+1][0],
                                              undersat_gas_tables_[is+1][item],
                                              r);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }




} // namespace Opm
