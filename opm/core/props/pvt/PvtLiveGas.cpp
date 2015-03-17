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
#include <opm/core/props/pvt/PvtLiveGas.hpp>
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
    PvtLiveGas::PvtLiveGas(const std::vector<Opm::PvtgTable>& pvtgTables)
    {
        int numTables = pvtgTables.size();
        saturated_gas_table_.resize(numTables);
        undersat_gas_tables_.resize(numTables);

        for (int pvtTableIdx = 0; pvtTableIdx < numTables; ++pvtTableIdx) {
            const Opm::PvtgTable& pvtgTable = pvtgTables[pvtTableIdx];

            // GAS, PVTG
            // saturated_gas_table_[pvtTableIdx].resize(4);
            // Adding one extra line to PVTG to store 1./(Bg*mu_g)
            saturated_gas_table_[pvtTableIdx].resize(5);
            saturated_gas_table_[pvtTableIdx][0] = pvtgTable.getOuterTable()->getPressureColumn(); // Pressure
            saturated_gas_table_[pvtTableIdx][1] = pvtgTable.getOuterTable()->getGasFormationFactorColumn(); // Bg
            saturated_gas_table_[pvtTableIdx][2] = pvtgTable.getOuterTable()->getGasViscosityColumn(); // mu_g
            // The number of the columns
            int nRows = saturated_gas_table_[pvtTableIdx][2].size();
            saturated_gas_table_[pvtTableIdx][3].resize(nRows); // allocate memory for 1/(Bg*mu_g)
            saturated_gas_table_[pvtTableIdx][4] = pvtgTable.getOuterTable()->getOilSolubilityColumn(); // Rv

            int sz = pvtgTable.getOuterTable()->numRows();
            undersat_gas_tables_[pvtTableIdx].resize(sz);
            for (int i=0; i<sz; ++i) {
                const auto &undersatTable = *pvtgTable.getInnerTable(i);

                // undersat_gas_tables_[pvtTableIdx][i].resize(3);
                undersat_gas_tables_[pvtTableIdx][i].resize(4);
                undersat_gas_tables_[pvtTableIdx][i][0] = undersatTable.getOilSolubilityColumn(); // Rv
                undersat_gas_tables_[pvtTableIdx][i][1] = undersatTable.getGasFormationFactorColumn(); // Bg
                undersat_gas_tables_[pvtTableIdx][i][2] = undersatTable.getGasViscosityColumn();  // mu_g
                int nUndersatRows = undersat_gas_tables_[pvtTableIdx][i][2].size();
                undersat_gas_tables_[pvtTableIdx][i][3].resize(nUndersatRows); // allocate memory for 1/(Bg*mu_g)
            }

            // Bg -> 1/Bg
            for (int i=0; i<sz; ++i) {
                saturated_gas_table_[pvtTableIdx][3][i] = 1.0 / (saturated_gas_table_[pvtTableIdx][1][i]
                                                               * saturated_gas_table_[pvtTableIdx][2][i]); // 1/(Bg*mu_g)
                saturated_gas_table_[pvtTableIdx][1][i] = 1.0 / saturated_gas_table_[pvtTableIdx][1][i];
                for (size_t j=0; j<undersat_gas_tables_[pvtTableIdx][i][1].size(); ++j) {
                    undersat_gas_tables_[pvtTableIdx][i][3][j] = 1.0 / (undersat_gas_tables_[pvtTableIdx][i][1][j]
                                                                      * undersat_gas_tables_[pvtTableIdx][i][2][j]); // 1/(Bg*mu_g)
                    undersat_gas_tables_[pvtTableIdx][i][1][j] = 1.0 / undersat_gas_tables_[pvtTableIdx][i][1][j];
                }
            }
        }
    }

    // Destructor
     PvtLiveGas::~PvtLiveGas()
    {
    }


    void PvtLiveGas::mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* /*T*/,
                        const double* z,
                        double* output_mu) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double inverseB = miscible_gas(p[i], z + num_phases_*i, getTableIndex_(pvtRegionIdx, i), 1, false);
            double inverseBMu = miscible_gas(p[i], z + num_phases_*i, getTableIndex_(pvtRegionIdx, i), 3, false);

            output_mu[i] = inverseB / inverseBMu;
        }
    }

    /// Viscosity and its p and r derivatives as a function of p, T and r.
    void PvtLiveGas::mu(const int /*n*/,
                        const int* /*pvtRegionIdx*/,
                        const double* /*p*/,
                        const double* /*T*/,
                              const double* /*r*/,
                              double* /*output_mu*/,
                              double* /*output_dmudp*/,
                              double* /*output_dmudr*/) const
    {
        OPM_THROW(std::runtime_error, "The new fluid interface not yet implemented");
    }

    /// Viscosity and its p and r derivatives as a function of p, T and r.
    void PvtLiveGas::mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* /*T*/,
                               const double* r,
                               const PhasePresence* cond,
                               double* output_mu,
                               double* output_dmudp,
                               double* output_dmudr) const
    {
        for (int i = 0; i < n; ++i) {
            const PhasePresence& cnd = cond[i];
            int tableIdx = getTableIndex_(pvtRegionIdx, i);

            double inverseBMu = miscible_gas(p[i], r[i], cnd, tableIdx, 3, 0);
            double inverseB = miscible_gas(p[i], r[i], cnd, tableIdx, 1, 0);

            output_mu[i] = inverseB / inverseBMu;

            double dinverseBmudp = miscible_gas(p[i], r[i], cnd, tableIdx, 3, 1);
            double dinverseBdp = miscible_gas(p[i], r[i], cnd, tableIdx, 1, 1);

            output_dmudp[i] = (inverseBMu * dinverseBdp - inverseB * dinverseBmudp)
                              / (inverseBMu * inverseBMu);

            double dinverseBmudr = miscible_gas(p[i], r[i], cnd, tableIdx, 3, 2);
            double dinverseBdr = miscible_gas(p[i], r[i], cnd, tableIdx, 1, 2);

            output_dmudr[i] = (inverseBMu * dinverseBdr - inverseB * dinverseBmudr)
                              / (inverseBMu * inverseBMu);

        }

    }


    /// Formation volume factor as a function of p, T and z.
    void PvtLiveGas::B(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* /*T*/,
                             const double* z,
                             double* output_B) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            output_B[i] = evalB(p[i], z + num_phases_*i, getTableIndex_(pvtRegionIdx, i));
        }

    }


    /// Formation volume factor and p-derivative as functions of p and z.
    void PvtLiveGas::dBdp(const int n,
                          const int* pvtRegionIdx,
                          const double* p,
                          const double* /*T*/,
                                const double* z,
                                double* output_B,
                                double* output_dBdp) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            evalBDeriv(p[i], z + num_phases_*i, getTableIndex_(pvtRegionIdx, i), output_B[i], output_dBdp[i]);
        }
    }

    /// The inverse of the formation volume factor b = 1 / B, and its p and r derivatives as a function of p, T and r.
    void PvtLiveGas::b(const int /*n*/,
                       const int* /*pvtRegionIdx*/,
                       const double* /*p*/,
                       const double* /*T*/,
                             const double* /*r*/,
                             double* /*output_b*/,
                             double* /*output_dbdp*/,
                             double* /*output_dbdr*/) const

    {
        OPM_THROW(std::runtime_error, "The new fluid interface not yet implemented");
    }

    /// The inverse of the formation volume factor b = 1 / B, and its p and r derivatives as a function of p, T and r.
    void PvtLiveGas::b(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* /*T*/,
                          const double* r,
                          const PhasePresence* cond,
                          double* output_b,
                          double* output_dbdp,
                          double* output_dbdr) const

    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    const PhasePresence& cnd = cond[i];

                    int tableIdx = getTableIndex_(pvtRegionIdx, i);
                    output_b[i] = miscible_gas(p[i], r[i], cnd, tableIdx, 1, 0);
                    output_dbdp[i] = miscible_gas(p[i], r[i], cnd, tableIdx, 1, 1);
                    output_dbdr[i] = miscible_gas(p[i], r[i], cnd, tableIdx, 1, 2);

                }
    }

    /// Gas resolution and its derivatives at bublepoint as a function of p.
    void PvtLiveGas::rvSat(const int n,
                           const int* pvtRegionIdx,
                             const double* p,
                             double* output_rvSat,
                             double* output_drvSatdp) const
    {
        for (int i = 0; i < n; ++i) {
            int pvtTableIdx = getTableIndex_(pvtRegionIdx, i);
            output_rvSat[i] = linearInterpolation(saturated_gas_table_[pvtTableIdx][0],
                    saturated_gas_table_[pvtTableIdx][4],p[i]);
            output_drvSatdp[i] = linearInterpolationDerivative(saturated_gas_table_[pvtTableIdx][0],
                    saturated_gas_table_[pvtTableIdx][4],p[i]);

        }
    }

    void PvtLiveGas::rsSat(const int n,
                           const int* /*pvtRegionIdx*/,
                             const double* /*p*/,
                             double* output_rsSat,
                             double* output_drsSatdp) const
    {
        std::fill(output_rsSat, output_rsSat + n, 0.0);
        std::fill(output_drsSatdp, output_drsSatdp + n, 0.0);
    }

    /// Solution factor as a function of p and z.
    void PvtLiveGas::R(const int n,
                       const int* pvtRegionIdx,
                             const double* p,
                             const double* z,
                             double* output_R) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            output_R[i] = evalR(p[i], z + num_phases_*i, getTableIndex_(pvtRegionIdx, i));
        }

    }


    /// Solution factor and p-derivative as functions of p and z.
    void PvtLiveGas::dRdp(const int n,
                          const int* pvtRegionIdx,
                                const double* p,
                                const double* z,
                                double* output_R,
                                double* output_dRdp) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            evalRDeriv(p[i], z + num_phases_*i, getTableIndex_(pvtRegionIdx, i), output_R[i], output_dRdp[i]);
        }
    }


    // ---- Private methods ----

    double PvtLiveGas::evalB(const double press, const double* surfvol, int pvtTableIdx) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            // To handle no-gas case.
            return 1.0;
        }
        return 1.0/miscible_gas(press, surfvol, pvtTableIdx, 1, false);
    }

    void PvtLiveGas::evalBDeriv(const double press, const double* surfvol, int pvtTableIdx,
                                double& Bval, double& dBdpval) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            // To handle no-gas case.
            Bval = 1.0;
            dBdpval = 0.0;
            return;
        }
        Bval = evalB(press, surfvol, pvtTableIdx);
        dBdpval =  -Bval*Bval*miscible_gas(press, surfvol, pvtTableIdx, 1, true);
    }

    double PvtLiveGas::evalR(const double press, const double* surfvol, int pvtTableIdx) const
    {
        if (surfvol[phase_pos_[Liquid]] == 0.0) {
            // To handle no-gas case.
            return 0.0;
        }
        double satR = linearInterpolation(saturated_gas_table_[pvtTableIdx][0],
                                             saturated_gas_table_[pvtTableIdx][4], press);
        double maxR = surfvol[phase_pos_[Liquid]]/surfvol[phase_pos_[Vapour]];
        if (satR < maxR ) {
            // Saturated case
            return satR;
        } else {
            // Undersaturated case
            return maxR;
        }
    }

    void PvtLiveGas::evalRDeriv(const double press, const double* surfvol, int pvtTableIdx,
                                double& Rval, double& dRdpval) const
    {
        if (surfvol[phase_pos_[Liquid]] == 0.0) {
            // To handle no-gas case.
            Rval = 0.0;
            dRdpval = 0.0;
            return;
        }
        double satR = linearInterpolation(saturated_gas_table_[pvtTableIdx][0],
                                             saturated_gas_table_[pvtTableIdx][4], press);
        double maxR = surfvol[phase_pos_[Liquid]]/surfvol[phase_pos_[Vapour]];
        if (satR < maxR ) {
            // Saturated case
            Rval = satR;
            dRdpval = linearInterpolationDerivative(saturated_gas_table_[pvtTableIdx][0],
                                            saturated_gas_table_[pvtTableIdx][4],
                                            press);
        } else {
            // Undersaturated case
            Rval = maxR;
            dRdpval = 0.0;
        }
    }

    double PvtLiveGas::miscible_gas(const double press,
                                    const double* surfvol,
                                    const int pvtTableIdx,
                                    const int item,
                                    const bool deriv) const
    {
        const std::vector<std::vector<double> > &saturatedGasTable =
            saturated_gas_table_[pvtTableIdx];
        const std::vector<std::vector<std::vector<double> > > &undersatGasTables =
            undersat_gas_tables_[pvtTableIdx];

        int section;
        double Rval = linearInterpolation(saturatedGasTable[0],
                                                saturatedGasTable[4], press,
                                                section);
        double maxR = surfvol[phase_pos_[Liquid]]/surfvol[phase_pos_[Vapour]];
        if (deriv) {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolationDerivative(saturatedGasTable[0],
                                                saturatedGasTable[item],
                                                press);
            } else {  // Undersaturated case
                int is = section;
                if (undersatGasTables[is][0].size() < 2) {
                    double val = (saturatedGasTable[item][is+1]
                                  - saturatedGasTable[item][is]) /
                        (saturatedGasTable[0][is+1] -
                         saturatedGasTable[0][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersatGasTables[is][0],
                                              undersatGasTables[is][item],
                                              maxR);
                double val2 =
                    linearInterpolation(undersatGasTables[is+1][0],
                                              undersatGasTables[is+1][item],
                                              maxR);
                double val = (val2 - val1)/
                    (saturatedGasTable[0][is+1] - saturatedGasTable[0][is]);
                return val;
            }
        } else {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolation(saturatedGasTable[0],
                                                 saturatedGasTable[item],
                                                 press);
            } else {  // Undersaturated case
                int is = section;
                // Extrapolate from first table section
                if (is == 0 && press < saturatedGasTable[0][0]) {
                    return linearInterpolation(undersatGasTables[0][0],
                                                     undersatGasTables[0][item],
                                                     maxR);
                }

                // Extrapolate from last table section
                int ltp = saturatedGasTable[0].size() - 1;
                if (is+1 == ltp && press > saturatedGasTable[0][ltp]) {
                    return linearInterpolation(undersatGasTables[ltp][0],
                                                     undersatGasTables[ltp][item],
                                                     maxR);
                }

                // Interpolate between table sections
                double w = (press - saturatedGasTable[0][is]) /
                    (saturatedGasTable[0][is+1] -
                     saturatedGasTable[0][is]);
                if (undersatGasTables[is][0].size() < 2) {
                    double val = saturatedGasTable[item][is] +
                        w*(saturatedGasTable[item][is+1] -
                           saturatedGasTable[item][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersatGasTables[is][0],
                                              undersatGasTables[is][item],
                                              maxR);
                double val2 =
                    linearInterpolation(undersatGasTables[is+1][0],
                                              undersatGasTables[is+1][item],
                                              maxR);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }

    double PvtLiveGas::miscible_gas(const double press,
                                    const double r,
                                    const PhasePresence& cond,
                                    const int pvtTableIdx,
                                    const int item,
                                    const int deriv) const
    {
        const std::vector<std::vector<double> > &saturatedGasTable =
            saturated_gas_table_[pvtTableIdx];
        const std::vector<std::vector<std::vector<double> > > &undersatGasTables =
            undersat_gas_tables_[pvtTableIdx];

        const bool isSat = cond.hasFreeOil();

        // Derivative w.r.t p
        if (deriv == 1) {
            if (isSat) {  // Saturated case
                return linearInterpolationDerivative(saturatedGasTable[0],
                                                saturatedGasTable[item],
                                                press);
            } else {  // Undersaturated case
                int is = tableIndex(saturatedGasTable[0], press);
                if (undersatGasTables[is][0].size() < 2) {
                    double val = (saturatedGasTable[item][is+1]
                                  - saturatedGasTable[item][is]) /
                        (saturatedGasTable[0][is+1] -
                         saturatedGasTable[0][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersatGasTables[is][0],
                                              undersatGasTables[is][item],
                                              r);
                double val2 =
                    linearInterpolation(undersatGasTables[is+1][0],
                                              undersatGasTables[is+1][item],
                                              r);
                double val = (val2 - val1)/
                    (saturatedGasTable[0][is+1] - saturatedGasTable[0][is]);
                return val;
            }
        } else if (deriv == 2){
            if (isSat) {
                return 0;
            } else {
                int is = tableIndex(saturatedGasTable[0], press);
                double w = (press - saturatedGasTable[0][is]) /
                    (saturatedGasTable[0][is+1] - saturatedGasTable[0][is]);
                assert(undersatGasTables[is][0].size() >= 2);
                assert(undersatGasTables[is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersatGasTables[is][0],
                                              undersatGasTables[is][item],
                                              r);
                double val2 =
                    linearInterpolationDerivative(undersatGasTables[is+1][0],
                                              undersatGasTables[is+1][item],
                                              r);

                double val = val1 + w * (val2 - val1);
                return val;

            }
        } else {
            if (isSat) {  // Saturated case
                return linearInterpolation(saturatedGasTable[0],
                                                 saturatedGasTable[item],
                                                 press);
            } else {  // Undersaturated case
                int is = tableIndex(saturatedGasTable[0], press);
                // Extrapolate from first table section
                if (is == 0 && press < saturatedGasTable[0][0]) {
                    return linearInterpolation(undersatGasTables[0][0],
                                                     undersatGasTables[0][item],
                                                     r);
                }

                // Extrapolate from last table section
                //int ltp = saturatedGasTable[0].size() - 1;
                //if (is+1 == ltp && press > saturatedGasTable[0][ltp]) {
                //    return linearInterpolation(undersatGasTables[ltp][0],
                //                                    undersatGasTables[ltp][item],
                //                                    r);
                //}

                // Interpolate between table sections
                double w = (press - saturatedGasTable[0][is]) /
                    (saturatedGasTable[0][is+1] -
                     saturatedGasTable[0][is]);
                if (undersatGasTables[is][0].size() < 2) {
                    double val = saturatedGasTable[item][is] +
                        w*(saturatedGasTable[item][is+1] -
                           saturatedGasTable[item][is]);
                    return val;
                }
                double val1 =
                    linearInterpolation(undersatGasTables[is][0],
                                              undersatGasTables[is][item],
                                              r);
                double val2 =
                    linearInterpolation(undersatGasTables[is+1][0],
                                              undersatGasTables[is+1][item],
                                              r);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }




} // namespace Opm
