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
#include <opm/core/props/pvt/PvtLiveOil.hpp>
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
    PvtLiveOil::PvtLiveOil(const std::vector<Opm::PvtoTable>& pvtoTables)
    {
        int numTables = pvtoTables.size();
        saturated_oil_table_.resize(numTables);
        undersat_oil_tables_.resize(numTables);

        for (int pvtTableIdx = 0; pvtTableIdx < numTables; ++pvtTableIdx) {
            Opm::PvtoTable pvtoTable = pvtoTables[pvtTableIdx];

            const auto saturatedPvto = pvtoTable.getOuterTable();

            // OIL, PVTO
            // saturated_oil_table_[pvtTableIdx].resize(4);
            // adding a extra colummn to the PVTO to store 1/(B*mu)
            saturated_oil_table_[pvtTableIdx].resize(5);
            const int sz = saturatedPvto->numRows();
            // for (int k=0; k<4; ++k) {
            for (int k=0; k<5; ++k) {
                saturated_oil_table_[pvtTableIdx][k].resize(sz);
            }
            for (int i=0; i<sz; ++i) {
                saturated_oil_table_[pvtTableIdx][0][i] = saturatedPvto->getPressureColumn()[i]; // p
                saturated_oil_table_[pvtTableIdx][1][i] = 1.0/saturatedPvto->getOilFormationFactorColumn()[i]; // 1/Bo
                saturated_oil_table_[pvtTableIdx][2][i] = saturatedPvto->getOilViscosityColumn()[i]; // mu_o
                saturated_oil_table_[pvtTableIdx][3][i] = 1.0 / (saturatedPvto->getOilFormationFactorColumn()[i]
                                                               * saturatedPvto->getOilViscosityColumn()[i]); // 1/(Bo*mu_o)
                saturated_oil_table_[pvtTableIdx][4][i] = saturatedPvto->getGasSolubilityColumn()[i]; // Rs
            }

            undersat_oil_tables_[pvtTableIdx].resize(sz);
            for (int i=0; i<sz; ++i) {
                const auto undersaturatedPvto = pvtoTable.getInnerTable(i);

                // undersat_oil_tables_[pvtTableIdx][i].resize(3);
                // adding a extra colummn to the PVTO to store 1/(B*mu)
                undersat_oil_tables_[pvtTableIdx][i].resize(4);
                int tsize = undersaturatedPvto->numRows();
                undersat_oil_tables_[pvtTableIdx][i][0].resize(tsize);
                undersat_oil_tables_[pvtTableIdx][i][1].resize(tsize);
                undersat_oil_tables_[pvtTableIdx][i][2].resize(tsize);
                undersat_oil_tables_[pvtTableIdx][i][3].resize(tsize);
                for (int j=0; j<tsize; ++j) {
                    undersat_oil_tables_[pvtTableIdx][i][0][j] = undersaturatedPvto->getPressureColumn()[j];  // p
                    undersat_oil_tables_[pvtTableIdx][i][1][j] = 1.0/undersaturatedPvto->getOilFormationFactorColumn()[j];  // 1/Bo
                    undersat_oil_tables_[pvtTableIdx][i][2][j] = undersaturatedPvto->getOilViscosityColumn()[j];  // mu_o
                    undersat_oil_tables_[pvtTableIdx][i][3][j] = 1.0 / (undersaturatedPvto->getOilFormationFactorColumn()[j] *
                                                                        undersaturatedPvto->getOilViscosityColumn()[j]);  // 1/(Bo*mu_o)
                }
            }

            // Complete undersaturated tables by extrapolating from existing data
            // as is done in Eclipse and Mrst
            // TODO: check if the following formulations applying to 1/(Bo*mu_o)
            int iNext = -1;
            for (int i=0; i<sz; ++i) {
                // Skip records already containing undersaturated data
                if (undersat_oil_tables_[pvtTableIdx][i][0].size() > 1) {
                    continue;
                }
                // Look ahead for next record containing undersaturated data
                if (iNext < i) {
                    iNext = i+1;
                    while (iNext<sz && undersat_oil_tables_[pvtTableIdx][iNext][0].size() < 2) {
                        ++iNext;
                    }
                    if (iNext == sz) OPM_THROW(std::runtime_error,"Unable to complete undersaturated table.");
                }
                // Add undersaturated data to current record while maintaining compressibility and viscosibility
                // TODO: How to add 1/(B*mu) in this way?
                typedef std::vector<std::vector<std::vector<double> > >::size_type sz_t;
                for (sz_t j=1; j<undersat_oil_tables_[pvtTableIdx][iNext][0].size(); ++j) {
                    double diffPressure = undersat_oil_tables_[pvtTableIdx][iNext][0][j]-undersat_oil_tables_[pvtTableIdx][iNext][0][j-1];
                    double pressure = undersat_oil_tables_[pvtTableIdx][i][0].back()+diffPressure;
                    undersat_oil_tables_[pvtTableIdx][i][0].push_back(pressure);
                    double compr = (1.0/undersat_oil_tables_[pvtTableIdx][iNext][1][j]-1.0/undersat_oil_tables_[pvtTableIdx][iNext][1][j-1])
                        / (0.5*(1.0/undersat_oil_tables_[pvtTableIdx][iNext][1][j]+1.0/undersat_oil_tables_[pvtTableIdx][iNext][1][j-1]));
                    double B = (1.0/undersat_oil_tables_[pvtTableIdx][i][1].back())*(1.0+0.5*compr)/(1.0-0.5*compr);
                    undersat_oil_tables_[pvtTableIdx][i][1].push_back(1.0/B);
                    double visc = (undersat_oil_tables_[pvtTableIdx][iNext][2][j]-undersat_oil_tables_[pvtTableIdx][iNext][2][j-1])
                        / (0.5*(undersat_oil_tables_[pvtTableIdx][iNext][2][j]+undersat_oil_tables_[pvtTableIdx][iNext][2][j-1]));
                    double mu = (undersat_oil_tables_[pvtTableIdx][i][2].back())*(1.0+0.5*visc)/(1.0-0.5*visc);
                    undersat_oil_tables_[pvtTableIdx][i][2].push_back(mu);

                    // A try to expolate the 1/BMu with the expolated mu and B
                    double inverseBMu = 1.0 / (B*mu);
                    undersat_oil_tables_[pvtTableIdx][i][3].push_back(inverseBMu);

                }
            }
        }
    }

    /// Destructor.
    PvtLiveOil::~PvtLiveOil()
    {
    }


    /// Viscosity as a function of p, T and z.
    void PvtLiveOil::mu(const int n,
                        const int* pvtTableIdx,
                        const double* p,
                        const double* /*T*/,
                              const double* z,
                              double* output_mu) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int tableIdx = getTableIndex_(pvtTableIdx, i);

            double inverseB = miscible_oil(p[i], z + num_phases_*i, tableIdx, 1, false);
            double inverseBMu = miscible_oil(p[i], z + num_phases_*i, tableIdx, 3, false);

            output_mu[i] = inverseB / inverseBMu; 
        }
    }

    /// Viscosity and its p and r derivatives as a function of p, T and r.
    void PvtLiveOil::mu(const int n,
                        const int* pvtTableIdx,
                        const double* p,
                        const double* /*T*/,
                               const double* r,
                               double* output_mu,
                               double* output_dmudp,
                               double* output_dmudr) const
    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    int tableIdx = getTableIndex_(pvtTableIdx, i);

                    double inverseBMu = miscible_oil(p[i], r[i], tableIdx, 3, 0);
                    double inverseB = miscible_oil(p[i], r[i], tableIdx, 1, 0);

                    output_mu[i] = inverseB / inverseBMu;

                    double dinverseBmudp = miscible_oil(p[i], r[i], tableIdx, 3, 1);
                    double dinverseBdp = miscible_oil(p[i], r[i], tableIdx, 1, 1);

                    output_dmudp[i] = (inverseBMu * dinverseBdp - inverseB * dinverseBmudp)
                                      / (inverseBMu * inverseBMu);


                    double dinverseBmudr = miscible_oil(p[i], r[i], tableIdx, 3, 2);
                    double dinverseBdr = miscible_oil(p[i], r[i], tableIdx, 1, 2);

                    output_dmudr[i] = (inverseBMu * dinverseBdr - inverseB * dinverseBmudr)
                                      / (inverseBMu * inverseBMu);

                }
    }

    /// Viscosity and its p and r derivatives as a function of p, T and r.
    void PvtLiveOil::mu(const int n,
                        const int* pvtTableIdx,
                        const double* p,
                        const double* /*T*/,
                               const double* r,
                               const PhasePresence* cond,
                               double* output_mu,
                               double* output_dmudp,
                               double* output_dmudr) const
    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    int tableIdx = getTableIndex_(pvtTableIdx, i);
                    const PhasePresence& cnd = cond[i];

                    double inverseBMu = miscible_oil(p[i], r[i], cnd, tableIdx, 3, 0);
                    double inverseB = miscible_oil(p[i], r[i], cnd, tableIdx, 1, 0);

                    output_mu[i] = inverseB / inverseBMu;

                    double dinverseBmudp = miscible_oil(p[i], r[i], cnd, tableIdx, 3, 1);
                    double dinverseBdp = miscible_oil(p[i], r[i], cnd, tableIdx, 1, 1);

                    output_dmudp[i] = (inverseBMu * dinverseBdp - inverseB * dinverseBmudp)
                                      / (inverseBMu * inverseBMu);


                    double dinverseBmudr = miscible_oil(p[i], r[i], cnd, tableIdx, 3, 2);
                    double dinverseBdr = miscible_oil(p[i], r[i], cnd, tableIdx, 1, 2);

                    output_dmudr[i] = (inverseBMu * dinverseBdr - inverseB * dinverseBmudr)
                                      / (inverseBMu * inverseBMu);

                }
    }


    /// Formation volume factor as a function of p, T and z.
    void PvtLiveOil::B(const int n,
                       const int* pvtTableIdx,
                       const double* p,
                       const double* /*T*/,
                             const double* z,
                             double* output_B) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int tableIdx = getTableIndex_(pvtTableIdx, i);

            output_B[i] = evalB(tableIdx, p[i], z + num_phases_*i);
        }

    }


    /// Formation volume factor and p-derivative as functions of p, T and z.
    void PvtLiveOil::dBdp(const int n,
                          const int* pvtTableIdx,
                          const double* p,
                          const double* /*T*/,
                                const double* z,
                                double* output_B,
                                double* output_dBdp) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int tableIdx = getTableIndex_(pvtTableIdx, i);

            evalBDeriv(tableIdx, p[i], z + num_phases_*i, output_B[i], output_dBdp[i]);
        }
    }

    void PvtLiveOil::b(const int n,
                       const int* pvtTableIdx,
                       const double* p,
                       const double* /*T*/,
                          const double* r,
                          double* output_b,
                          double* output_dbdp,
                          double* output_dbdr) const

    {
        // #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    int tableIdx = getTableIndex_(pvtTableIdx, i);

                    output_b[i] = miscible_oil(p[i], r[i], tableIdx, 1, 0);
                    output_dbdp[i] = miscible_oil(p[i], r[i], tableIdx, 1, 1);
                    output_dbdr[i] = miscible_oil(p[i], r[i], tableIdx, 1, 2);

                }
    }

    void PvtLiveOil::b(const int n,
                       const int* pvtTableIdx,
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
                    int tableIdx = getTableIndex_(pvtTableIdx, i);

                    output_b[i] = miscible_oil(p[i], r[i], cnd, tableIdx, 1, 0);
                    output_dbdp[i] = miscible_oil(p[i], r[i], cnd, tableIdx, 1, 1);
                    output_dbdr[i] = miscible_oil(p[i], r[i], cnd, tableIdx, 1, 2);

                }
    }

    void PvtLiveOil::rsSat(const int n,
                           const int* pvtTableIdx,
                             const double* p,
                             double* output_rsSat,
                             double* output_drsSatdp) const
    {

        for (int i = 0; i < n; ++i) {
            int tableIdx = getTableIndex_(pvtTableIdx, i);
            output_rsSat[i] = linearInterpolation(saturated_oil_table_[tableIdx][0],
                    saturated_oil_table_[tableIdx][4],p[i]);
            output_drsSatdp[i] = linearInterpolationDerivative(saturated_oil_table_[tableIdx][0],
                    saturated_oil_table_[tableIdx][4],p[i]);

        }
    }

    void PvtLiveOil::rvSat(const int n,
                           const int* /*pvtTableIdx*/,
                             const double* /*p*/,
                             double* output_rvSat,
                             double* output_drvSatdp) const
    {
        std::fill(output_rvSat, output_rvSat + n, 0.0);
        std::fill(output_drvSatdp, output_drvSatdp + n, 0.0);
    }

    /// Solution factor as a function of p and z.
    void PvtLiveOil::R(const int n,
                       const int* pvtTableIdx,
                             const double* p,
                             const double* z,
                             double* output_R) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int tableIdx = getTableIndex_(pvtTableIdx, i);
            output_R[i] = evalR(tableIdx, p[i], z + num_phases_*i);
        }

    }


    /// Solution factor and p-derivative as functions of p and z.
    void PvtLiveOil::dRdp(const int n,
                          const int* pvtTableIdx,
                                const double* p,
                                const double* z,
                                double* output_R,
                                double* output_dRdp) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int tableIdx = getTableIndex_(pvtTableIdx, i);
            evalRDeriv(tableIdx, p[i], z + num_phases_*i, output_R[i], output_dRdp[i]);
        }
    }




    // ---- Private methods ----

    double PvtLiveOil::evalB(size_t pvtTableIdx,
                                   double press, const double* surfvol) const
    {
        // if (surfvol[phase_pos_[Liquid]] == 0.0) return 1.0; // To handle no-oil case.
        return 1.0/miscible_oil(press, surfvol, pvtTableIdx, 1, false);
    }


    void PvtLiveOil::evalBDeriv(size_t pvtTableIdx,
                                      const double press, const double* surfvol,
                                      double& Bval, double& dBdpval) const
    {
        Bval = evalB(pvtTableIdx, press, surfvol);
        dBdpval = -Bval*Bval*miscible_oil(press, surfvol, pvtTableIdx, 1, true);
    }

    double PvtLiveOil::evalR(size_t pvtTableIdx,
                                   double press, const double* surfvol) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            return 0.0;
        }
        double Rval = linearInterpolation(saturated_oil_table_[pvtTableIdx][0],
                                             saturated_oil_table_[pvtTableIdx][4], press);
        double maxR = surfvol[phase_pos_[Vapour]]/surfvol[phase_pos_[Liquid]];
        if (Rval < maxR ) {  // Saturated case
            return Rval;
        } else {
            return maxR;  // Undersaturated case
        }
    }

    void PvtLiveOil::evalRDeriv(size_t pvtTableIdx,
                                      const double press, const double* surfvol,
                                      double& Rval, double& dRdpval) const
    {
        if (surfvol[phase_pos_[Vapour]] == 0.0) {
            Rval = 0.0;
            dRdpval = 0.0;
            return;
        }
        Rval = linearInterpolation(saturated_oil_table_[pvtTableIdx][0],
                                   saturated_oil_table_[pvtTableIdx][4], press);
        double maxR = surfvol[phase_pos_[Vapour]]/surfvol[phase_pos_[Liquid]];
        if (Rval < maxR ) {
            // Saturated case
            dRdpval = linearInterpolationDerivative(saturated_oil_table_[pvtTableIdx][0],
                                                    saturated_oil_table_[pvtTableIdx][4],
                                                    press);
        } else {
            // Undersaturated case
            Rval = maxR;
            dRdpval = 0.0;
        }
    }


    // TODO: Check if this function need to be adapted to the 1/(B*mu) interpolation.
    double PvtLiveOil::miscible_oil(const double press,
                                    const double* surfvol,
                                    const int pvtTableIdx,
                                    const int item,
                                    const bool deriv) const
    {
        int section;
        double Rval = linearInterpolation(saturated_oil_table_[pvtTableIdx][0],
                                          saturated_oil_table_[pvtTableIdx][4],
                                          press, section);
        double maxR = (surfvol[phase_pos_[Liquid]] == 0.0) ? 0.0 : surfvol[phase_pos_[Vapour]]/surfvol[phase_pos_[Liquid]];
        if (deriv) {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolationDerivative(saturated_oil_table_[pvtTableIdx][0],
                                                     saturated_oil_table_[pvtTableIdx][item],
                                                     press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], maxR);
                double w = (maxR - saturated_oil_table_[pvtTableIdx][4][is]) /
                    (saturated_oil_table_[pvtTableIdx][4][is+1] - saturated_oil_table_[pvtTableIdx][4][is]);
                assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersat_oil_tables_[pvtTableIdx][is][0],
                                                  undersat_oil_tables_[pvtTableIdx][is][item],
                                                  press);
                double val2 =
                    linearInterpolationDerivative(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                                  undersat_oil_tables_[pvtTableIdx][is+1][item],
                                                  press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        } else {
            if (Rval < maxR ) {  // Saturated case
                return linearInterpolation(saturated_oil_table_[pvtTableIdx][0],
                                           saturated_oil_table_[pvtTableIdx][item],
                                           press);
            } else {  // Undersaturated case
                // Interpolate between table sections
                int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], maxR);
                double w = (maxR - saturated_oil_table_[pvtTableIdx][4][is]) /
                    (saturated_oil_table_[pvtTableIdx][4][is+1] - saturated_oil_table_[pvtTableIdx][4][is]);
                assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is][0],
                                        undersat_oil_tables_[pvtTableIdx][is][item],
                                        press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                        undersat_oil_tables_[pvtTableIdx][is+1][item],
                                        press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }

    double PvtLiveOil::miscible_oil(const double press,
                                    const double r,
                                    const int pvtTableIdx,
                                    const int item,
                                    const int deriv) const
    {
        int section;
        double Rval = linearInterpolation(saturated_oil_table_[pvtTableIdx][0],
                                          saturated_oil_table_[pvtTableIdx][4],
                                          press, section);
        // derivative with respect to frist component (pressure)
        if (deriv == 1) {
            if (Rval <= r ) {  // Saturated case
                return linearInterpolationDerivative(saturated_oil_table_[pvtTableIdx][0],
                                                     saturated_oil_table_[pvtTableIdx][item],
                                                     press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], r);
                double w = (r - saturated_oil_table_[pvtTableIdx][4][is]) /
                    (saturated_oil_table_[pvtTableIdx][4][is+1] - saturated_oil_table_[pvtTableIdx][4][is]);
                assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersat_oil_tables_[pvtTableIdx][is][0],
                                             undersat_oil_tables_[pvtTableIdx][is][item],
                                             press);
                double val2 =
                    linearInterpolationDerivative(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                             undersat_oil_tables_[pvtTableIdx][is+1][item],
                                             press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
            // derivative with respect to second component (r)
        } else if (deriv == 2)  {
            if (Rval <= r ) {  // Saturated case
                return 0;
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], r);
                assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is][0],
                                              undersat_oil_tables_[pvtTableIdx][is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                              undersat_oil_tables_[pvtTableIdx][is+1][item],
                                              press);

                double val = (val2 - val1)/(saturated_oil_table_[pvtTableIdx][4][is+1]-saturated_oil_table_[pvtTableIdx][4][is]);
                return val;
            }


            }    else {
            if (Rval <= r ) {  // Saturated case
                return linearInterpolation(saturated_oil_table_[pvtTableIdx][0],
                                           saturated_oil_table_[pvtTableIdx][item],
                                           press);
            } else {  // Undersaturated case
                // Interpolate between table sections
                int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], r);
                double w = (r - saturated_oil_table_[pvtTableIdx][4][is]) /
                    (saturated_oil_table_[pvtTableIdx][4][is+1] - saturated_oil_table_[pvtTableIdx][4][is]);
                assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is][0],
                                              undersat_oil_tables_[pvtTableIdx][is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                              undersat_oil_tables_[pvtTableIdx][is+1][item],
                                              press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
        }
    }

    double PvtLiveOil::miscible_oil(const double press,
                                    const double r,
                                    const PhasePresence& cond,
                                    const int pvtTableIdx,
                                    const int item,
                                    const int deriv) const
    {
        const bool isSat = cond.hasFreeGas();

        // derivative with respect to frist component (pressure)
        if (deriv == 1) {
            if (isSat) {  // Saturated case
                return linearInterpolationDerivative(saturated_oil_table_[pvtTableIdx][0],
                                                     saturated_oil_table_[pvtTableIdx][item],
                                                     press);
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], r);
                double w = (r - saturated_oil_table_[pvtTableIdx][4][is]) /
                    (saturated_oil_table_[pvtTableIdx][4][is+1] - saturated_oil_table_[pvtTableIdx][4][is]);
                assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                double val1 =
                    linearInterpolationDerivative(undersat_oil_tables_[pvtTableIdx][is][0],
                                             undersat_oil_tables_[pvtTableIdx][is][item],
                                             press);
                double val2 =
                    linearInterpolationDerivative(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                             undersat_oil_tables_[pvtTableIdx][is+1][item],
                                             press);
                double val = val1 + w*(val2 - val1);
                return val;
            }
            // derivative with respect to second component (r)
        } else if (deriv == 2)  {
            if (isSat) {  // Saturated case
                return 0;
            } else {  // Undersaturated case
                int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], r);
                assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                double val1 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is][0],
                                              undersat_oil_tables_[pvtTableIdx][is][item],
                                              press);
                double val2 =
                    linearInterpolation(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                              undersat_oil_tables_[pvtTableIdx][is+1][item],
                                              press);

                double val = (val2 - val1)/(saturated_oil_table_[pvtTableIdx][4][is+1]-saturated_oil_table_[pvtTableIdx][4][is]);
                return val;
            }


            }
            // no derivative
            else {
                if (isSat) {  // Saturated case
                    return linearInterpolation(saturated_oil_table_[pvtTableIdx][0],
                                               saturated_oil_table_[pvtTableIdx][item],
                                               press);
                } else {  // Undersaturated case
                    // Interpolate between table sections
                    int is = tableIndex(saturated_oil_table_[pvtTableIdx][4], r);
                    double w = (r - saturated_oil_table_[pvtTableIdx][4][is]) /
                        (saturated_oil_table_[pvtTableIdx][4][is+1] - saturated_oil_table_[pvtTableIdx][4][is]);
                    assert(undersat_oil_tables_[pvtTableIdx][is][0].size() >= 2);
                    assert(undersat_oil_tables_[pvtTableIdx][is+1][0].size() >= 2);
                    double val1 =
                        linearInterpolation(undersat_oil_tables_[pvtTableIdx][is][0],
                                                  undersat_oil_tables_[pvtTableIdx][is][item],
                                                  press);
                    double val2 =
                        linearInterpolation(undersat_oil_tables_[pvtTableIdx][is+1][0],
                                                  undersat_oil_tables_[pvtTableIdx][is+1][item],
                                                  press);
                    double val = val1 + w*(val2 - val1);
                    return val;
                }
            }
        }

} // namespace Opm
