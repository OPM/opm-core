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
#include <opm/core/props/pvt/PvtDeadSpline.hpp>
#include <opm/core/utility/buildUniformMonotoneTable.hpp>


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

    PvtDeadSpline::PvtDeadSpline()
    {}

    void PvtDeadSpline::initFromOil(const std::vector<Opm::PvdoTable>& pvdoTables,
                                    int numSamples)
    {
        int numRegions = pvdoTables.size();

        // resize the attributes of the object
        b_.resize(numRegions);
        viscosity_.resize(numRegions);

        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const Opm::PvdoTable& pvdoTable = pvdoTables[regionIdx];

            int numRows = pvdoTable.numRows();

            // Copy data
            const std::vector<double> &press = pvdoTable.getPressureColumn();
            const std::vector<double> &B = pvdoTable.getFormationFactorColumn();
            const std::vector<double> &visc = pvdoTable.getViscosityColumn();

            std::vector<double> B_inv(numRows);
            for (int i = 0; i < numRows; ++i) {
                B_inv[i] = 1.0 / B[i];
            }

            buildUniformMonotoneTable(press, B_inv, numSamples, b_[regionIdx]);
            buildUniformMonotoneTable(press, visc, numSamples, viscosity_[regionIdx]);
        }
    }

    void PvtDeadSpline::initFromGas(const std::vector<Opm::PvdgTable>& pvdgTables,
                                    int numSamples)
    {
        int numRegions = pvdgTables.size();

        // resize the attributes of the object
        b_.resize(numRegions);
        viscosity_.resize(numRegions);

        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const Opm::PvdgTable& pvdgTable = pvdgTables[regionIdx];

            int numRows = pvdgTable.numRows();

            // Copy data
            const std::vector<double> &press = pvdgTable.getPressureColumn();
            const std::vector<double> &B = pvdgTable.getFormationFactorColumn();
            const std::vector<double> &visc = pvdgTable.getViscosityColumn();

            std::vector<double> B_inv(numRows);
            for (int i = 0; i < numRows; ++i) {
                B_inv[i] = 1.0 / B[i];
            }

            buildUniformMonotoneTable(press, B_inv, numSamples, b_[regionIdx]);
            buildUniformMonotoneTable(press, visc, numSamples, viscosity_[regionIdx]);
        }
    }

    // Destructor
    PvtDeadSpline::~PvtDeadSpline()
    {
    }



    void PvtDeadSpline::mu(const int n,
                           const int* pvtTableIdx,
                           const double* p,
                           const double* /*T*/,
                                 const double* /*z*/,
                                 double* output_mu) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int regionIdx = getTableIndex_(pvtTableIdx, i);
            output_mu[i] = viscosity_[regionIdx](p[i]);
        }
    }

    void PvtDeadSpline::mu(const int n,
                           const int* pvtTableIdx,
                           const double* p,
                           const double* /*T*/,
                           const double* /*r*/,
                           double* output_mu,
                           double* output_dmudp,
                           double* output_dmudr) const
    {
        // #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int regionIdx = getTableIndex_(pvtTableIdx, i);
            output_mu[i] = viscosity_[regionIdx](p[i]);
            output_dmudp[i] = viscosity_[regionIdx].derivative(p[i]);
        }
        std::fill(output_dmudr, output_dmudr + n, 0.0);
    }

    void PvtDeadSpline::mu(const int n,
                           const int* pvtTableIdx,
                           const double* p,
                           const double* /*T*/,
                           const double* /*r*/,
                           const PhasePresence* /*cond*/,
                           double* output_mu,
                           double* output_dmudp,
                           double* output_dmudr) const
    {
    // #pragma omp parallel for

        for (int i = 0; i < n; ++i) {
            int regionIdx = getTableIndex_(pvtTableIdx, i);
            output_mu[i] = viscosity_[regionIdx](p[i]);
            output_dmudp[i] = viscosity_[regionIdx].derivative(p[i]);
        }
        std::fill(output_dmudr, output_dmudr + n, 0.0);
    }

    void PvtDeadSpline::B(const int n,
                          const int* pvtTableIdx,
                          const double* p,
                          const double* /*T*/,
                                const double* /*z*/,
                                double* output_B) const
    {
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int regionIdx = getTableIndex_(pvtTableIdx, i);
            output_B[i] = 1.0/b_[regionIdx](p[i]);
        }
    }

    void PvtDeadSpline::dBdp(const int n,
                             const int* pvtTableIdx,
                             const double* p,
                             const double* T,
                             const double* /*z*/,
                             double* output_B,
                             double* output_dBdp) const
    {
        B(n, pvtTableIdx, p, T, 0, output_B);
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int regionIdx = getTableIndex_(pvtTableIdx, i);
            double Bg = output_B[i];
            output_dBdp[i] = -Bg*Bg*b_[regionIdx].derivative(p[i]);
        }
    }

    void PvtDeadSpline::b(const int n,
                          const int* pvtTableIdx,
                          const double* p,
                          const double* T,
                              const double* /*r*/,
                              double* output_b,
                              double* output_dbdp,
                              double* output_dbdr) const
    {
    // #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int regionIdx = getTableIndex_(pvtTableIdx, i);
            output_b[i] = b_[regionIdx](p[i]);
            output_dbdp[i] = b_[regionIdx].derivative(p[i]);
        }
        std::fill(output_dbdr, output_dbdr + n, 0.0);
    }

    void PvtDeadSpline::b(const int n,
                          const int* pvtTableIdx,
                          const double* p,
                          const double* T,
                          const double* /*r*/,
                          const PhasePresence* /*cond*/,
                          double* output_b,
                          double* output_dbdp,
                          double* output_dbdr) const
    {
    // #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int regionIdx = getTableIndex_(pvtTableIdx, i);
            output_b[i] = b_[regionIdx](p[i]);
            output_dbdp[i] = b_[regionIdx].derivative(p[i]);
        }
        std::fill(output_dbdr, output_dbdr + n, 0.0);
    }

    void PvtDeadSpline::rsSat(const int n,
                              const int* /*pvtTableIdx*/,
                             const double* /*p*/,
                             double* output_rsSat,
                             double* output_drsSatdp) const
    {
        std::fill(output_rsSat, output_rsSat + n, 0.0);
        std::fill(output_drsSatdp, output_drsSatdp + n, 0.0);
    }

    void PvtDeadSpline::rvSat(const int n,
                              const int* /*pvtTableIdx*/,
                             const double* /*p*/,
                             double* output_rvSat,
                             double* output_drvSatdp) const
    {
        std::fill(output_rvSat, output_rvSat + n, 0.0);
        std::fill(output_drvSatdp, output_drvSatdp + n, 0.0);
    }

    void PvtDeadSpline::R(const int n,
                          const int* /*pvtTableIdx*/,
                          const double* /*p*/,
                          const double* /*z*/,
                          double* output_R) const
    {
        std::fill(output_R, output_R + n, 0.0);
    }

    void PvtDeadSpline::dRdp(const int n,
                             const int* /*pvtTableIdx*/,
                             const double* /*p*/,
                             const double* /*z*/,
                             double* output_R,
                             double* output_dRdp) const
    {
        std::fill(output_R, output_R + n, 0.0);
        std::fill(output_dRdp, output_dRdp + n, 0.0);
    }

}
