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

#include <opm/core/fluid/SaturationPropsFromDeck.hpp>
#include <opm/core/fluid/blackoil/phaseUsageFromDeck.hpp>
#include <opm/core/utility/buildUniformMonotoneTable.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{

    /// Default constructor.
    SaturationPropsFromDeck::SaturationPropsFromDeck()
    {
    }

    /// Initialize from deck.
    void SaturationPropsFromDeck::init(const EclipseGridParser& deck,
                                       const std::vector<int>& global_cell)
    {
        phase_usage_ = phaseUsageFromDeck(deck);

        // Extract input data.
        // Oil phase should be active.
        if (!phase_usage_.phase_used[Liquid]) {
            THROW("SaturationPropsFromDeck::init()   --  oil phase must be active.");
        }

        // Obtain SATNUM, if it exists, and create cell_to_func_.
        // Otherwise, let the cell_to_func_ mapping be just empty.
        int satfuncs_expected = 1;
        if (deck.hasField("SATNUM")) {
            const std::vector<int>& satnum = deck.getIntegerValue("SATNUM");
            satfuncs_expected = *std::max_element(satnum.begin(), satnum.end());
            int num_cells = global_cell.size();
            cell_to_func_.resize(num_cells);
            for (int cell = 0; cell < num_cells; ++cell) {
                cell_to_func_[cell] = satnum[global_cell[cell]] - 1;
            }
        }

        // Find number of tables, check for consistency.
        enum { Uninitialized = -1 };
        int num_tables = Uninitialized;
        if (phase_usage_.phase_used[Aqua]) {
            const SWOF::table_t& swof_table = deck.getSWOF().swof_;
            num_tables = swof_table.size();
            if (num_tables < satfuncs_expected) {
                THROW("Found " << num_tables << " SWOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
        }
        if (phase_usage_.phase_used[Vapour]) {
            const SGOF::table_t& sgof_table = deck.getSGOF().sgof_;
            int num_sgof_tables = sgof_table.size();
            if (num_sgof_tables < satfuncs_expected) {
                THROW("Found " << num_tables << " SGOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
            if (num_tables == Uninitialized) {
                num_tables = num_sgof_tables;
            } else if (num_tables != num_sgof_tables) {
                THROW("Inconsistent number of tables in SWOF and SGOF.");
            }
        }

        // Initialize tables.
        satfuncset_.resize(num_tables);
        for (int table = 0; table < num_tables; ++table) {
            satfuncset_[table].init(deck, table, phase_usage_);
        }
    }




    /// \return   P, the number of phases.
    int SaturationPropsFromDeck::numPhases() const
    {
	return phase_usage_.num_phases;
    }




    /// Relative permeability.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dkr_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void SaturationPropsFromDeck::relperm(const int n,
                                          const double* s,
                                          const int* cells,
                                          double* kr,
                                          double* dkrds) const
    {
        const int np = phase_usage_.num_phases;
        if (dkrds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i);
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalKr(s + np*i, kr + np*i);
            }
        }
    }




    /// Capillary pressure.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
    /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dpc_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void SaturationPropsFromDeck::capPress(const int n,
                                           const double* s,
                                           const int* cells,
                                           double* pc,
                                           double* dpcds) const
    {
        const int np = phase_usage_.num_phases;
        if (dpcds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalPcDeriv(s + np*i, pc + np*i, dpcds + np*np*i);
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalPc(s + np*i, pc + np*i);
            }
        }
    }




    /// Obtain the range of allowable saturation values.
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    void SaturationPropsFromDeck::satRange(const int n,
                                           const int* cells,
					   double* smin,
					   double* smax) const
    {
	const int np = phase_usage_.num_phases;
	for (int i = 0; i < n; ++i) {
	    for (int p = 0; p < np; ++p) {
		smin[np*i + p] = funcForCell(cells[i]).smin_[p];
		smax[np*i + p] = funcForCell(cells[i]).smax_[p];
	    }
	}
    }


    // Map the cell number to the correct function set.
    const SaturationPropsFromDeck::SatFuncSet&
    SaturationPropsFromDeck::funcForCell(const int cell) const
    {
        return cell_to_func_.empty() ? satfuncset_[0] : satfuncset_[cell_to_func_[cell]];
    }



    // ----------- Methods of SatFuncSet below -----------


    void SaturationPropsFromDeck::SatFuncSet::init(const EclipseGridParser& deck,
                                                   const int table_num,
                                                   const PhaseUsage phase_usg)
    {
        phase_usage = phase_usg;
        const int samples = 200;
        double swco = 0.0;
        double swmax = 1.0;
        if (phase_usage.phase_used[Aqua]) {
            const SWOF::table_t& swof_table = deck.getSWOF().swof_;
            const std::vector<double>& sw = swof_table[table_num][0];
            const std::vector<double>& krw = swof_table[table_num][1];
            const std::vector<double>& krow = swof_table[table_num][2];
            const std::vector<double>& pcow = swof_table[table_num][3];
            buildUniformMonotoneTable(sw, krw,  samples, krw_);
            buildUniformMonotoneTable(sw, krow, samples, krow_);
            buildUniformMonotoneTable(sw, pcow, samples, pcow_);
            krocw_ = krow[0]; // At connate water -> ecl. SWOF
            swco = sw[0];
            smin_[phase_usage.phase_pos[Aqua]] = sw[0];
            swmax = sw.back();
            smax_[phase_usage.phase_pos[Aqua]] = sw.back();
        }
        if (phase_usage.phase_used[Vapour]) {
            const SGOF::table_t& sgof_table = deck.getSGOF().sgof_;
            const std::vector<double>& sg = sgof_table[table_num][0];
            const std::vector<double>& krg = sgof_table[table_num][1];
            const std::vector<double>& krog = sgof_table[table_num][2];
            const std::vector<double>& pcog = sgof_table[table_num][3];
            buildUniformMonotoneTable(sg, krg,  samples, krg_);
            buildUniformMonotoneTable(sg, krog, samples, krog_);
            buildUniformMonotoneTable(sg, pcog, samples, pcog_);
            smin_[phase_usage.phase_pos[Vapour]] = sg[0];
            if (std::fabs(sg.back() + swco - 1.0) > 1e-3) {
                THROW("Gas maximum saturation in SGOF table = " << sg.back() <<
                      ", should equal (1.0 - connate water sat) = " << (1.0 - swco));
            }
            smax_[phase_usage.phase_pos[Vapour]] = sg.back();
        }
        // These only consider water min/max sats. Consider gas sats?
        smin_[phase_usage.phase_pos[Liquid]] = 1.0 - swmax;
        smax_[phase_usage.phase_pos[Liquid]] = 1.0 - swco;
    }


    void SaturationPropsFromDeck::SatFuncSet::evalKr(const double* s, double* kr) const
    {
        if (phase_usage.num_phases == 3) {
            // Stone-II relative permeability model.
            double sw = s[Aqua];
            double sg = s[Vapour];
            double krw = krw_(sw);
            double krg = krg_(sg);
            double krow = krow_(sw + sg); // = 1 - so
            double krog = krog_(sg);      // = 1 - so - sw
            double krocw = krocw_;
            kr[Aqua] = krw;
            kr[Vapour] = krg;
            kr[Liquid] = krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[Liquid] < 0.0) {
                kr[Liquid] = 0.0;
            }
            return;
        }
        // We have a two-phase situation. We know that oil is active.
        if (phase_usage.phase_used[Aqua]) {
            int wpos = phase_usage.phase_pos[Aqua];
            int opos = phase_usage.phase_pos[Liquid];
            double sw = s[wpos];
            double krw = krw_(sw);
            double krow = krow_(sw);
            kr[wpos] = krw;
            kr[opos] = krow;
        } else {
            ASSERT(phase_usage.phase_used[Vapour]);
            int gpos = phase_usage.phase_pos[Vapour];
            int opos = phase_usage.phase_pos[Liquid];
            double sg = s[gpos];
            double krg = krg_(sg);
            double krog = krog_(sg);
            kr[gpos] = krg;
            kr[opos] = krog;
        }
    }


    void SaturationPropsFromDeck::SatFuncSet::evalKrDeriv(const double* s, double* kr, double* dkrds) const
    {
        const int np = phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            // Stone-II relative permeability model.
            double sw = s[Aqua];
            double sg = s[Vapour];
            double krw = krw_(sw);
            double dkrww = krw_.derivative(sw);
            double krg = krg_(sg);
            double dkrgg = krg_.derivative(sg);
            double krow = krow_(sw + sg);
            double dkrow = krow_.derivative(sw + sg);
            double krog = krog_(sg);
            double dkrog = krog_.derivative(sg);
            double krocw = krocw_;
            kr[Aqua] = krw;
            kr[Vapour] = krg;
            kr[Liquid] = krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[Liquid] < 0.0) {
                kr[Liquid] = 0.0;
            }
            dkrds[Aqua + Aqua*np] = dkrww;
            dkrds[Vapour + Vapour*np] = dkrgg;
            dkrds[Liquid + Aqua*np] = krocw*((dkrow/krocw + dkrww)*(krog/krocw + krg) - dkrww);
            dkrds[Liquid + Vapour*np] = krocw*((krow/krocw + krw)*(dkrog/krocw + dkrgg) - dkrgg)
                + krocw*((dkrow/krocw + krw)*(krog/krocw + krg) - dkrgg);
            return;
        }
        // We have a two-phase situation. We know that oil is active.
        if (phase_usage.phase_used[Aqua]) {
            int wpos = phase_usage.phase_pos[Aqua];
            int opos = phase_usage.phase_pos[Liquid];
            double sw = s[wpos];
            double krw = krw_(sw);
            double dkrww = krw_.derivative(sw);
            double krow = krow_(sw);
            double dkrow = krow_.derivative(sw);
            kr[wpos] = krw;
            kr[opos] = krow;
            dkrds[wpos + wpos*np] = dkrww;
            dkrds[opos + wpos*np] = dkrow; // Row opos, column wpos, fortran order.
        } else {
            ASSERT(phase_usage.phase_used[Vapour]);
            int gpos = phase_usage.phase_pos[Vapour];
            int opos = phase_usage.phase_pos[Liquid];
            double sg = s[gpos];
            double krg = krg_(sg);
            double dkrgg = krg_.derivative(sg);
            double krog = krog_(sg);
            double dkrog = krog_.derivative(sg);
            kr[gpos] = krg;
            kr[opos] = krog;
            dkrds[gpos + gpos*np] = dkrgg;
            dkrds[opos + gpos*np] = dkrog;
        }

    }


    void SaturationPropsFromDeck::SatFuncSet::evalPc(const double* s, double* pc) const
    {
        pc[phase_usage.phase_pos[Liquid]] = 0.0;
        if (phase_usage.phase_used[Aqua]) {
            int pos = phase_usage.phase_pos[Aqua];
            pc[pos] = pcow_(s[pos]);
        }
        if (phase_usage.phase_used[Vapour]) {
            int pos = phase_usage.phase_pos[Vapour];
            pc[pos] = pcog_(s[pos]);
        }
    }

    void SaturationPropsFromDeck::SatFuncSet::evalPcDeriv(const double* s, double* pc, double* dpcds) const
    {
	// The problem of determining three-phase capillary pressures
	// is very hard experimentally, usually one extends two-phase
	// data (as for relative permeability).
	// In our approach the derivative matrix is quite sparse, only
	// the diagonal elements corresponding to non-oil phases are
	// (potentially) nonzero.
        const int np = phase_usage.num_phases;
	std::fill(dpcds, dpcds + np*np, 0.0);
        pc[phase_usage.phase_pos[Liquid]] = 0.0;
        if (phase_usage.phase_used[Aqua]) {
            int pos = phase_usage.phase_pos[Aqua];
            pc[pos] = pcow_(s[pos]);
            dpcds[np*pos + pos] = pcow_.derivative(s[pos]);
        }
        if (phase_usage.phase_used[Vapour]) {
            int pos = phase_usage.phase_pos[Vapour];
            pc[pos] = pcog_(s[pos]);
            dpcds[np*pos + pos] = pcog_.derivative(s[pos]);
        }
    }




} // namespace Opm


