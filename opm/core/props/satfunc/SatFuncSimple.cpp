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
#include <opm/core/props/satfunc/SatFuncSimple.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/satfunc/SaturationPropsFromDeck.hpp>
#include <opm/core/grid.h>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/utility/buildUniformMonotoneTable.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{




    void SatFuncSimpleUniform::init(const EclipseGridParser& deck,
                             const int table_num,
                             const PhaseUsage phase_usg,
                             const int samples)
    {
        phase_usage = phase_usg;
        double swco = 0.0;
        double swmax = 1.0;
        if (phase_usage.phase_used[Aqua]) {
            const SWOF::table_t& swof_table = deck.getSWOF().swof_;
            const std::vector<double>& sw = swof_table[table_num][0];
            const std::vector<double>& krw = swof_table[table_num][1];
            const std::vector<double>& krow = swof_table[table_num][2];
            const std::vector<double>& pcow = swof_table[table_num][3];
            if (krw.front() != 0.0 || krow.back() != 0.0) {
                OPM_THROW(std::runtime_error, "Error SWOF data - non-zero krw(swco) and/or krow(1-sor)");
            }

            // Extend the tables with constant values such that the
            // derivatives at the endpoints are zero
            int n = sw.size();
            std::vector<double> sw_ex(n+2);
            std::vector<double> krw_ex(n+2);
            std::vector<double> krow_ex(n+2);
            std::vector<double> pcow_ex(n+2);

            SatFuncSimpleUniform::ExtendTable(sw,sw_ex,1);
            SatFuncSimpleUniform::ExtendTable(krw,krw_ex,0);
            SatFuncSimpleUniform::ExtendTable(krow,krow_ex,0);
            SatFuncSimpleUniform::ExtendTable(pcow,pcow_ex,0);

            buildUniformMonotoneTable(sw_ex, krw_ex,  samples, krw_);
            buildUniformMonotoneTable(sw_ex, krow_ex, samples, krow_);
            buildUniformMonotoneTable(sw_ex, pcow_ex, samples, pcow_);

            krocw_ = krow[0]; // At connate water -> ecl. SWOF
            swco = sw[0];
            smin_[phase_usage.phase_pos[Aqua]] = sw[0];
            swmax = sw.back();
            smax_[phase_usage.phase_pos[Aqua]] = sw.back();

            krwmax_ = krw.back();
            kromax_ = krow.front();
            swcr_ = swmax;
            sowcr_ = 1.0 - swco;
            krwr_ = krw.back();
            krorw_ = krow.front();
            for (std::vector<double>::size_type i=1; i<sw.size(); ++i) {
                if (krw[i]> 0.0) {
                   swcr_ = sw[i-1];
                   krorw_ = krow[i-1];
                   break;
                }
            }
            for (std::vector<double>::size_type i=sw.size()-1; i>=1; --i) {
                if (krow[i-1]> 0.0) {
                   sowcr_ = 1.0 - sw[i];
                   krwr_ = krw[i];
                   break;
                }
            }
        }
        if (phase_usage.phase_used[Vapour]) {
            const SGOF::table_t& sgof_table = deck.getSGOF().sgof_;
            const std::vector<double>& sg = sgof_table[table_num][0];
            const std::vector<double>& krg = sgof_table[table_num][1];
            const std::vector<double>& krog = sgof_table[table_num][2];
            const std::vector<double>& pcog = sgof_table[table_num][3];

            // Extend the tables with constant values such that the
            // derivatives at the endpoints are zero
            int n = sg.size();
            std::vector<double> sg_ex(n+2);
            std::vector<double> krg_ex(n+2);
            std::vector<double> krog_ex(n+2);
            std::vector<double> pcog_ex(n+2);

            SatFuncSimpleUniform::ExtendTable(sg,sg_ex,1);
            SatFuncSimpleUniform::ExtendTable(krg,krg_ex,0);
            SatFuncSimpleUniform::ExtendTable(krog,krog_ex,0);
            SatFuncSimpleUniform::ExtendTable(pcog,pcog_ex,0);

            buildUniformMonotoneTable(sg_ex, krg_ex,  samples, krg_);
            buildUniformMonotoneTable(sg_ex, krog_ex, samples, krog_);
            buildUniformMonotoneTable(sg_ex, pcog_ex, samples, pcog_);
            smin_[phase_usage.phase_pos[Vapour]] = sg[0];
            if (std::fabs(sg.back() + swco - 1.0) > 1e-3) {
                OPM_THROW(std::runtime_error, "Gas maximum saturation in SGOF table = " << sg.back() <<
                      ", should equal (1.0 - connate water sat) = " << (1.0 - swco));
            }
            smax_[phase_usage.phase_pos[Vapour]] = sg.back();
        }
        // These only consider water min/max sats. Consider gas sats?
        smin_[phase_usage.phase_pos[Liquid]] = 1.0 - swmax;
        smax_[phase_usage.phase_pos[Liquid]] = 1.0 - swco;
    }



    void SatFuncSimpleUniform::ExtendTable(const std::vector<double>& xv,
                                           std::vector<double>& xv_ex,
                                           double pm) const
    {
        int n = xv.size();
        xv_ex[0] = xv[0]-pm;
        xv_ex[n+1] = xv[n-1]+pm;
        for (int i=0; i<n; i++)
        {
            xv_ex[i+1] = xv[i];

        }

    }

    void SatFuncSimpleUniform::evalKr(const double* s, double* kr) const
    {
        if (phase_usage.num_phases == 3) {
            // A simplified relative permeability model.
            double sw = s[Aqua];
            double sg = s[Vapour];
            double krw = krw_(sw);
            double krg = krg_(sg);
            double krow = krow_(sw + sg); // = 1 - so
            // double krog = krog_(sg);      // = 1 - so - sw
            // double krocw = krocw_;
            kr[Aqua] = krw;
            kr[Vapour] = krg;
            kr[Liquid] = krow;
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
            double so = s[opos];
            double krow = krow_(1.0-so);
            kr[wpos] = krw;
            kr[opos] = krow;
        } else {
            assert(phase_usage.phase_used[Vapour]);
            int gpos = phase_usage.phase_pos[Vapour];
            int opos = phase_usage.phase_pos[Liquid];
            double sg = s[gpos];
            double krg = krg_(sg);
            double krog = krog_(sg);
            kr[gpos] = krg;
            kr[opos] = krog;
        }
    }


    void SatFuncSimpleUniform::evalKrDeriv(const double* s, double* kr, double* dkrds) const
    {
        const int np = phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            // A simplified relative permeability model.
            double sw = s[Aqua];
            double sg = s[Vapour];
            double krw = krw_(sw);
            double dkrww = krw_.derivative(sw);
            double krg = krg_(sg);
            double dkrgg = krg_.derivative(sg);
            double krow = krow_(sw + sg);
            double dkrow = krow_.derivative(sw + sg);
            // double krog = krog_(sg);
            // double dkrog = krog_.derivative(sg);
            // double krocw = krocw_;
            kr[Aqua] = krw;
            kr[Vapour] = krg;
            kr[Liquid] = krow;
            //krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[Liquid] < 0.0) {
                kr[Liquid] = 0.0;
            }
            dkrds[Aqua + Aqua*np] = dkrww;
            dkrds[Vapour + Vapour*np] = dkrgg;
            //dkrds[Liquid + Aqua*np] = dkrow;
            dkrds[Liquid + Liquid*np] = -dkrow;
                    //krocw*((dkrow/krocw + dkrww)*(krog/krocw + krg) - dkrww);
            dkrds[Liquid + Vapour*np] = 0.0;
                    //krocw*((krow/krocw + krw)*(dkrog/krocw + dkrgg) - dkrgg)
                    //+ krocw*((dkrow/krocw + krw)*(krog/krocw + krg) - dkrgg);
            return;
        }
        // We have a two-phase situation. We know that oil is active.
        if (phase_usage.phase_used[Aqua]) {
            int wpos = phase_usage.phase_pos[Aqua];
            int opos = phase_usage.phase_pos[Liquid];
            double sw = s[wpos];
            double krw = krw_(sw);
            double dkrww = krw_.derivative(sw);
            double so = s[opos];
            double krow = krow_(1.0-so);
            double dkrow = krow_.derivative(1.0-so);
            kr[wpos] = krw;
            kr[opos] = krow;
            dkrds[wpos + wpos*np] = dkrww;
            dkrds[opos + wpos*np] = dkrow; // Row opos, column wpos, fortran order.
        } else {
            assert(phase_usage.phase_used[Vapour]);
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


    void SatFuncSimpleUniform::evalPc(const double* s, double* pc) const
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

    void SatFuncSimpleUniform::evalPcDeriv(const double* s, double* pc, double* dpcds) const
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





    // ====== Methods for SatFuncSimpleNonuniform ======



    void SatFuncSimpleNonuniform::ExtendTable(const std::vector<double>& xv,
                                           std::vector<double>& xv_ex,
                                           double pm) const
    {
        int n = xv.size();
        xv_ex[0] = xv[0]-pm;
        xv_ex[n+1] = xv[n-1]+pm;
        for (int i=0; i<n; i++)
        {
            xv_ex[i+1] = xv[i];

        }

    }


    void SatFuncSimpleNonuniform::init(const EclipseGridParser& deck,
                                       const int table_num,
                                       const PhaseUsage phase_usg,
                                       const int /*samples*/)
    {
        phase_usage = phase_usg;
        double swco = 0.0;
        double swmax = 1.0;
        if (phase_usage.phase_used[Aqua]) {
            const SWOF::table_t& swof_table = deck.getSWOF().swof_;
            const std::vector<double>& sw = swof_table[table_num][0];
            const std::vector<double>& krw = swof_table[table_num][1];
            const std::vector<double>& krow = swof_table[table_num][2];
            const std::vector<double>& pcow = swof_table[table_num][3];
            if (krw.front() != 0.0 || krow.back() != 0.0) {
                OPM_THROW(std::runtime_error, "Error SWOF data - non-zero krw(swco) and/or krow(1-sor)");
            }

            // Extend the tables with constant values such that the
            // derivatives at the endpoints are zero
            int n = sw.size();
            std::vector<double> sw_ex(n+2);
            std::vector<double> krw_ex(n+2);
            std::vector<double> krow_ex(n+2);
            std::vector<double> pcow_ex(n+2);

            SatFuncSimpleNonuniform::ExtendTable(sw,sw_ex,1);
            SatFuncSimpleNonuniform::ExtendTable(krw,krw_ex,0);
            SatFuncSimpleNonuniform::ExtendTable(krow,krow_ex,0);
            SatFuncSimpleNonuniform::ExtendTable(pcow,pcow_ex,0);

            krw_ = NonuniformTableLinear<double>(sw_ex, krw_ex);
            krow_ = NonuniformTableLinear<double>(sw_ex, krow_ex);
            pcow_ = NonuniformTableLinear<double>(sw_ex, pcow_ex);
            krocw_ = krow[0]; // At connate water -> ecl. SWOF
            swco = sw[0];
            smin_[phase_usage.phase_pos[Aqua]] = sw[0];
            swmax = sw.back();
            smax_[phase_usage.phase_pos[Aqua]] = sw.back();

            krwmax_ = krw.back();
            kromax_ = krow.front();
            swcr_ = swmax;
            sowcr_ = 1.0 - swco;
            krwr_ = krw.back();
            krorw_ = krow.front();
            for (std::vector<double>::size_type i=1; i<sw.size(); ++i) {
                if (krw[i]> 0.0) {
                   swcr_ = sw[i-1];
                   krorw_ = krow[i-1];
                   break;
                }
            }
            for (std::vector<double>::size_type i=sw.size()-1; i>=1; --i) {
                if (krow[i-1]> 0.0) {
                   sowcr_ = 1.0 - sw[i];
                   krwr_ = krw[i];
                   break;
                }
            }

        }
        if (phase_usage.phase_used[Vapour]) {
            const SGOF::table_t& sgof_table = deck.getSGOF().sgof_;
            const std::vector<double>& sg = sgof_table[table_num][0];
            const std::vector<double>& krg = sgof_table[table_num][1];
            const std::vector<double>& krog = sgof_table[table_num][2];
            const std::vector<double>& pcog = sgof_table[table_num][3];

            // Extend the tables with constant values such that the
            // derivatives at the endpoints are zero
            int n = sg.size();
            std::vector<double> sg_ex(n+2);
            std::vector<double> krg_ex(n+2);
            std::vector<double> krog_ex(n+2);
            std::vector<double> pcog_ex(n+2);

            SatFuncSimpleNonuniform::ExtendTable(sg,sg_ex,1);
            SatFuncSimpleNonuniform::ExtendTable(krg,krg_ex,0);
            SatFuncSimpleNonuniform::ExtendTable(krog,krog_ex,0);
            SatFuncSimpleNonuniform::ExtendTable(pcog,pcog_ex,0);

            krg_ = NonuniformTableLinear<double>(sg_ex, krg_ex);
            krog_ = NonuniformTableLinear<double>(sg_ex, krog_ex);
            pcog_ = NonuniformTableLinear<double>(sg_ex, pcog_ex);
            smin_[phase_usage.phase_pos[Vapour]] = sg[0];
            if (std::fabs(sg.back() + swco - 1.0) > 1e-3) {
                OPM_THROW(std::runtime_error, "Gas maximum saturation in SGOF table = " << sg.back() <<
                      ", should equal (1.0 - connate water sat) = " << (1.0 - swco));
            }
            smax_[phase_usage.phase_pos[Vapour]] = sg.back();
        }
        // These only consider water min/max sats. Consider gas sats?
        smin_[phase_usage.phase_pos[Liquid]] = 1.0 - swmax;
        smax_[phase_usage.phase_pos[Liquid]] = 1.0 - swco;
    }


    void SatFuncSimpleNonuniform::evalKr(const double* s, double* kr) const
    {
        if (phase_usage.num_phases == 3) {
            // A simplified relative permeability model.
            double sw = s[Aqua];
            double sg = s[Vapour];
            double krw = krw_(sw);
            double krg = krg_(sg);
            double krow = krow_(sw + sg); // = 1 - so
            // double krog = krog_(sg);      // = 1 - so - sw
            // double krocw = krocw_;
            kr[Aqua] = krw;
            kr[Vapour] = krg;
            kr[Liquid] = krow;
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
            double so = s[opos];
            double krow = krow_(1.0-so);
            kr[wpos] = krw;
            kr[opos] = krow;
        } else {
            assert(phase_usage.phase_used[Vapour]);
            int gpos = phase_usage.phase_pos[Vapour];
            int opos = phase_usage.phase_pos[Liquid];
            double sg = s[gpos];
            double krg = krg_(sg);
            double krog = krog_(sg);
            kr[gpos] = krg;
            kr[opos] = krog;
        }
    }


    void SatFuncSimpleNonuniform::evalKrDeriv(const double* s, double* kr, double* dkrds) const
    {
        const int np = phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            // A simplified relative permeability model.
            double sw = s[Aqua];
            double sg = s[Vapour];
            double krw = krw_(sw);
            double dkrww = krw_.derivative(sw);
            double krg = krg_(sg);
            double dkrgg = krg_.derivative(sg);
            double krow = krow_(sw + sg);
            double dkrow = krow_.derivative(sw + sg);
            // double krog = krog_(sg);
            // double dkrog = krog_.derivative(sg);
            // double krocw = krocw_;
            kr[Aqua] = krw;
            kr[Vapour] = krg;
            kr[Liquid] = krow;
            //krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[Liquid] < 0.0) {
                kr[Liquid] = 0.0;
            }
            dkrds[Aqua + Aqua*np] = dkrww;
            dkrds[Vapour + Vapour*np] = dkrgg;
            //dkrds[Liquid + Aqua*np] = dkrow;
            dkrds[Liquid + Liquid*np] = -dkrow;
                    //krocw*((dkrow/krocw + dkrww)*(krog/krocw + krg) - dkrww);
            dkrds[Liquid + Vapour*np] = 0.0;
                    //krocw*((krow/krocw + krw)*(dkrog/krocw + dkrgg) - dkrgg)
                    //+ krocw*((dkrow/krocw + krw)*(krog/krocw + krg) - dkrgg);
            return;
        }
        // We have a two-phase situation. We know that oil is active.
        if (phase_usage.phase_used[Aqua]) {
            int wpos = phase_usage.phase_pos[Aqua];
            int opos = phase_usage.phase_pos[Liquid];
            double sw = s[wpos];
            double krw = krw_(sw);
            double dkrww = krw_.derivative(sw);
            double so = s[opos];
            double krow = krow_(1.0-so);
            double dkrow = krow_.derivative(1.0-so);
            kr[wpos] = krw;
            kr[opos] = krow;
            dkrds[wpos + wpos*np] = dkrww;
            dkrds[opos + wpos*np] = dkrow; // Row opos, column wpos, fortran order.
        } else {
            assert(phase_usage.phase_used[Vapour]);
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


    void SatFuncSimpleNonuniform::evalPc(const double* s, double* pc) const
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

    void SatFuncSimpleNonuniform::evalPcDeriv(const double* s, double* pc, double* dpcds) const
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
