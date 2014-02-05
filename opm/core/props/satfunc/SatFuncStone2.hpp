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
#ifndef SATFUNCSTONE2_HPP
#define SATFUNCSTONE2_HPP

#include <opm/core/props/satfunc/SatFuncBase.hpp>

namespace Opm
{
    template<class TableType>
    class SatFuncStone2 : public SatFuncBase<TableType>
    {
    public:
        void evalKr(const double* s, double* kr) const;
        void evalKrDeriv(const double* s, double* kr, double* dkrds) const;
        void evalPc(const double* s, double* pc) const;
        void evalPcDeriv(const double* s, double* pc, double* dpcds) const;

        void evalKr(const double* /*s*/, double* /*kr*/, const EPSTransforms* /*epst*/) const
        {OPM_THROW(std::runtime_error, "SatFuncStone2   --  need to be implemented ...");}
        void evalKr(const double* /*s*/, double* /*kr*/, const EPSTransforms* /*epst*/, const EPSTransforms* /*epst_hyst*/, const SatHyst* /*sat_hyst*/) const
        {OPM_THROW(std::runtime_error, "SatFuncStone2   --  need to be implemented ...");}
        void evalKrDeriv(const double* /*s*/, double* /*kr*/, double* /*dkrds*/, const EPSTransforms* /*epst*/) const
        {OPM_THROW(std::runtime_error, "SatFuncStone2   --  need to be implemented ...");}
        void evalKrDeriv(const double* /*s*/, double* /*kr*/, double* /*dkrds*/, const EPSTransforms* /*epst*/, const EPSTransforms* /*epst_hyst*/, const SatHyst* /*sat_hyst*/) const
        {OPM_THROW(std::runtime_error, "SatFuncStone2   --  need to be implemented ...");}
        void evalPc(const double* /*s*/, double* /*pc*/, const EPSTransforms* /*epst*/) const
        {OPM_THROW(std::runtime_error, "SatFuncStone2   --  need to be implemented ...");}
        void evalPcDeriv(const double* /*s*/, double* /*pc*/, double* /*dpcds*/, const EPSTransforms* /*epst*/) const
        {OPM_THROW(std::runtime_error, "SatFuncStone2   --  need to be implemented ...");}


    private:

    };

    typedef SatFuncStone2<UniformTableLinear<double> > SatFuncStone2Uniform;
    typedef SatFuncStone2<NonuniformTableLinear<double> > SatFuncStone2Nonuniform;

    template<class TableType>
    void SatFuncStone2<TableType>::evalKr(const double* s, double* kr) const
    {
        if (this->phase_usage.num_phases == 3) {
            // Stone-II relative permeability model.
            double sw = s[BlackoilPhases::Aqua];
            double sg = s[BlackoilPhases::Vapour];
            double krw = this->krw_(sw);
            double krg = this->krg_(sg);
            double krow = this->krow_(sw + sg); // = 1 - so
            double krog = this->krog_(sg);      // = 1 - so - sw
            double krocw = this->krocw_;
            kr[BlackoilPhases::Aqua] = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[BlackoilPhases::Liquid] < 0.0) {
                kr[BlackoilPhases::Liquid] = 0.0;
            }
            return;
        }
        // We have a two-phase situation. We know that oil is active.
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int opos = this->phase_usage.phase_pos[BlackoilPhases::Liquid];
            double sw = s[wpos];
            double krw = this->krw_(sw);
            double krow = this->krow_(sw);
            kr[wpos] = krw;
            kr[opos] = krow;
        } else {
            assert(this->phase_usage.phase_used[BlackoilPhases::Vapour]);
            int gpos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];
            int opos = this->phase_usage.phase_pos[BlackoilPhases::Liquid];
            double sg = s[gpos];
            double krg = this->krg_(sg);
            double krog = this->krog_(sg);
            kr[gpos] = krg;
            kr[opos] = krog;
        }
    }

    template<class TableType>
    void SatFuncStone2<TableType>::evalKrDeriv(const double* s, double* kr, double* dkrds) const
    {
        const int np = this->phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            // Stone-II relative permeability model.
            double sw = s[BlackoilPhases::Aqua];
            double sg = s[BlackoilPhases::Vapour];
            double krw = this->krw_(sw);
            double dkrww = this->krw_.derivative(sw);
            double krg = this->krg_(sg);
            double dkrgg = this->krg_.derivative(sg);
            double krow = this->krow_(sw + sg);
            double dkrow = this->krow_.derivative(sw + sg);
            double krog = this->krog_(sg);
            double dkrog = this->krog_.derivative(sg);
            double krocw = this->krocw_;
            kr[BlackoilPhases::Aqua] = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[BlackoilPhases::Liquid] < 0.0) {
                kr[BlackoilPhases::Liquid] = 0.0;
            }
            dkrds[BlackoilPhases::Aqua + BlackoilPhases::Aqua*np] = dkrww;
            dkrds[BlackoilPhases::Vapour + BlackoilPhases::Vapour*np] = dkrgg;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Aqua*np] = krocw*((dkrow/krocw + dkrww)*(krog/krocw + krg) - dkrww);
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Vapour*np] = krocw*((krow/krocw + krw)*(dkrog/krocw + dkrgg) - dkrgg)
                    + krocw*((dkrow/krocw + krw)*(krog/krocw + krg) - dkrgg);
            return;
        }
        // We have a two-phase situation. We know that oil is active.
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int opos = this->phase_usage.phase_pos[BlackoilPhases::Liquid];
            double sw = s[wpos];
            double krw = this->krw_(sw);
            double dkrww = this->krw_.derivative(sw);
            double krow = this->krow_(sw);
            double dkrow = this->krow_.derivative(sw);
            kr[wpos] = krw;
            kr[opos] = krow;
            dkrds[wpos + wpos*np] = dkrww;
            dkrds[opos + wpos*np] = dkrow; // Row opos, column wpos, fortran order.
        } else {
            assert(this->phase_usage.phase_used[BlackoilPhases::Vapour]);
            int gpos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];
            int opos = this->phase_usage.phase_pos[BlackoilPhases::Liquid];
            double sg = s[gpos];
            double krg = this->krg_(sg);
            double dkrgg = this->krg_.derivative(sg);
            double krog = this->krog_(sg);
            double dkrog = this->krog_.derivative(sg);
            kr[gpos] = krg;
            kr[opos] = krog;
            dkrds[gpos + gpos*np] = dkrgg;
            dkrds[opos + gpos*np] = dkrog;
        }

    }

    template<class TableType>
    void SatFuncStone2<TableType>::evalPc(const double* s, double* pc) const
    {
        pc[this->phase_usage.phase_pos[BlackoilPhases::Liquid]] = 0.0;
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            int pos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            pc[pos] = this->pcow_(s[pos]);
        }
        if (this->phase_usage.phase_used[BlackoilPhases::Vapour]) {
            int pos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];
            pc[pos] = this->pcog_(s[pos]);
        }
    }

    template<class TableType>
    void SatFuncStone2<TableType>::evalPcDeriv(const double* s, double* pc, double* dpcds) const
    {
        // The problem of determining three-phase capillary pressures
        // is very hard experimentally, usually one extends two-phase
        // data (as for relative permeability).
        // In our approach the derivative matrix is quite sparse, only
        // the diagonal elements corresponding to non-oil phases are
        // (potentially) nonzero.
        const int np = this->phase_usage.num_phases;
        std::fill(dpcds, dpcds + np*np, 0.0);
        pc[this->phase_usage.phase_pos[BlackoilPhases::Liquid]] = 0.0;
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            int pos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            pc[pos] = this->pcow_(s[pos]);
            dpcds[np*pos + pos] = this->pcow_.derivative(s[pos]);
        }
        if (this->phase_usage.phase_used[BlackoilPhases::Vapour]) {
            int pos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];
            pc[pos] = this->pcog_(s[pos]);
            dpcds[np*pos + pos] = this->pcog_.derivative(s[pos]);
        }
    }

} // namespace Opm
#endif // SATFUNCSTONE2_HPP
