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
#ifndef SATFUNCSIMPLE_HPP
#define SATFUNCSIMPLE_HPP

#include <opm/core/props/satfunc/SatFuncBase.hpp>

namespace Opm
{
    template<class TableType>
    class SatFuncSimple : public SatFuncBase<TableType>
    {
    public:
        void evalKr(const double* s, double* kr) const;
        void evalKrDeriv(const double* s, double* kr, double* dkrds) const;
        void evalPc(const double* s, double* pc) const;
        void evalPcDeriv(const double* s, double* pc, double* dpcds) const;

        void evalKr(const double* /* s */, double* /* kr */, const EPSTransforms* /* epst */) const
        {OPM_THROW(std::runtime_error, "SatFuncSimple   --  need to be implemented ...");}
        void evalKr(const double* /* s */, double* /* kr */, const EPSTransforms* /* epst */, const EPSTransforms* /* epst_hyst */, const SatHyst* /* sat_hyst */) const
        {OPM_THROW(std::runtime_error, "SatFuncSimple   --  need to be implemented ...");}
        void evalKrDeriv(const double* /* s */, double* /* kr */, double* /* dkrds */, const EPSTransforms* /* epst */) const;
        void evalKrDeriv(const double* /* s */, double* /* kr */, double* /* dkrds */, const EPSTransforms* /* epst */, const EPSTransforms* /* epst_hyst */, const SatHyst* /* sat_hyst */) const 
        {OPM_THROW(std::runtime_error, "SatFuncSimple   --  need to be implemented ...");}
        void evalPc(const double* /* s */, double* /* pc */, const EPSTransforms* /* epst */) const
        {OPM_THROW(std::runtime_error, "SatFuncSimple   --  need to be implemented ...");}
        void evalPcDeriv(const double* /* s */, double* /* pc */, double* /* dpcds */, const EPSTransforms* /* epst */) const
        {OPM_THROW(std::runtime_error, "SatFuncSimple   --  need to be implemented ...");}

    private:

    };
    
    typedef SatFuncSimple<UniformTableLinear<double> > SatFuncSimpleUniform;
    typedef SatFuncSimple<NonuniformTableLinear<double> > SatFuncSimpleNonuniform;

    template<class TableType>
    void SatFuncSimple<TableType>::evalKr(const double* s, double* kr) const
    {
        if (this->phase_usage.num_phases == 3) {
            // A simplified relative permeability model.
            double sw = s[BlackoilPhases::Aqua];
            double sg = s[BlackoilPhases::Vapour];
            double krw = this->krw_(sw);
            double krg = this->krg_(sg);
            double krow = this->krow_(sw + sg); // = 1 - so
            // double krog = krog_(sg);      // = 1 - so - sw
            // double krocw = krocw_;
            kr[BlackoilPhases::Aqua] = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = krow;
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
            double so = s[opos];
            double krow = this->krow_(1.0-so);
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
    void SatFuncSimple<TableType>::evalKrDeriv(const double* s, double* kr, double* dkrds) const
    {
        const int np = this->phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            // A simplified relative permeability model.
            double sw = s[BlackoilPhases::Aqua];
            double sg = s[BlackoilPhases::Vapour];
            double krw = this->krw_(sw);
            double dkrww = this->krw_.derivative(sw);
            double krg = this->krg_(sg);
            double dkrgg = this->krg_.derivative(sg);
            double krow = this->krow_(sw + sg);
            double dkrow = this->krow_.derivative(sw + sg);
            // double krog = krog_(sg);
            // double dkrog = krog_.derivative(sg);
            // double krocw = krocw_;
            kr[BlackoilPhases::Aqua] = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = krow;
            //krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[BlackoilPhases::Liquid] < 0.0) {
                kr[BlackoilPhases::Liquid] = 0.0;
            }
            dkrds[BlackoilPhases::Aqua + BlackoilPhases::Aqua*np] = dkrww;
            dkrds[BlackoilPhases::Vapour + BlackoilPhases::Vapour*np] = dkrgg;
            //dkrds[Liquid + Aqua*np] = dkrow;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Liquid*np] = -dkrow;
                    //krocw*((dkrow/krocw + dkrww)*(krog/krocw + krg) - dkrww);
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Vapour*np] = 0.0;
                    //krocw*((krow/krocw + krw)*(dkrog/krocw + dkrgg) - dkrgg)
                    //+ krocw*((dkrow/krocw + krw)*(krog/krocw + krg) - dkrgg);
            return;
        }
        // We have a two-phase situation. We know that oil is active.
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int opos = this->phase_usage.phase_pos[BlackoilPhases::Liquid];
            double sw = s[wpos];
            double krw = this->krw_(sw);
            double dkrww = this->krw_.derivative(sw);
            double so = s[opos];
            double krow = this->krow_(1.0-so);
            double dkrow = this->krow_.derivative(1.0-so);
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
    void SatFuncSimple<TableType>::evalKrDeriv(const double* s, double* kr, double* dkrds, const EPSTransforms* epst) const
    {
        const int np = this->phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int gpos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];
            // A simplified relative permeability model.
            // Define KR(s) = scaleKr(kr(scalSat(s)))
            // Thus KR'(s) = scaleKr'(kr(scaleSat(s)))*kr'((scaleSat(s))*scaleSat'(s)
            double _sw = epst->wat.scaleSat(s[wpos], 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double _dsdsw = epst->wat.scaleSatDeriv(s[wpos], 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double _sg = epst->gas.scaleSat(s[gpos], 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
            double _dsdsg = epst->gas.scaleSatDeriv(s[gpos], 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
            double _sow = epst->watoil.scaleSat(1.0-s[wpos]-s[gpos], 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double _dsdsow = epst->watoil.scaleSatDeriv(1.0-s[wpos]-s[gpos], 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            //double _sog = epst->gasoil.scaleSat(1.0-s[wpos]-s[gpos], 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            //double _dsdsog = epst->gasoil.scaleSatDeriv(1.0-s[wpos]-s[gpos], 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);

            double krw = epst->wat.scaleKr(s[wpos], this->krw_(_sw), this->krwr_);
            double dkrww = _dsdsw*epst->wat.scaleKrDeriv(s[wpos], this->krw_.derivative(_sw));
            double krg = epst->gas.scaleKr(s[gpos], this->krg_(_sg), this->krgr_);
            double dkrgg = _dsdsg*epst->gas.scaleKrDeriv(s[gpos], this->krg_.derivative(_sg));
            // TODO Check the arguments to the krow- and krog-tables below...  
            double krow = epst->watoil.scaleKr(1.0-s[wpos]-s[gpos], this->krow_(1.0-_sow-this->smin_[gpos]), this->krorw_);   // ????
            double dkrow = _dsdsow*epst->watoil.scaleKrDeriv(1.0-s[wpos]-s[gpos], this->krow_.derivative(1.0-_sow-this->smin_[gpos])); // ????
            //double krog = epst->gasoil.scaleKr(this->krog_(1.0-_sog-this->smin_[wpos]), 1.0-s[wpos]-s[gpos], this->krorg_);   // ????
            //double dkrog = _dsdsog*epst->gasoil.scaleKrDeriv(1.0-s[wpos]-s[gpos], this->krog_.derivative(1.0-_sog-this->smin_[wpos]));   // ????
            // double krocw = krocw_;
            kr[wpos] = krw;
            kr[gpos] = krg;
            kr[BlackoilPhases::Liquid] = krow;
            //krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
            if (kr[BlackoilPhases::Liquid] < 0.0) {
                kr[BlackoilPhases::Liquid] = 0.0;
            }
            dkrds[wpos + wpos*np] = dkrww;
            dkrds[gpos + gpos*np] = dkrgg;
            //dkrds[Liquid + Aqua*np] = dkrow;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Liquid*np] = -dkrow;
                    //krocw*((dkrow/krocw + dkrww)*(krog/krocw + krg) - dkrww);
            dkrds[BlackoilPhases::Liquid + gpos*np] = 0.0;
                    //krocw*((krow/krocw + krw)*(dkrog/krocw + dkrgg) - dkrgg)
                    //+ krocw*((dkrow/krocw + krw)*(krog/krocw + krg) - dkrgg);
            return;
        }
        OPM_THROW(std::runtime_error, "SatFuncSimple   --  need to be implemented ...");
        // We have a two-phase situation. We know that oil is active.
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int opos = this->phase_usage.phase_pos[BlackoilPhases::Liquid];
            double sw = s[wpos];
            double krw = this->krw_(sw);
            double dkrww = this->krw_.derivative(sw);
            double so = s[opos];
            double krow = this->krow_(1.0-so);
            double dkrow = this->krow_.derivative(1.0-so);
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
    void SatFuncSimple<TableType>::evalPc(const double* s, double* pc) const
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
    void SatFuncSimple<TableType>::evalPcDeriv(const double* s, double* pc, double* dpcds) const
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
#endif // SATFUNCSIMPLE_HPP
