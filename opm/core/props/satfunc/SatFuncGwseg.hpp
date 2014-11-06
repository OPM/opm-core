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
#ifndef SATFUNCGWSEG_HPP
#define SATFUNCGWSEG_HPP

#include <opm/core/props/satfunc/SatFuncBase.hpp>

namespace Opm
{
    template<class TableType>
    class SatFuncGwseg : public SatFuncBase<TableType>
    {
    public:
        void evalKr(const double* s, double* kr) const;
        void evalKrDeriv(const double* s, double* kr, double* dkrds) const;
        void evalPc(const double* s, double* pc) const;
        void evalPcDeriv(const double* s, double* pc, double* dpcds) const;

        void evalKr(const double* s, double* kr, const EPSTransforms* epst) const;
        void evalKr(const double* s, double* kr, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const;
        void evalKrDeriv(const double* s, double* kr, double* dkrds, const EPSTransforms* epst) const;
        void evalKrDeriv(const double* s, double* kr, double* dkrds, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const;
        void evalPc(const double* s, double* pc, const EPSTransforms* epst) const;
        void evalPcDeriv(const double* s, double* pc, double* dpcds, const EPSTransforms* epst) const;

    private:

    };
    
    typedef SatFuncGwseg<UniformTableLinear<double> > SatFuncGwsegUniform;
    typedef SatFuncGwseg<NonuniformTableLinear<double> > SatFuncGwsegNonuniform;

    template<class TableType>
    void SatFuncGwseg<TableType>::evalKr(const double* s, double* kr) const
    {
        if (this->phase_usage.num_phases == 3) {
            // Relative permeability model based on segregation of water
            // and gas, with oil present in both water and gas zones.
            double swco = this->smin_[this->phase_usage.phase_pos[BlackoilPhases::Aqua]];
            const double sw = std::max(s[BlackoilPhases::Aqua], swco);
            const double sg = s[BlackoilPhases::Vapour];
            const double eps = 1e-5;
            swco = std::min(swco,sw-eps);

            // xw and xg are the fractions occupied by water and gas zones.
            const double ssg = sw - swco + sg;
            const double xw = (sw - swco) / ssg;
            const double xg = 1 - xw;
            const double ssw = sg + sw;

            const double krw = this->krw_(sw);
            const double krg = this->krg_(sg);
            const double krow = this->krow_(ssw);
            const double krog = this->krog_(ssg);
            kr[BlackoilPhases::Aqua]   = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = xw*krow + xg*krog;
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
    void SatFuncGwseg<TableType>::evalKr(const double* s, double* kr, const EPSTransforms* epst) const
    {
        if (this->phase_usage.num_phases == 3) {     
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int gpos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];

            // Relative permeability model based on segregation of water
            // and gas, with oil present in both water and gas zones.

            // TODO  Also consider connate gas ...
            double _swco = this->smin_[this->phase_usage.phase_pos[BlackoilPhases::Aqua]];
            double swco = epst->wat.smin;
            const double sw = std::max(s[BlackoilPhases::Aqua], swco);
            const double sg = s[BlackoilPhases::Vapour];
            const double eps = 1e-6;
            swco = std::min(swco,sw-eps);
            const double ssw = sg + sw;
            const double ssg = std::max(sg + sw - swco, eps);
            const double d = ssg; // = sw - swco + sg (using 'd' for consistency with mrst docs).

            double ssow = 1.0-ssw;
            double ssog = 1.0-ssg-swco;
            double _sw = epst->wat.scaleSat(sw, 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double _sg = epst->gas.scaleSat(sg, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
            double _ssow = epst->watoil.scaleSat(ssow, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double _ssog = epst->gasoil.scaleSat(ssog, 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);

            const double krw = epst->wat.scaleKr(sw, this->krw_(_sw), this->krwr_);
            const double krg = epst->gas.scaleKr(sg, this->krg_(_sg), this->krgr_);
            const double krow = epst->watoil.scaleKr(ssow, this->krow_(1.0-_ssow), this->krorw_);
            const double krog = epst->gasoil.scaleKr(ssog, this->krog_(1.0-_ssog-_swco), this->krorg_);

            // xw and xg are the fractions occupied by water and gas zones.
            const double xw = (sw - swco) / d;
            const double xg = 1 - xw;
            kr[BlackoilPhases::Aqua]   = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = xw*krow + xg*krog;           
            return;
        }
        OPM_THROW(std::runtime_error, "SatFuncGwseg   --  need to be implemented ...");
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
    void SatFuncGwseg<TableType>::evalKr(const double* s, double* kr, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const
    {
        if (this->phase_usage.num_phases == 3) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int gpos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];

            // Relative permeability model based on segregation of water
            // and gas, with oil present in both water and gas zones.

            // TODO  Consider connate gas ...
            double _swco = this->smin_[this->phase_usage.phase_pos[BlackoilPhases::Aqua]];
            double swco = epst->wat.smin;
            const double sw = std::max(s[BlackoilPhases::Aqua], swco);
            const double sg = s[BlackoilPhases::Vapour];
            const double eps = 1e-6;
            swco = std::min(swco,sw-eps);
            const double ssw = sg + sw;
            const double ssg = std::max(sg + sw - swco, eps);
            const double d = ssg; // = sw - swco + sg (using 'd' for consistency with mrst docs).

            double ssow = 1.0-ssw;
            double ssog = 1.0-ssg-swco;
            
            // The code below corresponds to EHYSTR * 0 * * KR/ 
            //   - wettability properties water>oil>gas.
            //   - Carlsen hysteresis model for non-wetting (scanning=shifted_imb).  No hysteresis for wetting phase.
            // The imb-curve currently only differs from drainage curves via endpoint scaling ...
            
            // Water - use drainage curve only
            double _sw = epst->wat.scaleSat(sw, 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double krw = epst->wat.scaleKr(sw, this->krw_(_sw), this->krwr_);
            
            // Gas
            double krg;
            if (sg >= sat_hyst->sg_hyst) { // Drainage
                double _sg = epst->gas.scaleSat(sg, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
                krg = epst->gas.scaleKr(sg, this->krg_(_sg), this->krgr_);
            } else { // Imbibition
                double sg_shifted = sg + sat_hyst->sg_shift;
                double _sg =  epst_hyst->gas.scaleSat(sg_shifted, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
                krg = epst_hyst->gas.scaleKr(sg_shifted, this->krg_(_sg), this->krgr_);
            }       

            // Oil in water
            double krow;
            if (ssow >= sat_hyst->sow_hyst) { // Drainage
                double _ssow = epst->watoil.scaleSat(ssow, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
                krow = epst->watoil.scaleKr(ssow, this->krow_(1.0-_ssow), this->krorw_);
            } else { // Imbibition
                double ssow_shifted = ssow + sat_hyst->sow_shift;
                double _ssow = epst_hyst->watoil.scaleSat(ssow_shifted, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
                krow = epst_hyst->watoil.scaleKr(ssow_shifted, this->krow_(1.0-_ssow), this->krorw_);
            }
            
            // Oil in gas and connate water - use drainage curve only 
            double _ssog = epst->gasoil.scaleSat(ssog, 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double krog = epst->gasoil.scaleKr(ssog, this->krog_(1.0-_ssog-_swco), this->krorg_);

            // xw and xg are the fractions occupied by water and gas zones.
            const double xw = (sw - swco) / d;
            const double xg = 1 - xw;
            
            // relperms
            kr[BlackoilPhases::Aqua]   = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = xw*krow + xg*krog;

            return;       
        }
        OPM_THROW(std::runtime_error, "SatFuncGwseg   --  need to be implemented ...");
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
    void SatFuncGwseg<TableType>::evalKrDeriv(const double* s, double* kr, double* dkrds) const
    {
        const int np = this->phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            // Relative permeability model based on segregation of water
            // and gas, with oil present in both water and gas zones.

            double swco = this->smin_[this->phase_usage.phase_pos[BlackoilPhases::Aqua]];
            const double sw = std::max(s[BlackoilPhases::Aqua], swco);
            const double sg = s[BlackoilPhases::Vapour];
            const double eps = 1e-5;
            swco = std::min(swco,sw-eps);

            // xw and xg are the fractions occupied by water and gas zones.
            const double ssw = sg + sw;
            // d = ssg = sw - swco + sg (using 'd' for consistency with mrst docs).
            const double d = sg + sw - swco;
            const double xw = (sw - swco) / d;

            const double krw = this->krw_(sw);
            const double krg = this->krg_(sg);
            const double krow = this->krow_(ssw);
            const double krog = this->krog_(d);

            const double xg = 1 - xw;
            kr[BlackoilPhases::Aqua]   = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = xw*krow + xg*krog;

            // Derivatives.
            const double dkrww = this->krw_.derivative(sw);
            const double dkrgg = this->krg_.derivative(sg);
            const double dkrow = this->krow_.derivative(ssw);
            const double dkrog = this->krog_.derivative(d);
            dkrds[BlackoilPhases::Aqua   + BlackoilPhases::Aqua*np]   =  dkrww;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Aqua*np]   =  (xg/d)*krow + xw*dkrow - (xg/d)*krog + xg*dkrog;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Vapour*np] = -(xw/d)*krow + xw*dkrow + (xw/d)*krog + xg*dkrog;
            dkrds[BlackoilPhases::Vapour + BlackoilPhases::Vapour*np] =  dkrgg;
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
    void SatFuncGwseg<TableType>::evalKrDeriv(const double* s, double* kr, double* dkrds, const EPSTransforms* epst) const
    {
        const int np = this->phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int gpos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];

            // Relative permeability model based on segregation of water
            // and gas, with oil present in both water and gas zones.

            // TODO  Also consider connate gas ...
            double _swco = this->smin_[this->phase_usage.phase_pos[BlackoilPhases::Aqua]];
            double swco = epst->wat.smin;
            const double sw = std::max(s[BlackoilPhases::Aqua], swco);
            const double sg = s[BlackoilPhases::Vapour];
            const double eps = 1e-6;
            swco = std::min(swco,sw-eps);
            const double ssw = sg + sw;
            const double ssg = std::max(sg + sw - swco, eps);
            const double d = ssg; // = sw - swco + sg (using 'd' for consistency with mrst docs).

            double ssow = 1.0-ssw;
            double ssog = 1.0-ssg-swco;
            double _sw = epst->wat.scaleSat(sw, 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double _dsdsw = epst->wat.scaleSatDeriv(sw, 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double _sg = epst->gas.scaleSat(sg, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
            double _dsdsg = epst->gas.scaleSatDeriv(sg, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
            double _ssow = epst->watoil.scaleSat(ssow, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double _dsdssow = epst->watoil.scaleSatDeriv(ssow, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double _ssog = epst->gasoil.scaleSat(ssog, 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double _dsdssog = epst->gasoil.scaleSatDeriv(ssog, 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);

            const double krw = epst->wat.scaleKr(sw, this->krw_(_sw), this->krwr_);
            const double krg = epst->gas.scaleKr(sg, this->krg_(_sg), this->krgr_);
            const double krow = epst->watoil.scaleKr(ssow, this->krow_(std::max(1.0-_ssow,_swco)), this->krorw_);
            const double krog = epst->gasoil.scaleKr(ssog, this->krog_(std::max(1.0-_ssog-_swco,eps)), this->krorg_);

            // xw and xg are the fractions occupied by water and gas zones.
            const double xw = (sw - swco) / d;
            const double xg = 1 - xw;
            kr[BlackoilPhases::Aqua]   = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = xw*krow + xg*krog;

            // Derivatives.
            double dkrww = _dsdsw*epst->wat.scaleKrDeriv(sw, this->krw_.derivative(_sw));
            double dkrgg = _dsdsg*epst->gas.scaleKrDeriv(sg, this->krg_.derivative(_sg));
            double dkrow = _dsdssow*epst->watoil.scaleKrDeriv(ssow, this->krow_.derivative(std::max(1.0-_ssow,_swco)));
            double dkrog = _dsdssog*epst->gasoil.scaleKrDeriv(ssog, this->krog_.derivative(std::max(1.0-_ssog-_swco,eps)));           
            dkrds[BlackoilPhases::Aqua   + BlackoilPhases::Aqua*np]   =  dkrww;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Aqua*np]   =  (xg/d)*krow + xw*dkrow - (xg/d)*krog + xg*dkrog;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Vapour*np] = -(xw/d)*krow + xw*dkrow + (xw/d)*krog + xg*dkrog;
            dkrds[BlackoilPhases::Vapour + BlackoilPhases::Vapour*np] =  dkrgg;
            return;
        }
        OPM_THROW(std::runtime_error, "SatFuncGwseg   --  need to be implemented ...");
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
    void SatFuncGwseg<TableType>::evalKrDeriv(const double* s, double* kr, double* dkrds, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const
    {
        const int np = this->phase_usage.num_phases;
        std::fill(dkrds, dkrds + np*np, 0.0);

        if (np == 3) {
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            int gpos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];

            // Relative permeability model based on segregation of water
            // and gas, with oil present in both water and gas zones.

            // TODO  Also consider connate gas ...
            double _swco = this->smin_[this->phase_usage.phase_pos[BlackoilPhases::Aqua]];
            double swco = epst->wat.smin;
            const double sw = std::max(s[BlackoilPhases::Aqua], swco);
            const double sg = s[BlackoilPhases::Vapour];
            const double eps = 1e-6;
            swco = std::min(swco,sw-eps);
            const double ssw = sg + sw;
            const double ssg = std::max(sg + sw - swco, eps);
            const double d = ssg; // = sw - swco + sg (using 'd' for consistency with mrst docs).

            double ssow = 1.0-ssw;
            double ssog = 1.0-ssg-swco;
            
            // The code below corresponds to EHYSTR * 0 * * KR/ 
            //   - wettability properties water>oil>gas.
            //   - Carlsen hysteresis model for non-wetting (scanning=shifted_imb).  No hysteresis for wetting phase.
            // The imb-curve currently only differs from drainage curves via endpoint scaling ...
            
            // Water - use drainage curve only
            double _sw = epst->wat.scaleSat(sw, 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double _dsdsw = epst->wat.scaleSatDeriv(sw, 1.0-this->sowcr_-this->smin_[gpos], this->swcr_, this->smax_[wpos]);
            double krw = epst->wat.scaleKr(sw, this->krw_(_sw), this->krwr_);
            double dkrww = _dsdsw*epst->wat.scaleKrDeriv(sw, this->krw_.derivative(_sw));
            
            // Gas
            double krg, dkrgg;
            if (sg >= sat_hyst->sg_hyst) { // Drainage
                double _sg = epst->gas.scaleSat(sg, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
                double _dsdsg = epst->gas.scaleSatDeriv(sg, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
                krg = epst->gas.scaleKr(sg, this->krg_(_sg), this->krgr_);
                dkrgg = _dsdsg*epst->gas.scaleKrDeriv(sg, this->krg_.derivative(_sg));
            } else { // Imbibition
                double sg_shifted = sg + sat_hyst->sg_shift;
                double _sg =  epst_hyst->gas.scaleSat(sg_shifted, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
                double _dsdsg = epst_hyst->gas.scaleSatDeriv(sg_shifted, 1.0-this->sogcr_-this->smin_[wpos], this->sgcr_, this->smax_[gpos]);
                krg = epst_hyst->gas.scaleKr(sg_shifted, this->krg_(_sg), this->krgr_);
                dkrgg = _dsdsg*epst_hyst->gas.scaleKrDeriv(sg_shifted, this->krg_.derivative(_sg));
            }       

            // Oil in water
            double krow, dkrow;
            if (ssow >= sat_hyst->sow_hyst) { // Drainage
                double _ssow = epst->watoil.scaleSat(ssow, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
                double _dsdssow = epst->watoil.scaleSatDeriv(ssow, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
                krow = epst->watoil.scaleKr(ssow, this->krow_(std::max(1.0-_ssow,_swco)), this->krorw_);
                dkrow = _dsdssow*epst->watoil.scaleKrDeriv(ssow, this->krow_.derivative(std::max(1.0-_ssow,_swco)));
            } else { // Imbibition
                double ssow_shifted = ssow + sat_hyst->sow_shift;
                double _ssow = epst_hyst->watoil.scaleSat(ssow_shifted, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
                double _dsdssow = epst_hyst->watoil.scaleSatDeriv(ssow_shifted, 1.0-this->swcr_-this->smin_[gpos], this->sowcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
                krow = epst_hyst->watoil.scaleKr(ssow_shifted, this->krow_(std::max(1.0-_ssow,_swco)), this->krorw_);
                dkrow = _dsdssow*epst_hyst->watoil.scaleKrDeriv(ssow_shifted, this->krow_.derivative(std::max(1.0-_ssow,_swco)));
            }
            
            // Oil in gas and connate water - use drainage curve only 
            double _ssog = epst->gasoil.scaleSat(ssog, 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double _dsdssog = epst->gasoil.scaleSatDeriv(ssog, 1.0-this->sgcr_-this->smin_[wpos], this->sogcr_, 1.0-this->smin_[wpos]-this->smin_[gpos]);
            double krog = epst->gasoil.scaleKr(ssog, this->krog_(std::max(1.0-_ssog-_swco,eps)), this->krorg_);
            double dkrog = _dsdssog*epst->gasoil.scaleKrDeriv(ssog, this->krog_.derivative(std::max(1.0-_ssog-_swco,eps)));

            // xw and xg are the fractions occupied by water and gas zones.
            const double xw = (sw - swco) / d;
            const double xg = 1 - xw;
            
            // relperms
            kr[BlackoilPhases::Aqua]   = krw;
            kr[BlackoilPhases::Vapour] = krg;
            kr[BlackoilPhases::Liquid] = xw*krow + xg*krog;

            // Derivatives.      
            dkrds[BlackoilPhases::Aqua   + BlackoilPhases::Aqua*np]   =  dkrww;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Aqua*np]   =  (xg/d)*krow + xw*dkrow - (xg/d)*krog + xg*dkrog;
            dkrds[BlackoilPhases::Liquid + BlackoilPhases::Vapour*np] = -(xw/d)*krow + xw*dkrow + (xw/d)*krog + xg*dkrog;
            dkrds[BlackoilPhases::Vapour + BlackoilPhases::Vapour*np] =  dkrgg;
            return;
        }
        OPM_THROW(std::runtime_error, "SatFuncGwseg   --  need to be implemented ...");
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
    void SatFuncGwseg<TableType>::evalPc(const double* s, double* pc) const
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
    void SatFuncGwseg<TableType>::evalPc(const double* s, double* pc, const EPSTransforms* epst) const
    {
        pc[this->phase_usage.phase_pos[BlackoilPhases::Liquid]] = 0.0;
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            int pos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            double _sw = epst->wat.scaleSatPc(s[pos], this->smin_[pos], this->smax_[pos]);
            pc[pos] = epst->wat.pcFactor*this->pcow_(_sw);
        }
        if (this->phase_usage.phase_used[BlackoilPhases::Vapour]) {
            int pos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];
            double _sg = epst->gas.scaleSatPc(s[pos], this->smin_[pos], this->smax_[pos]);
            pc[pos] = epst->gas.pcFactor*this->pcog_(_sg);
        }
    }

    template<class TableType>
    void SatFuncGwseg<TableType>::evalPcDeriv(const double* s, double* pc, double* dpcds) const
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


    template<class TableType>
    void SatFuncGwseg<TableType>::evalPcDeriv(const double* s, double* pc, double* dpcds, const EPSTransforms* epst) const
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
            double _sw = epst->wat.scaleSatPc(s[pos], this->smin_[pos], this->smax_[pos]);
            pc[pos] = epst->wat.pcFactor*this->pcow_(s[pos]);
            double _dsdsw = epst->wat.scaleSatDerivPc(s[pos], this->smin_[pos], this->smax_[pos]);
            dpcds[np*pos + pos] = epst->wat.pcFactor*_dsdsw*this->pcow_.derivative(_sw);
        }
        if (this->phase_usage.phase_used[BlackoilPhases::Vapour]) {
            int pos = this->phase_usage.phase_pos[BlackoilPhases::Vapour];
            double _sg = epst->gas.scaleSatPc(s[pos], this->smin_[pos], this->smax_[pos]);
            pc[pos] = epst->gas.pcFactor*this->pcog_(_sg);
            double _dsdsg = epst->gas.scaleSatDerivPc(s[pos], this->smin_[pos], this->smax_[pos]);
            dpcds[np*pos + pos] = epst->gas.pcFactor*_dsdsg*this->pcog_.derivative(_sg);
        }
    }


} // namespace Opm
#endif // SATFUNCGWSEG_HPP
