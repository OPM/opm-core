/*
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
#ifndef SATFUNCBASE_HPP
#define SATFUNCBASE_HPP

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/utility/buildUniformMonotoneTable.hpp>
#include <opm/core/utility/NonuniformTableLinear.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <vector>

namespace Opm
{
    // Transforms for saturation table scaling
    struct EPSTransforms {
        struct Transform {
            bool doNotScale;
            bool do_3pt;
            double smin;
            double scr;
            double sr;
            double smax;
            double slope1;
            double slope2;
            double scaleSat(double ss, double s_r, double s_cr, double s_max) const;
            double scaleSatInv(double s, double s_r, double s_cr, double s_max) const;
            double scaleSatDeriv(double s, double s_r, double s_cr, double s_max) const; // Returns scaleSat'(s)
            double scaleSatPc(double s, double s_min, double s_max) const;
            double scaleSatDerivPc(double s, double s_min, double s_max) const; // Returns scaleSatPc'(s)
            bool doKrMax;
            bool doKrCrit;
            bool doSatInterp;
            double krsr;
            double krmax;
            double krSlopeMax;
            double krSlopeCrit;
            double scaleKr(double s, double kr, double krsr_tab) const;
            double scaleKrDeriv(double s, double krDeriv) const;   // Returns scaleKr'(kr(scaleSat(s)))*kr'((scaleSat(s))
            void printMe(std::ostream & out);
        };

        Transform wat;
        Transform watoil;
        Transform gas;
        Transform gasoil;
    };
    
    // Hysteresis
    struct SatHyst {
        double sg_hyst;
        double sg_shift;
        double sow_hyst;
        double sow_shift;
        void printMe(std::ostream & out);
    };


    template <class TableType>
    class SatFuncBase : public BlackoilPhases
    {
    public:
        void init(const EclipseGridParser& deck,
                  const int table_num,
                  const PhaseUsage phase_usg,
                  const int samples);
        void updateSatHyst(const double* s,
                           const EPSTransforms* epst, 
                           const EPSTransforms* epst_hyst, 
                           SatHyst* sat_hyst) const;
        double smin_[PhaseUsage::MaxNumPhases];
        double smax_[PhaseUsage::MaxNumPhases];
        double krwmax_; // Max water relperm
        double krgmax_; // Max gas relperm
        double kromax_; // Max oil relperm
        double swcr_;   // Critical water saturation.
        double sgcr_;   // Critical gas saturation.
        double krwr_;   // Water relperm at critical oil-in-water saturation.
        double krgr_;   // Gas relperm at critical oil-in-gas saturation.
        double sowcr_;  // Critical oil-in-water saturation.
        double sogcr_;  // Critical oil-in-gas-and-connate-water saturation.
        double krorw_;  // Oil relperm at critical water saturation.
        double krorg_;  // Oil relperm at critical gas saturation.

    protected:
        PhaseUsage phase_usage; // A copy of the outer class' phase_usage_.
        TableType krw_;
        TableType krow_;
        TableType pcow_;
        TableType krg_;
        TableType krog_;
        TableType pcog_;
        double krocw_; // = krow_(s_wc)

    private:
        void extendTable(const std::vector<double>& xv,
                         std::vector<double>& xv_ex,
                         double pm) const;
        void initializeTableType(TableType& table,
                                 const std::vector<double>& arg,
                                 const std::vector<double>& value,
                                 const int samples);
    };

    template <class TableType>
    void SatFuncBase<TableType>::init(const EclipseGridParser& deck,
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

            extendTable(sw,sw_ex,1);
            extendTable(krw,krw_ex,0);
            extendTable(krow,krow_ex,0);
            extendTable(pcow,pcow_ex,0);

            initializeTableType(krw_,sw_ex, krw_ex, samples);
            initializeTableType(krow_,sw_ex, krow_ex, samples);
            initializeTableType(pcow_,sw_ex, pcow_ex, samples);

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

            extendTable(sg,sg_ex,1);
            extendTable(krg,krg_ex,0);
            extendTable(krog,krog_ex,0);
            extendTable(pcog,pcog_ex,0);

            initializeTableType(krg_,sg_ex, krg_ex, samples);
            initializeTableType(krog_,sg_ex, krog_ex, samples);
            initializeTableType(pcog_,sg_ex, pcog_ex, samples);

            smin_[phase_usage.phase_pos[Vapour]] = sg[0];
            if (std::fabs(sg.back() + swco - 1.0) > 1e-3) {
                OPM_THROW(std::runtime_error, "Gas maximum saturation in SGOF table = " << sg.back() <<
                      ", should equal (1.0 - connate water sat) = " << (1.0 - swco));
            }
            smax_[phase_usage.phase_pos[Vapour]] = sg.back();
            smin_[phase_usage.phase_pos[Vapour]] = sg.front();
            krgmax_ = krg.back();

            sgcr_ = sg.front();
            sogcr_ = 1.0 - sg.back();
            krgr_ = krg.back();
            krorg_ = krg.front();
            for (std::vector<double>::size_type i=1; i<sg.size(); ++i) {
                if (krg[i]> 0.0) {
                   sgcr_ = sg[i-1];
                   krorg_ = krog[i-1];
                   break;
                }
            }
            for (std::vector<double>::size_type i=sg.size()-1; i>=1; --i) {
                if (krog[i-1]> 0.0) {
                   sogcr_ = 1.0 - sg[i];
                   krgr_ = krg[i];
                   break;
                }
            }

        }

        if (phase_usage.phase_used[Vapour] && phase_usage.phase_used[Aqua]) {
            sowcr_ -= smin_[phase_usage.phase_pos[Vapour]];
            sogcr_ -= smin_[phase_usage.phase_pos[Aqua]];
            smin_[phase_usage.phase_pos[Liquid]] = 0.0;
            smax_[phase_usage.phase_pos[Liquid]] = 1.0 - smin_[phase_usage.phase_pos[Aqua]]
                                                       - smin_[phase_usage.phase_pos[Vapour]];  // First entry in SGOF-table supposed to be zero anyway ...
        } else if (phase_usage.phase_used[Aqua]) {
            smin_[phase_usage.phase_pos[Liquid]] = 1.0 - smax_[phase_usage.phase_pos[Aqua]];
            smax_[phase_usage.phase_pos[Liquid]] = 1.0 - smin_[phase_usage.phase_pos[Aqua]];
        } else if (phase_usage.phase_used[Vapour]) {
            smin_[phase_usage.phase_pos[Liquid]] = 1.0 - smax_[phase_usage.phase_pos[Vapour]];
            smax_[phase_usage.phase_pos[Liquid]] = 1.0 - smin_[phase_usage.phase_pos[Vapour]];
        }
    }

    template <class TableType>
    void SatFuncBase<TableType>::updateSatHyst(const double* s,
                                               const EPSTransforms* epst, 
                                               const EPSTransforms* epst_hyst, 
                                               SatHyst* sat_hyst) const
    {    
        if (phase_usage.phase_used[Aqua] && phase_usage.phase_used[Vapour]) { //Water/Oil/Gas
            int opos = phase_usage.phase_pos[BlackoilPhases::Liquid];
            int gpos = phase_usage.phase_pos[BlackoilPhases::Vapour];
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            if (s[opos] > sat_hyst->sow_hyst)
            {
                sat_hyst->sow_hyst = s[opos];
                double _sow_hyst = epst->watoil.scaleSat(sat_hyst->sow_hyst, 1.0-swcr_-smin_[gpos], sowcr_, 1.0-smin_[wpos]-smin_[gpos]);
                double sow_hyst_shifted = epst_hyst->watoil.scaleSatInv(_sow_hyst, 1.0-swcr_-smin_[gpos], sowcr_, 1.0-smin_[wpos]-smin_[gpos]);
                sat_hyst->sow_shift = sow_hyst_shifted - sat_hyst->sow_hyst;
            }
            if (s[gpos] > sat_hyst->sg_hyst)
            {
                sat_hyst->sg_hyst = s[gpos];
                double _sg_hyst = epst->gas.scaleSat(sat_hyst->sg_hyst, 1.0-sogcr_-smin_[wpos], sgcr_, smax_[gpos]);
                double sg_hyst_shifted = epst_hyst->gas.scaleSatInv(_sg_hyst, 1.0-sogcr_-smin_[wpos], sgcr_, smax_[gpos]);
                sat_hyst->sg_shift = sg_hyst_shifted - sat_hyst->sg_hyst;
            }
        } else if (phase_usage.phase_used[Aqua]) { //Water/oil
            int opos = phase_usage.phase_pos[BlackoilPhases::Liquid];
            int wpos = this->phase_usage.phase_pos[BlackoilPhases::Aqua];
            if (s[opos] > sat_hyst->sow_hyst)
            {
                sat_hyst->sow_hyst = s[opos];
                double _sow_hyst = epst->watoil.scaleSat(sat_hyst->sow_hyst, 1.0-swcr_, sowcr_, 1.0-smin_[wpos]);
                double sow_hyst_shifted = epst_hyst->watoil.scaleSatInv(_sow_hyst, 1.0-swcr_, sowcr_, 1.0-smin_[wpos]);
                sat_hyst->sow_shift = sow_hyst_shifted - sat_hyst->sow_hyst;
            }
        } else if (phase_usage.phase_used[Vapour]) {//Gas/Oil
            int gpos = phase_usage.phase_pos[BlackoilPhases::Vapour];
            if (s[gpos] > sat_hyst->sg_hyst)
            {
                sat_hyst->sg_hyst = s[gpos];
                double _sg_hyst = epst->gas.scaleSat(sat_hyst->sg_hyst, 1.0-sogcr_, sgcr_, smax_[gpos]);
                double sg_hyst_shifted = epst_hyst->gas.scaleSatInv(_sg_hyst, 1.0-sogcr_, sgcr_, smax_[gpos]);
                sat_hyst->sg_shift = sg_hyst_shifted - sat_hyst->sg_hyst;
            }
        }
    }

    template <class TableType>
    void SatFuncBase<TableType>::extendTable(const std::vector<double>& xv,
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

} // namespace Opm
#endif // SATFUNCBASE_HPP
