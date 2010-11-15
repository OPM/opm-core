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

#ifndef OPM_FLUIDMATRIXINTERACTIONBLACKOIL_HEADER_INCLUDED
#define OPM_FLUIDMATRIXINTERACTIONBLACKOIL_HEADER_INCLUDED

#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/Units.hpp>
#include <dune/porsol/common/UniformTableLinear.hpp>
#include <dune/porsol/common/buildUniformMonotoneTable.hpp>
#include <iostream>
#include <stdexcept>

namespace Opm
{

// Forward declaration needed by associated parameters class.
template <class ScalarT, class ParamsT>
class FluidMatrixInteractionBlackoil;

template <class ScalarT>
class FluidMatrixInteractionBlackoilParams
{
public:
    typedef ScalarT Scalar;
    void init(const Dune::EclipseGridParser& ep)
    {
        // Extract input data.
        const Dune::SGOF::table_t& sgof_table = ep.getSGOF().sgof_;
        const Dune::SWOF::table_t& swof_table = ep.getSWOF().swof_;
        if (sgof_table.size() != 1 || swof_table.size() != 1) {
            std::cerr << "We must have exactly one SWOF and one SGOF table (at the moment).\n";
            throw std::logic_error("Not implemented");
        }
        const std::vector<double>& sw = swof_table[0][0];
        const std::vector<double>& krw = swof_table[0][1];
        const std::vector<double>& krow = swof_table[0][2];
        const std::vector<double>& pcow_raw = swof_table[0][3];
        const std::vector<double>& sg = sgof_table[0][0];
        const std::vector<double>& krg = sgof_table[0][1];
        const std::vector<double>& krog = sgof_table[0][2];
        const std::vector<double>& pcog_raw = sgof_table[0][3];

        // Create tables for krw, krow, krg and krog.
        int samples = 200;
        buildUniformMonotoneTable(sw, krw,  samples, krw_);
        buildUniformMonotoneTable(sw, krow, samples, krow_);
        buildUniformMonotoneTable(sg, krg,  samples, krg_);
        buildUniformMonotoneTable(sg, krog, samples, krog_);
        for (int i = 0; i < int(krw.size()); ++i) {
            if (krw[i] > 0.0) {
                krocw_ = krow[i];
                break;
            }
        }

        // Create tables for pcow and pcog.
        // We must convert the pressures depending on units.
        double pressure_unit = ep.units().pressure;
        int numw = sw.size();
        std::vector<double> pcow(numw);
        for (int i = 0; i < numw; ++i) {
            pcow[i] = Dune::unit::convert::from(pcow_raw[i], pressure_unit);
        }
        int numg = sg.size();
        std::vector<double> pcog(numg);
        for (int i = 0; i < numg; ++i) {
            pcog[i] = Dune::unit::convert::from(pcog_raw[i], pressure_unit);
        }
        buildUniformMonotoneTable(sw, pcow, samples, pcow_);
        buildUniformMonotoneTable(sg, pcog, samples, pcog_);
    }

private:
    template <class S, class P>
    friend class FluidMatrixInteractionBlackoil;

    Dune::utils::UniformTableLinear<Scalar> krw_;
    Dune::utils::UniformTableLinear<Scalar> krow_;
    Dune::utils::UniformTableLinear<Scalar> pcow_;
    Dune::utils::UniformTableLinear<Scalar> krg_;
    Dune::utils::UniformTableLinear<Scalar> krog_;
    Dune::utils::UniformTableLinear<Scalar> pcog_;
    Scalar krocw_; // = krow_(s_wc)
};


/*!
 * \ingroup material
 *
 * \brief Capillary pressures and relative permeabilities for a black oil system.
 */
template <class ScalarT, class ParamsT = FluidMatrixInteractionBlackoilParams<ScalarT> >
class FluidMatrixInteractionBlackoil
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;
    enum { numPhases = 3 };
    enum PhaseIndex { Aqua = 0, Vapour = 1, Liquid = 2 };

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param Swe Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    template <class pcContainerT, class SatContainerT>
    static void pC(pcContainerT &pc,
                   const Params &params, 
                   const SatContainerT &saturations,
                   Scalar /*temperature*/)
    {
        Scalar sw = saturations[Aqua];
        Scalar sg = saturations[Vapour];
        pc[Liquid] = 0.0;
        pc[Aqua] = params.pcow_(sw);
        pc[Vapour] = params.pcog_(sg);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     \f]
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    template <class SatContainerT, class pcContainerT>
    static void S(SatContainerT &saturations,
                  const Params &params, 
                  const pcContainerT &pc,
                  Scalar /*temperature*/)
    {
        std::cerr << "FluidMatrixInteractionBlackoil::S() is not implemented yet\n";
        throw std::logic_error("Not implemented");
    }


    /*!
     * \brief The relative permeability of all phases.
     */
    template <class krContainerT, class SatContainerT>
    static void kr(krContainerT &kr,
                   const Params &params, 
                   const SatContainerT &saturations,
                   Scalar /*temperature*/)
    {
        // Stone-II relative permeability model.
        Scalar sw = saturations[Aqua];
        Scalar sg = saturations[Vapour];
        Scalar so = saturations[Liquid];
        Scalar krw = params.krw_(sw);
        Scalar krg = params.krg_(sg);
        Scalar krow = params.krow_(so);
        Scalar krog = params.krog_(sg);
        Scalar krocw = params.krocw_;
        kr[Aqua] = krw;
        kr[Vapour] = krg;
        kr[Liquid] = krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
        if (kr[Liquid] < 0.0) {
            kr[Liquid] = 0.0;
        }
    }
};

} // namespace Opm




#endif // OPM_FLUIDMATRIXINTERACTIONBLACKOIL_HEADER_INCLUDED
