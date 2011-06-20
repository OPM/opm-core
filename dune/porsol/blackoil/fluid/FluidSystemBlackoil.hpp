/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
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

#ifndef OPM_FLUIDSYSTEMBLACKOIL_HEADER_INCLUDED
#define OPM_FLUIDSYSTEMBLACKOIL_HEADER_INCLUDED

#include "BlackoilPVT.hpp"
#include "BlackoilDefs.hpp"
#include <dune/porsol/common/Matrix.hpp>
#include <dune/common/EclipseGridParser.hpp>
#include <stdexcept>

namespace Opm
{

// Forward declaration needed by associated parameters classes.
template <class ParamT>
class FluidSystemBlackoil;

// Parameters for the black oil fluid system.
class FluidSystemBlackoilParameters
{
public:
    void init(const Dune::EclipseGridParser& ep)
    {
        pvt_.init(ep);
    }
private:
    template <class P>
    friend class FluidSystemBlackoil;
    BlackoilPVT pvt_;
};


/*!
 * \brief A black oil fluid system.
 */
template <class ParamsT = FluidSystemBlackoilParameters>
class FluidSystemBlackoil : public BlackoilDefs
{
public:
    typedef ParamsT Params;

    /*!
     * \brief Initialize system from input.
     */
    static void init(const Dune::EclipseGridParser& ep)
    {
        params().init(ep);
    }

    /*!
     * \brief Return a human-readable phase name.
     */
    static const char* phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case Aqua: return "aqua";
        case Vapour: return "vapour";
        case Liquid: return "liquid";
        default: throw std::logic_error("No such phase.");
        }
    }

    /*!
     * \brief Return a human-readable component name.
     */
    static const char* componentName(int compIdx)
    {
        switch (compIdx) {
        case Water: return "water";
        case Gas: return "gas";
        case Oil: return "oil";
        default: throw std::logic_error("No such component.");
        }
    }


    /*!
     * \brief Return the molar mass of a component [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        throw std::logic_error("molarMass() not implemented.");
    }


    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        throw std::logic_error("phaseDensity() not implemented.");
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its viscosity.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
    {
        throw std::logic_error("phaseViscosity() not implemented.");
    }


    /*!
     * \brief Assuming the surface volumes and the pressures of all
     *        phases are known, compute everything except relperm and
     *        mobility.
     */
    template <class FluidState>
    static void computeEquilibrium(FluidState& fluid_state)
    {
        // Get B and R factors.
        const PhaseVec& p = fluid_state.phase_pressure_;
        const CompVec& z = fluid_state.surface_volume_;
        PhaseVec& B = fluid_state.formation_volume_factor_;
        B[Aqua]   = params().pvt_.B(p[Aqua],   z, Aqua);
        B[Vapour] = params().pvt_.B(p[Vapour], z, Vapour);
        B[Liquid] = params().pvt_.B(p[Liquid], z, Liquid);
        PhaseVec& R = fluid_state.solution_factor_; 
        R[Vapour] = params().pvt_.R(p[Vapour], z, Vapour);
        R[Liquid] = params().pvt_.R(p[Liquid], z, Liquid);
        // Set the A matrix (A = RB^{-1})
        Dune::SharedFortranMatrix A(numComponents, numPhases, &fluid_state.phase_to_comp_[0][0]);
        zero(A);
        A(Water, Aqua) = 1.0/B[Aqua];
        A(Gas, Vapour) = 1.0/B[Vapour];
        A(Gas, Liquid) = R[Liquid]/B[Liquid];
        A(Oil, Vapour) = R[Vapour]/B[Vapour];
        A(Oil, Liquid) = 1.0/B[Liquid];

        // Update phase volumes. This is the same as multiplying with A^{-1}
        PhaseVec& u = fluid_state.phase_volume_density_;
        double detR = 1.0 - R[Vapour]*R[Liquid];
        u[Aqua] = B[Aqua]*z[Water];
        u[Vapour] = B[Vapour]*(z[Gas] - R[Liquid]*z[Oil])/detR;
        u[Liquid] = B[Liquid]*(z[Oil] - R[Vapour]*z[Gas])/detR;
        fluid_state.total_phase_volume_density_ = u[Aqua] + u[Vapour] + u[Liquid];

        // Update saturations.
        double sumu = fluid_state.total_phase_volume_density_;
        PhaseVec& s = fluid_state.saturation_;
        for (int phase = 0; phase < numPhases; ++phase) {
            s[phase] = u[phase]/sumu;
        }

        // Compute viscosities.
        PhaseVec& mu = fluid_state.viscosity_;
        mu[Aqua]   = params().pvt_.getViscosity(p[Aqua],   z, Aqua);
        mu[Vapour] = params().pvt_.getViscosity(p[Vapour], z, Vapour);
        mu[Liquid] = params().pvt_.getViscosity(p[Liquid], z, Liquid);

        // Phase compressibilities.
        PhaseVec& cp = fluid_state.phase_compressibility_;
        double dB[3];
        dB[Aqua]   = params().pvt_.dBdp(p[Aqua],   z, Aqua);
        dB[Vapour] = params().pvt_.dBdp(p[Vapour], z, Vapour);
        dB[Liquid] = params().pvt_.dBdp(p[Liquid], z, Liquid);
        double dR[3]; // Only using 2 of them, though.
        dR[Vapour] = params().pvt_.dRdp(p[Vapour], z, Vapour);
        dR[Liquid] = params().pvt_.dRdp(p[Liquid], z, Liquid);
        // Set the derivative of the A matrix (A = RB^{-1})
        double data_for_dA[numComponents*numPhases];
        Dune::SharedFortranMatrix dA(numComponents, numPhases, data_for_dA);
        zero(dA);
        dA(Water, Aqua) = -dB[Aqua]/(B[Aqua]*B[Aqua]);
        dA(Gas, Vapour) = -dB[Vapour]/(B[Vapour]*B[Vapour]);
        dA(Oil, Liquid) = -dB[Liquid]/(B[Liquid]*B[Liquid]); // Different order than above.
        dA(Gas, Liquid) = dA(Oil, Liquid)*R[Liquid] + dR[Liquid]/B[Liquid];
        dA(Oil, Vapour) = dA(Gas, Vapour)*R[Vapour] + dR[Vapour]/B[Vapour];
        double data_for_Ai[numComponents*numPhases];
        Dune::SharedFortranMatrix Ai(numComponents, numPhases, data_for_Ai);
        // Ai = A; // This does not make a deep copy.
        std::copy(A.data(), A.data() + numComponents*numPhases, Ai.data());
        Dune::invert(Ai);
        double data_for_C[numComponents*numPhases];
        Dune::SharedFortranMatrix C(numComponents, numPhases, data_for_C);
        Dune::prod(Ai, dA, C);
        //CompVec ones(1.0);
        //cp = Dune::prod(C, ones); // Probably C' and not C; we want phasewise compressibilities:
        cp[Aqua] = C(Water, Aqua);
        cp[Liquid] = C(Oil, Liquid) + C(Gas, Liquid);
        cp[Vapour] = C(Gas, Vapour) + C(Oil, Vapour);
        fluid_state.total_compressibility_ = cp*s;

        // Experimental term.
        PhaseVec tmp = prod(Ai, prod(dA, prod(Ai, z)));
        fluid_state.experimental_term_ = tmp[Aqua] + tmp[Liquid] + tmp[Gas];
    }

    /*!
     * \brief Assuming the surface volumes and the pressures of all
     *        phases are known, compute everything except relperm and
     *        mobility.
     */
    template <class ManyFluidStates>
    static void computeManyEquilibria(ManyFluidStates& fluid_state)
    {
        // Get B and R factors, viscosities. Vectorized at lower level.
        const std::vector<PhaseVec>& pv = fluid_state.phase_pressure_;
        const std::vector<CompVec>& zv = fluid_state.surface_volume_;
        std::vector<PhaseVec>& Bv = fluid_state.formation_volume_factor_;
        std::vector<PhaseVec> dBv;
        params().pvt_.dBdp(pv, zv, Bv, dBv);
        std::vector<PhaseVec>& Rv = fluid_state.solution_factor_;
        std::vector<PhaseVec> dRv;
        params().pvt_.dRdp(pv, zv, Rv, dRv);
        std::vector<PhaseVec>& muv = fluid_state.viscosity_;
        params().pvt_.getViscosity(pv, zv, muv);

        // The rest is vectorized in this function.
        int num = pv.size();
        fluid_state.phase_to_comp_.resize(num);
        fluid_state.phase_volume_density_.resize(num);
        fluid_state.total_phase_volume_density_.resize(num);
        fluid_state.phase_to_comp_.resize(num);
        fluid_state.saturation_.resize(num);
        fluid_state.phase_compressibility_.resize(num);
        fluid_state.total_compressibility_.resize(num);
        fluid_state.experimental_term_.resize(num);
#pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            // Convenience vars.
            const PhaseVec& B = Bv[i];
            const PhaseVec& dB = dBv[i];
            const PhaseVec& R = Rv[i];
            const PhaseVec& dR = dRv[i];
            const CompVec& z = zv[i];

            // Set the A matrix (A = RB^{-1})
            Dune::SharedFortranMatrix A(numComponents, numPhases, &fluid_state.phase_to_comp_[i][0][0]);
            zero(A);
            A(Water, Aqua) = 1.0/B[Aqua];
            A(Gas, Vapour) = 1.0/B[Vapour];
            A(Gas, Liquid) = R[Liquid]/B[Liquid];
            A(Oil, Vapour) = R[Vapour]/B[Vapour];
            A(Oil, Liquid) = 1.0/B[Liquid];

            // Update phase volumes. This is the same as multiplying with A^{-1}
            PhaseVec& u = fluid_state.phase_volume_density_[i];
            double detR = 1.0 - R[Vapour]*R[Liquid];
            u[Aqua] = B[Aqua]*z[Water];
            u[Vapour] = B[Vapour]*(z[Gas] - R[Liquid]*z[Oil])/detR;
            u[Liquid] = B[Liquid]*(z[Oil] - R[Vapour]*z[Gas])/detR;
            fluid_state.total_phase_volume_density_[i] = u[Aqua] + u[Vapour] + u[Liquid];

            // Update saturations.
            double sumu = fluid_state.total_phase_volume_density_[i];
            PhaseVec& s = fluid_state.saturation_[i];
            for (int phase = 0; phase < numPhases; ++phase) {
                s[phase] = u[phase]/sumu;
            }

            // Phase compressibilities.
            PhaseVec& cp = fluid_state.phase_compressibility_[i];
            // Set the derivative of the A matrix (A = RB^{-1})
            double data_for_dA[numComponents*numPhases];
            Dune::SharedFortranMatrix dA(numComponents, numPhases, data_for_dA);
            zero(dA);
            dA(Water, Aqua) = -dB[Aqua]/(B[Aqua]*B[Aqua]);
            dA(Gas, Vapour) = -dB[Vapour]/(B[Vapour]*B[Vapour]);
            dA(Oil, Liquid) = -dB[Liquid]/(B[Liquid]*B[Liquid]); // Different order than above.
            dA(Gas, Liquid) = dA(Oil, Liquid)*R[Liquid] + dR[Liquid]/B[Liquid];
            dA(Oil, Vapour) = dA(Gas, Vapour)*R[Vapour] + dR[Vapour]/B[Vapour];
            double data_for_Ai[numComponents*numPhases];
            Dune::SharedFortranMatrix Ai(numComponents, numPhases, data_for_Ai);
            // Ai = A; // This does not make a deep copy.
            std::copy(A.data(), A.data() + numComponents*numPhases, Ai.data());
            Dune::invert(Ai);
            double data_for_C[numComponents*numPhases];
            Dune::SharedFortranMatrix C(numComponents, numPhases, data_for_C);
            Dune::prod(Ai, dA, C);
            //CompVec ones(1.0);
            //cp = Dune::prod(C, ones); // Probably C' and not C; we want phasewise compressibilities:
            cp[Aqua] = C(Water, Aqua);
            cp[Liquid] = C(Oil, Liquid) + C(Gas, Liquid);
            cp[Vapour] = C(Gas, Vapour) + C(Oil, Vapour);
            fluid_state.total_compressibility_[i] = cp*s;

            // Experimental term.
            PhaseVec tmp = prod(Ai, prod(dA, prod(Ai, z)));
            fluid_state.experimental_term_[i] = tmp[Aqua] + tmp[Liquid] + tmp[Gas];
        }
    }

    /*!
     * \brief Returns the activity coefficient of a component in a
     *        phase.
     *
     * We define the activity coefficent \f$\gamma_{\alpha,\kappa}\f$
     * of component \f$\kappa\f$ by the following equation:
     *  \f[ f_\kappa = p_\alpha \gamma_{\alpha,\kappa} \f]
     * where \f$f_\kappa\f$  is the component's fugacity and \f$p_\alpha\f$
     * is the phase' pressure
     *
     * For liquids with very low miscibility this boils down to the
     * inverse Henry constant for the solutes and the partial pressure
     * for the solvent.
     *
     * For ideal gases this is equivalent to the gas pressure, in real
     * gases it is the gas pressure times the component's fugacity
     * coefficient.
     */
    template <class FluidState>
    static Scalar activityCoeff(int phaseIdx,
                                int compIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &state)
    {
        throw std::logic_error("activityCoeff() not implemented.");
        return 0.0;
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficent for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature,
                            Scalar pressure,
                            const FluidState &fluidState)
    {
        throw std::logic_error("diffCoeff() not implemented.");
        return 0.0;
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseEnthalpy(int phaseIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &fluidState)
    {
        throw std::logic_error("phaseEnthalpy() not implemented.");
        return 0.0;
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific internal energy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseInternalEnergy(int phaseIdx,
                                      Scalar temperature,
                                      Scalar pressure,
                                      const FluidState &fluidState)
    {
        throw std::logic_error("phaseInternalEnergy() not implemented.");
        return 0.0;
    }

private:

    static Params& params()
    {
        static Params params; // \TODO: Replace singleton here by something more thread-robust?
        return params;
    }
};

} // namespace Opm

#endif // OPM_FLUIDSYSTEMBLACKOIL_HEADER_INCLUDED
