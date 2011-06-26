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

#ifndef OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED
#define OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED


#include "BlackoilDefs.hpp"


namespace Opm
{
    /*!
     * \brief Fluid state for a black oil model.
     */
    struct FluidStateBlackoil : public BlackoilDefs
    {
        Scalar temperature_;
        CompVec surface_volume_;
        PhaseVec phase_pressure_;
        PhaseVec phase_volume_density_;
        Scalar total_phase_volume_density_;
        PhaseVec formation_volume_factor_;
        PhaseVec solution_factor_;
        PhaseToCompMatrix phase_to_comp_; // RB^{-1} in Fortran ordering
        PhaseVec saturation_;
        PhaseVec phase_compressibility_;
        Scalar total_compressibility_;
        Scalar experimental_term_;
        PhaseVec viscosity_;
        PhaseVec relperm_;
        PhaseJacobian drelperm_;
        PhaseVec mobility_;
        PhaseJacobian dmobility_;
    };

    /*!
     * \brief Multiple fluid states for a black oil model.
     */
    struct ManyFluidStatesBlackoil : public BlackoilDefs
    {
        std::vector<Scalar> temperature_;
        std::vector<CompVec> surface_volume_;
        std::vector<PhaseVec> phase_pressure_;
        std::vector<PhaseVec> phase_volume_density_;
        std::vector<Scalar> total_phase_volume_density_;
        std::vector<PhaseVec> formation_volume_factor_;
        std::vector<PhaseVec> solution_factor_;
        std::vector<PhaseToCompMatrix> phase_to_comp_; // RB^{-1} in Fortran ordering
        std::vector<PhaseVec> saturation_;
        std::vector<PhaseVec> phase_compressibility_;
        std::vector<Scalar> total_compressibility_;
        std::vector<Scalar> experimental_term_;
        std::vector<PhaseVec> viscosity_;
        std::vector<PhaseVec> relperm_;
        std::vector<PhaseJacobian> drelperm_;
        std::vector<PhaseVec> mobility_;
        std::vector<PhaseJacobian> dmobility_;
    };

    /*!
     * \brief Multiple fluid states for a black oil model.
     */
    struct AllFluidStates : public BlackoilDefs
    {
        // Canonical state variables.
        std::vector<CompVec> surface_volume_density;         // z
        std::vector<PhaseVec> phase_pressure;                // p

        // Variables from PVT functions.
        std::vector<PhaseVec> formation_volume_factor;       // B
        std::vector<PhaseVec> formation_volume_factor_deriv; // dB/dp
        std::vector<PhaseVec> solution_factor;               // R
        std::vector<PhaseVec> solution_factor_deriv;         // dR/dp
        std::vector<PhaseVec> viscosity;                     // mu

        // Variables computed from PVT data.
        // The A matrices are all in Fortran order (or, equivalently,
        // we store the transposes).
        std::vector<PhaseToCompMatrix> state_matrix;         // A' = (RB^{-1})'
        std::vector<PhaseVec> phase_volume_density;          // u
        std::vector<Scalar> total_phase_volume_density;      // sum(u)
        std::vector<PhaseVec> saturation;                    // s = u/sum(u)
        std::vector<PhaseVec> phase_compressibility;         // c
        std::vector<Scalar> total_compressibility;           // cT
        std::vector<Scalar> experimental_term;               // ex = sum(Ai*dA*Ai*z)

        // Variables computed from saturation.
        std::vector<PhaseVec> relperm;                       // kr
        std::vector<PhaseJacobian> relperm_deriv;            // dkr/ds
        std::vector<PhaseVec> mobility;                      // lambda
        std::vector<PhaseJacobian> mobility_deriv;           // dlambda/ds
    };

} // end namespace Opm


#endif // OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED
