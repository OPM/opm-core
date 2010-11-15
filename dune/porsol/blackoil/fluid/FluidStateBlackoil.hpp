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


#include "FluidSystemBlackoil.hpp"


namespace Opm
{
    /*!
     * \brief Calculates the phase state from the primary variables in the
     *        blackoil model.
     */
    struct FluidStateBlackoil : public BlackoilDefs
    {
        typedef FluidSystemBlackoil FluidSystem;
        typedef double Scalar;
        typedef std::tr1::array<Scalar, numComponents + numPhases> PrimaryVariables; // Surface volumes and phase pressures.
        typedef FluidMatrixInteractionBlackoil<Scalar> MaterialLaw;
        typedef typename MaterialLaw::Params MaterialLawParams;

        Scalar temperature_;
        Scalar surface_volume_[numComponents];
        Scalar phase_pressure_[numPhases];
        Scalar phase_volume_[numPhases];
        Scalar saturation_[numPhases];


        /*!
         * \brief Update the phase state from the primary variables.
         */
        void update(const PrimaryVariables& primary_vars,
                    const MaterialLawParams& material_params,
                    Scalar temperature)
        {
            // Set the temperature.
            temperature_ = temperature;

            // Set the surface volumes.
            for (int i = 0; i < numComponents; ++i) {
                surface_volume_[i] = primary_vars[i];
            }

            // Set the phase pressures.
            for (int i = 0; i < numPhases; ++i) {
                phase_pressure_[i] = primary_vars[numComponents + i];
            }

            // Compute phase volumes by the fluid system rules.
            FluidSystem::computeEquilibrium(*this);
        }
    };

} // end namespace Opm


#endif // OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED
