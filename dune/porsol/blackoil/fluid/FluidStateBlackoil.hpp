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
     * \brief Fluid states for a black oil model.
     */
    struct FluidStateBlackoil : public BlackoilDefs
    {
        Scalar temperature_;
        CompVec surface_volume_;
        PhaseVec phase_pressure_;
        PhaseVec phase_volume_;
        Scalar total_phase_volume_;
        Scalar phase_to_comp_[numPhases*numComponents];
        PhaseVec saturation_;
        PhaseVec phase_compressibility_;
        Scalar total_compressibility_;
        PhaseVec viscosity_;
        PhaseVec relperm_;
        PhaseVec mobility_;
    };

} // end namespace Opm


#endif // OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED
