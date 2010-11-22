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

#ifndef OPM_BLACKOILFLUID_HEADER_INCLUDED
#define OPM_BLACKOILFLUID_HEADER_INCLUDED


#include <dune/porsol/blackoil/fluid/FluidMatrixInteractionBlackoil.hpp>
#include <dune/porsol/blackoil/fluid/FluidSystemBlackoil.hpp>
#include <dune/porsol/blackoil/fluid/FluidStateBlackoil.hpp>
#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/fvector.hh>


namespace Opm
{

    class BlackoilFluid
    {
    public:
        enum { numPhases = 3 };
        enum { numComponents = 3 };
        typedef Dune::FieldVector<double, numPhases> PhaseVec;
        typedef Dune::FieldVector<double, numComponents> CompVec;

        void init(const Dune::EclipseGridParser& parser)
        {
            fmi_params_.init(parser);
            FluidSystemBlackoil<>::init(parser);
        }

        FluidStateBlackoil computeState(PhaseVec phase_pressure, CompVec z)
        {
            FluidStateBlackoil state;
            state.temperature_ = 300;
            state.phase_pressure_ = phase_pressure;
            state.surface_volume_ = z;
            FluidSystemBlackoil<>::computeEquilibrium(state); // Sets everything but relperm and mobility.
            FluidMatrixInteractionBlackoil<double>::kr(state.relperm_, fmi_params_, state.saturation_, state.temperature_);
            for (int phase = 0; phase < numPhases; ++phase) {
                state.mobility_[phase] = state.relperm_[phase]/state.viscosity_[phase];
            }
            return state;
        }

    private:
        FluidMatrixInteractionBlackoilParams<double> fmi_params_;
    };

} // namespace Opm

#endif // OPM_BLACKOILFLUID_HEADER_INCLUDED
