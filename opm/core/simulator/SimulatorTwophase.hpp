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

#ifndef OPM_SIMULATORTWOPHASE_HEADER_INCLUDED
#define OPM_SIMULATORTWOPHASE_HEADER_INCLUDED

#include <boost/scoped_ptr.hpp>

struct UnstructuredGrid;
struct Wells;

namespace Opm
{
    namespace parameter { class ParameterGroup; }
    class SimulatorTimer;
    class TwophaseState;
    class WellState;

    /// Class collecting all necessary components for a two-phase simulation.
    class SimulatorTwophase
    {
    public:
        /// Default constructor.
        SimulatorTwophase();

        /// Initialise from parameters.
        void init(const parameter::ParameterGroup& param);

        /// Run the simulation.
        /// \param[in]  timer       governs the requested reporting timesteps
        /// \param[in]  wells       data structure for wells
        /// \param[out] state       state of reservoir: pressure, fluxes
        /// \param[out] well_state  state of wells: bhp, perforation rates
        void run(const SimulatorTimer& timer,
                 const Wells& wells,
                 TwophaseState& state,
                 WellState& well_state);

    private:
        class Impl;
        boost::scoped_ptr<Impl> pimpl_;
    };

} // namespace Opm

#endif // OPM_SIMULATORTWOPHASE_HEADER_INCLUDED
