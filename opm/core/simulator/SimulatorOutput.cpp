/*
  Copyright (c) 2013 Uni Research AS

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

#include "SimulatorOutput.hpp"

// we need complete definitions for these types
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/io/OutputWriter.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>

#include <numeric> // partial_sum

using namespace Opm;

SimulatorOutputBase::SimulatorOutputBase (
        const parameter::ParameterGroup& params,
        std::shared_ptr <EclipseGridParser> parser,
        std::shared_ptr <SimulatorTimer> timer,
        std::shared_ptr <BlackoilState> state,
        std::shared_ptr <WellState> wellState)

    // store all parameters passed into the object, making them curried
    // parameters to the writeOutput function.
    : parser_         (parser   )
    , timer_          (timer    )
    , reservoirState_ (state    )
    , wellState_      (wellState)

    // process parameters into a writer. we don't setup a new chain in
    // every timestep!
    , writer_ (std::move (OutputWriter::create (params, parser_)))

    // always start from the first timestep
    , next_ (0) {

    // make a list of times to dump. since the original list are relative
    // timesteps, we make a list of accumulated such to compare with
    // current time.
    const std::vector <double>& tstep = parser->getTSTEP ().tstep_;
    times_.resize (tstep.size (), 0.);
    std::partial_sum (tstep.begin(), tstep.end(), times_.begin());

    // write the static initialization files, even before simulation starts
    writer_->writeInit (*timer);
}

SimulatorOutputBase::operator std::function <void ()> () {
    // return (a pointer to) the writeOutput() function as an object
    // which can be passed to the event available from the simulator
    return std::bind (&SimulatorOutputBase::writeOutput, std::ref (*this));
}

void
SimulatorOutputBase::writeOutput () {
    // write output for all timesteps that have passed (usually only
    // one, unless you have entered the same time twice, or the simulator
    // doesn't honor the TSTEP settings)
    for (;
         next_ < times_.size () && times_[next_] <= timer_->currentTime ();
         ++next_) {

        // make sure the simulator has spilled all necessary internal
        // state. notice that this calls *our* sync, which is overridden
        // in the template companion to call the simulator
        sync ();

        // relay the request to the handlers (setup in the constructor
        // from parameters)
        writer_->writeTimeStep (*timer_, *reservoirState_, *wellState_);
    }
}

void
SimulatorOutputBase::sync () {
    // no-op in base class (overridden by simulator-specific template)
}
