//===========================================================================
//
// File: StopWatch.hpp
//
// Created: Thu Jul  2 23:04:17 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_STOPWATCH_HEADER
#define OPENRS_STOPWATCH_HEADER

#include <boost/date_time/posix_time/posix_time.hpp>

namespace Opm
{

    namespace time
    {

	class StopWatch
	{
	public:
	    /// Default constructor. Before the StopWatch is start()-ed,
	    /// it is an error to call anything other than start().
	    StopWatch();

	    /// Starts the StopWatch. It is always legal to call
	    /// start(), even if not stop()-ped.
	    void start();
	    /// Stops the StopWatch. The watch no longer runs, until
	    /// restarted by a call to start().
	    void stop();

	    /// \ret the number of running seconds that have passed
	    /// since last call to start(), secsSinceLast() or
	    /// secsSinceStart()
	    double secsSinceLast();
	    /// \ret the number of running seconds that have passed
	    /// since last call to start().
	    double secsSinceStart();

	private:
	    enum StopWatchState { UnStarted, Running, Stopped };

	    StopWatchState state_;
	    boost::posix_time::ptime start_time_;
	    boost::posix_time::ptime last_time_;
	    boost::posix_time::ptime stop_time_;
	};

    } // namespace time

} // namespace Opm

#endif // OPENRS_STOPWATCH_HEADER
