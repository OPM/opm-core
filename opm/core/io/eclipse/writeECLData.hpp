/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

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

#ifndef OPM_WRITEECLDATA_HEADER_INCLUDED
#define OPM_WRITEECLDATA_HEADER_INCLUDED


#include <opm/core/utility/DataMap.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <string>

struct UnstructuredGrid;

namespace Opm
{

  // ECLIPSE output for general grids.
  void writeECLData(const UnstructuredGrid& grid,
                    const DataMap& data,
                    const int current_step,
                    const double current_time,
                    const boost::posix_time::ptime& current_date_time,
                    const std::string& output_dir,
                    const std::string& base_name);

} 

#endif 
