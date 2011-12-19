//===========================================================================
//
// File: Parameter.cpp
//
// Created: Tue Jun  2 19:18:25 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <string>
#include <opm/core/utility/parameters/Parameter.hpp>

namespace Dune {
    namespace parameter {
	std::string
        correct_parameter_tag(const ParameterMapItem& item)
        {
	    std::string tag = item.getTag();
	    if (tag != ID_xmltag__param) {
		std::string error = "The XML tag was '" +
                                    tag + "' but should be '" +
                                    ID_xmltag__param + "'.\n";
		return error;
	    } else {
		return "";
	    }
	}

	std::string
        correct_type(const Parameter& parameter,
                     const std::string& param_type)
        {
	    std::string type = parameter.getType();
	    if ( (type != param_type) &&
                 (type != ID_param_type__cmdline) ) {
		std::string error = "The data was of type '" + type +
                                    "' but should be of type '" +
                                    param_type + "'.\n";
		return error;
	    } else {
		return "";
	    }
	}
    } // namespace parameter
} // namespace Dune
