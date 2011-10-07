//===========================================================================
//
// File: ParameterXML.hpp
//
// Created: Tue Jun  2 18:55:59 2009
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

#ifndef OPENRS_PARAMETERXML_HEADER
#define OPENRS_PARAMETERXML_HEADER

#include <dune/common/param/ParameterGroup.hpp>

namespace Dune {
    namespace parameter {
	/// \brief This function fills a ParameterGroup with a  XML file.
	///
	/// NOTE: The function throws if there is an error during reading
	/// of the XML file.
	///
	/// \retval pg is the ParameterGroup to be filled.
	/// \param filename is the name of an XML file.
	void fill_xml(ParameterGroup& pg, const std::string filename);
    } // namespace parameter
} // namespace Dune

#endif // OPENRS_PARAMETERXML_HEADER
