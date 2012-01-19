//===========================================================================
//                                                                           
// File: EclipseGridInspector.h                                              
//                                                                           
// Created: Mon Jun  2 09:46:08 2008                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//                                                                           
// Revision: $Id: EclipseGridInspector.h,v 1.2 2008/08/18 14:16:12 atgeirr Exp $
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

#ifndef SINTEF_ECLIPSEGRIDINSPECTOR_HEADER
#define SINTEF_ECLIPSEGRIDINSPECTOR_HEADER

#include <vector>
#include <tr1/array>


namespace Opm
{

/**
   @brief A class for inspecting the contents of an eclipse file.

   Given an EclipseGridParser that has successfully read en Eclipse
   .grdecl-file, this class may be used to answer certain queries about
   its contents.

   @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
   @date 2008/06/02 09:46:08
*/

class EclipseGridParser;

class EclipseGridInspector
{
public:
    /// Constructor taking a parser as argument.
    /// The parser must already have read an Eclipse file.
    EclipseGridInspector(const EclipseGridParser& parser);

    /// Assuming that the pillars are vertical, compute the
    /// volume of the cell given by logical coordinates (i, j, k).
    double cellVolumeVerticalPillars(int i, int j, int k) const;

    /// Assuming that the pillars are vertical, compute the
    /// volume of the cell given by the cell index
    double cellVolumeVerticalPillars(int cell_idx) const;

    /// Compute the average dip in x- and y-direction of the
    /// cell tops and bottoms relative to the xy-plane
    std::pair<double,double> cellDips(int i, int j, int k) const;
    std::pair<double,double> cellDips(int cell_idx) const;

    // Convert global cell index to logical ijk-coordinates
    std::tr1::array<int, 3> cellIdxToLogicalCoords(int cell_idx) const;
    
    /// Returns a vector with the outer limits of grid (in the grid's unit).
    /// The vector contains [xmin, xmax, ymin, ymax, zmin, zmax], as 
    /// read from COORDS and ZCORN
    std::tr1::array<double, 6> getGridLimits() const;

    /// Returns the extent of the logical cartesian grid
    /// as number of cells in the (i, j, k) directions.
    std::tr1::array<int, 3> gridSize() const;

    /// Returns the eight z-values associated with a given cell.
    /// The ordering is such that i runs fastest. That is, with
    /// L = low and H = high:
    /// {LLL, HLL, LHL, HHL, LLH, HLH, LHH, HHH }.
    std::tr1::array<double, 8> cellZvals(int i, int j, int k) const;

private:
    const EclipseGridParser& parser_;
    int logical_gridsize_[3];
    void checkLogicalCoords(int i, int j, int k) const;
};

} // namespace Opm

#endif // SINTEF_ECLIPSEGRIDINSPECTOR_HEADER

