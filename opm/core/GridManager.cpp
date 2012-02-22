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

#include <opm/core/GridManager.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/cart_grid.h>
#include <opm/core/utility/cpgpreprocess/cgridinterface.h>



namespace Opm
{



    /// Construct a 3d corner-point grid from a deck.
    GridManager::GridManager(const Opm::EclipseGridParser& deck)
    {
	// Extract data from deck.
	const std::vector<double>& zcorn = deck.getFloatingPointValue("ZCORN");
	const std::vector<double>& coord = deck.getFloatingPointValue("COORD");
	const std::vector<int>& actnum = deck.getIntegerValue("ACTNUM");
	std::vector<int> dims;
	if (deck.hasField("DIMENS")) {
	    dims = deck.getIntegerValue("DIMENS");
	} else if (deck.hasField("SPECGRID")) {
	    dims = deck.getSPECGRID().dimensions;
	} else {
	    THROW("Deck must have either DIMENS or SPECGRID.");
	}

	// Collect in input struct for preprocessing.
	struct grdecl grdecl;
	grdecl.zcorn = &zcorn[0];
	grdecl.coord = &coord[0];
	grdecl.actnum = &actnum[0];
	grdecl.dims[0] = dims[0];
	grdecl.dims[1] = dims[1];
	grdecl.dims[2] = dims[2];

	// Process and compute.
	ug_ = preprocess(&grdecl, 0.0);
	if (!ug_) {
	    THROW("Failed to construct grid.");
	}
	compute_geometry(ug_);
    }




    /// Construct a 2d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny)
    {
	ug_ = create_cart_grid_2d(nx, ny);
	if (!ug_) {
	    THROW("Failed to construct grid.");
	}
    }




    /// Construct a 3d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny, int nz)
    {
	ug_ = create_cart_grid_3d(nx, ny, nz);
	if (!ug_) {
	    THROW("Failed to construct grid.");
	}
    }




    /// Construct a 3d cartesian grid with cells of size [dx, dy, dz].
    GridManager::GridManager(int nx, int ny, int nz,
		double dx, double dy, double dz)
    {
	ug_ = create_hexa_grid_3d(nx, ny, nz, dx, dy, dz);
	if (!ug_) {
	    THROW("Failed to construct grid.");
	}
    }




    /// Destructor.
    GridManager::~GridManager()
    {
	free_grid(ug_);
    }




    /// Access the managed UnstructuredGrid.
    /// The method is named similarly to c_str() in std::string,
    /// to make it clear that we are returning a C-compatible struct.
    const UnstructuredGrid* GridManager::c_grid() const
    {
	return ug_;
    }




} // namespace Opm
