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

#include <opm/core/grid/GridManager.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/cart_grid.h>
#include <opm/core/grid/cornerpoint_grid.h>
#include <algorithm>
#include <numeric>



namespace Opm
{



    /// Construct a 3d corner-point grid from a deck.
    GridManager::GridManager(const Opm::EclipseGridParser& deck)
    {
        // We accept two different ways to specify the grid.
        //    1. Corner point format.
        //       Requires ZCORN, COORDS, DIMENS or SPECGRID, optionally
        //       ACTNUM, optionally MAPAXES.
        //       For this format, we will verify that DXV, DYV, DZV,
        //       DEPTHZ and TOPS are not present.
        //    2. Tensor grid format.
        //       Requires DXV, DYV, DZV, optionally DEPTHZ or TOPS.
        //       For this format, we will verify that ZCORN, COORDS
        //       and ACTNUM are not present.
        //       Note that for TOPS, we only allow a uniform vector of values.

        if (deck.hasField("ZCORN") && deck.hasField("COORD")) {
            initFromDeckCornerpoint(deck);
        } else if (deck.hasField("DXV") && deck.hasField("DYV") && deck.hasField("DZV")) {
            initFromDeckTensorgrid(deck);
        } else {
            THROW("Could not initialize grid from deck. "
                  "Need either ZCORN + COORD or DXV + DYV + DZV keywords.");
        }
    }




    /// Construct a 2d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny)
    {
        ug_ = create_grid_cart2d(nx, ny, 1.0, 1.0);
        if (!ug_) {
            THROW("Failed to construct grid.");
        }
    }

    GridManager::GridManager(int nx, int ny,double dx, double dy)
    {
        ug_ = create_grid_cart2d(nx, ny, dx, dy);
        if (!ug_) {
            THROW("Failed to construct grid.");
        }
    }


    /// Construct a 3d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny, int nz)
    {
        ug_ = create_grid_cart3d(nx, ny, nz);
        if (!ug_) {
            THROW("Failed to construct grid.");
        }
    }




    /// Construct a 3d cartesian grid with cells of size [dx, dy, dz].
    GridManager::GridManager(int nx, int ny, int nz,
                             double dx, double dy, double dz)
    {
        ug_ = create_grid_hexa3d(nx, ny, nz, dx, dy, dz);
        if (!ug_) {
            THROW("Failed to construct grid.");
        }
    }




    /// Construct a grid from an input file.
    /// The file format used is currently undocumented,
    /// and is therefore only suited for internal use.
    GridManager::GridManager(const std::string& input_filename)
    {
        ug_ = read_grid(input_filename.c_str());
        if (!ug_) {
            THROW("Failed to read grid from file " << input_filename);
        }
    }




    /// Destructor.
    GridManager::~GridManager()
    {
        destroy_grid(ug_);
    }




    /// Access the managed UnstructuredGrid.
    /// The method is named similarly to c_str() in std::string,
    /// to make it clear that we are returning a C-compatible struct.
    const UnstructuredGrid* GridManager::c_grid() const
    {
        return ug_;
    }



    // Construct corner-point grid from deck.
    void GridManager::initFromDeckCornerpoint(const Opm::EclipseGridParser& deck)
    {
        // Extract data from deck.
        // Collect in input struct for preprocessing.
        struct grdecl grdecl = deck.get_grdecl();

        // Process grid.
        ug_ = create_grid_cornerpoint(&grdecl, 0.0);
        if (!ug_) {
            THROW("Failed to construct grid.");
        }
    }


    namespace
    {
        std::vector<double> coordsFromDeltas(const std::vector<double>& deltas)
        {
            std::vector<double> coords(deltas.size() + 1);
            coords[0] = 0.0;
            std::partial_sum(deltas.begin(), deltas.end(), coords.begin() + 1);
            return coords;
        }
    } // anonymous namespace


    // Construct tensor grid from deck.
    void GridManager::initFromDeckTensorgrid(const Opm::EclipseGridParser& deck)
    {
        // Extract logical cartesian size.
        std::vector<int> dims;
        if (deck.hasField("DIMENS")) {
            dims = deck.getIntegerValue("DIMENS");
        } else if (deck.hasField("SPECGRID")) {
            dims = deck.getSPECGRID().dimensions;
        } else {
            THROW("Deck must have either DIMENS or SPECGRID.");
        }

        // Extract coordinates (or offsets from top, in case of z).
        const std::vector<double>& dxv = deck.getFloatingPointValue("DXV");
        const std::vector<double>& dyv = deck.getFloatingPointValue("DYV");
        const std::vector<double>& dzv = deck.getFloatingPointValue("DZV");
        std::vector<double> x = coordsFromDeltas(dxv);
        std::vector<double> y = coordsFromDeltas(dyv);
        std::vector<double> z = coordsFromDeltas(dzv);

        // Check that number of cells given are consistent with DIMENS/SPECGRID.
        if (dims[0] != int(dxv.size())) {
            THROW("Number of DXV data points do not match DIMENS or SPECGRID.");
        }
        if (dims[1] != int(dyv.size())) {
            THROW("Number of DYV data points do not match DIMENS or SPECGRID.");
        }
        if (dims[2] != int(dzv.size())) {
            THROW("Number of DZV data points do not match DIMENS or SPECGRID.");
        }

        // Extract top corner depths, if available.
        const double* top_depths = 0;
        std::vector<double> top_depths_vec;
        if (deck.hasField("DEPTHZ")) {
            const std::vector<double>& depthz = deck.getFloatingPointValue("DEPTHZ");
            if (depthz.size() != x.size()*y.size()) {
                THROW("Incorrect size of DEPTHZ: " << depthz.size());
            }
            top_depths = &depthz[0];
        } else if (deck.hasField("TOPS")) {
            // We only support constant values for TOPS.
            // It is not 100% clear how we best can deal with
            // varying TOPS (stair-stepping grid, or not).
            const std::vector<double>& tops = deck.getFloatingPointValue("TOPS");
            if (std::count(tops.begin(), tops.end(), tops[0]) != int(tops.size())) {
                THROW("We do not support nonuniform TOPS, please use ZCORN/COORDS instead.");
            }
            top_depths_vec.resize(x.size()*y.size(), tops[0]);
            top_depths = &top_depths_vec[0];
        }

        // Construct grid.
        ug_ = create_grid_tensor3d(dxv.size(), dyv.size(), dzv.size(),
                                   &x[0], &y[0], &z[0], top_depths);
        if (!ug_) {
            THROW("Failed to construct grid.");
        }
    }



} // namespace Opm
