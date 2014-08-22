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

#include "config.h"

#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/cart_grid.h>
#include <opm/core/grid/cornerpoint_grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <array>
#include <algorithm>
#include <numeric>

namespace Opm
{
    /// Construct a 3d corner-point grid from a deck.
    GridManager::GridManager(Opm::EclipseGridConstPtr eclipseGrid)
        : ug_(0)
    {
        initFromEclipseGrid(eclipseGrid);
    }


    GridManager::GridManager(Opm::DeckConstPtr deck)
        : ug_(0)
    {
        auto eclipseGrid = std::make_shared<const Opm::EclipseGrid>(deck);
        initFromEclipseGrid(eclipseGrid);
    }



    /// Construct a 2d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny)
    {
        ug_ = create_grid_cart2d(nx, ny, 1.0, 1.0);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }

    GridManager::GridManager(int nx, int ny,double dx, double dy)
    {
        ug_ = create_grid_cart2d(nx, ny, dx, dy);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }


    /// Construct a 3d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny, int nz)
    {
        ug_ = create_grid_cart3d(nx, ny, nz);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }




    /// Construct a 3d cartesian grid with cells of size [dx, dy, dz].
    GridManager::GridManager(int nx, int ny, int nz,
                             double dx, double dy, double dz)
    {
        ug_ = create_grid_hexa3d(nx, ny, nz, dx, dy, dz);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }




    /// Construct a grid from an input file.
    /// The file format used is currently undocumented,
    /// and is therefore only suited for internal use.
    GridManager::GridManager(const std::string& input_filename)
    {
        ug_ = read_grid(input_filename.c_str());
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to read grid from file " << input_filename);
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




    // Construct corner-point grid from EclipseGrid.
    void GridManager::initFromEclipseGrid(Opm::EclipseGridConstPtr eclipseGrid)
    {
        struct grdecl g;
        std::vector<int> actnum;
        std::vector<double> coord;
        std::vector<double> zcorn;
        std::vector<double> mapaxes;

        g.dims[0] = eclipseGrid->getNX();
        g.dims[1] = eclipseGrid->getNY();
        g.dims[2] = eclipseGrid->getNZ();

        eclipseGrid->exportMAPAXES( mapaxes );
        eclipseGrid->exportCOORD( coord );
        eclipseGrid->exportZCORN( zcorn );
        eclipseGrid->exportACTNUM( actnum );

        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.data();
        g.mapaxes = mapaxes.data();

        const double z_tolerance = eclipseGrid->isPinchActive() ?
            eclipseGrid->getPinchThresholdThickness() : 0.0;
        ug_ = create_grid_cornerpoint(&g, z_tolerance);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }




    void GridManager::createGrdecl(Opm::DeckConstPtr deck, struct grdecl &grdecl)
    {
        // Extract data from deck.
        const std::vector<double>& zcorn = deck->getKeyword("ZCORN")->getSIDoubleData();
        const std::vector<double>& coord = deck->getKeyword("COORD")->getSIDoubleData();
        const int* actnum = NULL;
        if (deck->hasKeyword("ACTNUM")) {
            actnum = &(deck->getKeyword("ACTNUM")->getIntData()[0]);
        }

        std::array<int, 3> dims;
        if (deck->hasKeyword("DIMENS")) {
            Opm::DeckKeywordConstPtr dimensKeyword = deck->getKeyword("DIMENS");
            dims[0] = dimensKeyword->getRecord(0)->getItem(0)->getInt(0);
            dims[1] = dimensKeyword->getRecord(0)->getItem(1)->getInt(0);
            dims[2] = dimensKeyword->getRecord(0)->getItem(2)->getInt(0);
        } else if (deck->hasKeyword("SPECGRID")) {
            Opm::DeckKeywordConstPtr specgridKeyword = deck->getKeyword("SPECGRID");
            dims[0] = specgridKeyword->getRecord(0)->getItem(0)->getInt(0);
            dims[1] = specgridKeyword->getRecord(0)->getItem(1)->getInt(0);
            dims[2] = specgridKeyword->getRecord(0)->getItem(2)->getInt(0);
        } else {
            OPM_THROW(std::runtime_error, "Deck must have either DIMENS or SPECGRID.");
        }

        // Collect in input struct for preprocessing.

        grdecl.zcorn = &zcorn[0];
        grdecl.coord = &coord[0];
        grdecl.actnum = actnum;
        grdecl.dims[0] = dims[0];
        grdecl.dims[1] = dims[1];
        grdecl.dims[2] = dims[2];

        if (deck->hasKeyword("MAPAXES")) {
            Opm::DeckKeywordConstPtr mapaxesKeyword = deck->getKeyword("MAPAXES");
            Opm::DeckRecordConstPtr mapaxesRecord = mapaxesKeyword->getRecord(0);

            // memleak alert: here we need to make sure that C code
            // can properly take ownership of the grdecl.mapaxes
            // object. if it is not freed, it will result in a
            // memleak...
            double *cWtfMapaxes = static_cast<double*>(malloc(sizeof(double)*mapaxesRecord->size()));
            for (unsigned i = 0; i < mapaxesRecord->size(); ++i)
                cWtfMapaxes[i] = mapaxesRecord->getItem(i)->getSIDouble(0);
            grdecl.mapaxes = cWtfMapaxes;
        } else
            grdecl.mapaxes = NULL;

    }

} // namespace Opm
