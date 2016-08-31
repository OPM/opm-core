/*
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services.
  Copyright 2014 Statoil AS
  Copyright 2015 NTNU

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
#include <opm/core/grid/GridHelpers.hpp>
namespace Opm
{
namespace UgGridHelpers
{
int numCells(const UnstructuredGrid& grid)
{
    return grid.number_of_cells;
}

int numFaces(const UnstructuredGrid& grid)
{
    return grid.number_of_faces;
}
int dimensions(const UnstructuredGrid& grid)
{
    return grid.dimensions;
}
int numCellFaces(const UnstructuredGrid& grid)
{
    return grid.cell_facepos[grid.number_of_cells];
}

const int* globalCell(const UnstructuredGrid& grid)
{
    return grid.global_cell;
}

const int* cartDims(const UnstructuredGrid& grid)
{
    return grid.cartdims;
}
                    
const double* beginCellCentroids(const UnstructuredGrid& grid)
{
    return grid.cell_centroids;
}

double cellCenterDepth(const UnstructuredGrid& grid, int cell_index) 
{
    // This method is an alternative to the method cellCentroidCoordinate(...) below.
    // The cell center depth is computed as a raw average of cell corner depths.
    // For cornerpoint grids, this is likely to give slightly different depths that seem
    // to agree with eclipse.
    assert(grid.dimensions == 3);
    const int nd = 3; // Assuming 3-dimensional grid ...
    const int nv = 8; // Assuming 2*4 vertices ...
    double zz = 0.0;
    // Traverse the bottom and top cell-face
    for (int i=grid.cell_facepos[cell_index+1]-2; i<grid.cell_facepos[cell_index+1]; ++i) {
        // Traverse the vertices associated with each face
        assert(grid.face_nodepos[grid.cell_faces[i]+1] - grid.face_nodepos[grid.cell_faces[i]] == nv/2);
        for (int j=grid.face_nodepos[grid.cell_faces[i]]; j<grid.face_nodepos[grid.cell_faces[i]+1]; ++j) {
            zz += (grid.node_coordinates+nd*(grid.face_nodes[j]))[nd-1];
        }
    }
    return zz/nv;
}

double cellCentroidCoordinate(const UnstructuredGrid& grid, int cell_index,
                                 int coordinate)
{
    return grid.cell_centroids[grid.dimensions*cell_index+coordinate];
}

const double*
cellCentroid(const UnstructuredGrid& grid, int cell_index)
{
    return grid.cell_centroids+(cell_index*grid.dimensions);
}

const double* beginCellVolumes(const UnstructuredGrid& grid)
{
    return grid.cell_volumes;
}
const double* endCellVolumes(const UnstructuredGrid& grid)
{
    return grid.cell_volumes+numCells(grid);
}

const double* beginFaceCentroids(const UnstructuredGrid& grid)
{
    return grid.face_centroids;
}

const double* faceCentroid(const UnstructuredGrid& grid, int face_index)
{
    return grid.face_centroids+face_index*grid.dimensions;
}

const double* faceNormal(const UnstructuredGrid& grid, int face_index)
{
    return grid.face_normals+face_index*grid.dimensions;
}

double faceArea(const UnstructuredGrid& grid, int face_index)
{
    return grid.face_areas[face_index];
}

int faceTag(const UnstructuredGrid& grid,
            boost::iterator_range<const int*>::const_iterator face)
{
    return grid.cell_facetag[face-cell2Faces(grid)[0].begin()];
}

SparseTableView cell2Faces(const UnstructuredGrid& grid)
{
    return SparseTableView(grid.cell_faces, grid.cell_facepos, numCells(grid));
}

SparseTableView face2Vertices(const UnstructuredGrid& grid)
{
    return SparseTableView(grid.face_nodes, grid.face_nodepos, numFaces(grid));
}

const double* vertexCoordinates(const UnstructuredGrid& grid, int index)
{
    return grid.node_coordinates+dimensions(grid)*index;
}

double cellVolume(const UnstructuredGrid& grid, int cell_index)
{
    return grid.cell_volumes[cell_index];
}

FaceCellTraits<UnstructuredGrid>::Type faceCells(const UnstructuredGrid& grid)
{
    return FaceCellsProxy(grid);
}


Opm::EclipseGrid createEclipseGrid(const UnstructuredGrid& grid, const Opm::EclipseGrid& inputGrid ) {
    const int * dims = UgGridHelpers::cartDims( grid );

    if ((inputGrid.getNX( ) == static_cast<size_t>(dims[0])) &&
        (inputGrid.getNY( ) == static_cast<size_t>(dims[1])) &&
        (inputGrid.getNZ( ) == static_cast<size_t>(dims[2]))) {
        std::vector<int> updatedACTNUM;
        const int* global_cell = UgGridHelpers::globalCell( grid );

        if (global_cell) {
            updatedACTNUM.assign( inputGrid.getCartesianSize( ) , 0 );
            for (int c = 0; c < numCells( grid ); c++) {
                updatedACTNUM[global_cell[c]] = 1;
            }
        }

        return Opm::EclipseGrid( inputGrid, grid.zcorn, updatedACTNUM );
    } else {
        throw std::invalid_argument("Size mismatch - dimensions of inputGrid argument and current UnstructuredGrid instance disagree");
    }
}

}
}
