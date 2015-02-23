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
}
}
