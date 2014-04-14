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

SparseTableView cell2Faces(const UnstructuredGrid& grid)
{
    return SparseTableView(grid.cell_faces, grid.cell_facepos, numCells(grid));
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
