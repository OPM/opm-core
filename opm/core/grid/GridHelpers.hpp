/*
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2014 Statoil AS
  Copyright 2015

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
#ifndef OPM_CORE_GRIDHELPERS_HEADER_INCLUDED
#define OPM_CORE_GRIDHELPERS_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/range/iterator_range.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>


namespace Opm
{
namespace UgGridHelpers
{

/// \brief Allows viewing a sparse table consisting out of C-array
///
/// This class can be used to convert two int array (like they are
/// in UnstructuredGrid for representing the cell to faces mapping
/// as a sparse table object.
class SparseTableView
{
public:
    class IntRange : public boost::iterator_range<const int*>
    {
    public:
        typedef boost::iterator_range<const int*> BaseRowType;
        typedef BaseRowType::size_type size_type;
        typedef int value_type;

        IntRange(const int* start_arg, const int* end_arg)
            : BaseRowType(start_arg, end_arg)
        {}
    };
    /// \brief The type of the roww.
    typedef boost::iterator_range<const int*> row_type;

    /// \brief Creates a sparse table view
    /// \param data The array with data of the table.
    /// \param offset The offsets of the rows. Row i starts
    ///               at offset[i] and ends a offset[i+1]
    /// \param size   The number of entries/rows of the table
    SparseTableView(int* data, int *offset, std::size_t size_arg)
        : data_(data), offset_(offset), size_(size_arg)
    {}

    /// \brief Get a row of the the table.
    /// \param row The row index.
    /// \return The corresponding row.
    row_type operator[](std::size_t row) const
    {
        assert(row<=size());
        return row_type(data_ + offset_[row], data_ + offset_[row+1]);
    }
    
    /// \brief Get the size of the table.
    /// \return the number rows.
    std::size_t size() const
    {
        return size_;
    }
    
    /// \brief Get the number of non-zero entries.
    std::size_t noEntries() const
    {
        return offset_[size_];
    }

private:
    /// \brief The array with data of the table.
    const int* data_;
    /// \brief offset The offsets of the rows. 
    ///
    /// Row i starts at offset[i] and ends a offset[i+1]
    const int* offset_;
    /// \brief The size, i.e. the number of rows.
    std::size_t size_;
};

/// \brief Get the number of cells of a grid.
int numCells(const UnstructuredGrid& grid);

/// \brief Get the number of faces of a grid.
int numFaces(const UnstructuredGrid& grid);

/// \brief Get the dimensions of a grid
int dimensions(const UnstructuredGrid& grid);

/// \brief Get the number of faces, where each face counts as many times as there are adjacent faces
int numCellFaces(const UnstructuredGrid& grid);

/// \brief Get the cartesion dimension of the underlying structured grid.
const int* cartDims(const UnstructuredGrid& grid);

/// \brief Get the local to global index mapping.
///
/// The global index is the index of the active cell
/// in the underlying structured grid.
const int* globalCell(const UnstructuredGrid& grid);

/// \brief Traits of the cell centroids of a grid.
///
/// This class exports two types: IteratorType, the type of the iterator
/// over the cell centroids, and the ValueTpe, the type of the cell centroid.
/// \tpatam G The type of the grid.
template<class G>
struct CellCentroidTraits
{
};

template<>
struct CellCentroidTraits<UnstructuredGrid>
{
    typedef const double* IteratorType;
    typedef const double* ValueType;
};

/// \brief Get an iterator over the cell centroids positioned at the first cell.
///
/// The return type needs to be usable with the functions increment, and
/// getCoordinate.
CellCentroidTraits<UnstructuredGrid>::IteratorType
beginCellCentroids(const UnstructuredGrid& grid);


/// \brief Get vertical position of cell center ("zcorn" average.)
/// \brief grid The grid.
/// \brief cell_index The index of the specific cell.
double cellCenterDepth(const UnstructuredGrid& grid, int cell_index);


/// \brief Get a coordinate of a specific cell centroid.
/// \brief grid The grid.
/// \brief cell_index The index of the specific cell.
/// \breif coordinate The coordinate index.
double cellCentroidCoordinate(const UnstructuredGrid& grid, int cell_index,
                                 int coordinate);


/// \brief Get the centroid of a cell.
/// \param grid The grid whose cell centroid we query.
/// \param cell_index The index of the corresponding cell.
const double* cellCentroid(const UnstructuredGrid& grid, int cell_index);


/// \brief Get the volume of a cell.
/// \param grid The grid the cell belongs to.
/// \param cell_index The index of the cell.
double cellVolume(const UnstructuredGrid& grid, int cell_index);

/// \brief The mapping of the grid type to type of the iterator over
/// the cell volumes.
///
/// The value of the mapping is stored in nested type IteratorType
/// \tparam T The type of the grid.
template<class T>
struct CellVolumeIteratorTraits
{
};

template<>
struct CellVolumeIteratorTraits<UnstructuredGrid>
{
    typedef const double* IteratorType;
};

    /**
       Will create an EclipseGrid representation (i.e. based on
       ZCORN and COORD) of the current UnstructuredGrid
       instance. When creating the UnstructuredGrid the detailed
       cornerpoint information is discarded, and it is difficult
       to go backwards to recreated ZCORN and COORD.

       The current implementation is based on retaining a copy of the
       zcorn keyword after the Minpvprocessor has modified it.

       We then create a new EclipseGrid instance based on the original
       input grid, but we "replace" the ZCORN and ACTNUM keywords with
       the updated versions.

       If the tolerance in the call to create_grid_cornerpoint( ) is
       finite the grid processing code might collapse cells, the z
       coordinate transformations from this process will *not* be
       correctly represented in the EclipseGrid created by this
       method.
    */

/// \brief Construct an EclipseGrid instance based on the inputGrid, with modifications to
/// zcorn and actnum from the dune UnstructuredGrid.
Opm::EclipseGrid createEclipseGrid(const UnstructuredGrid& grid, const Opm::EclipseGrid& inputGrid );

/// \brief Get an iterator over the cell volumes of a grid positioned at the first cell.
const double* beginCellVolumes(const UnstructuredGrid& grid);

/// \brief Get an iterator over the cell volumes of a grid positioned after the last cell.
const double* endCellVolumes(const UnstructuredGrid& grid);

/// \brief Get the cell centroid of a face.
/// \param grid The grid whose cell centroid we query.
/// \param face_index The index of the corresponding face.
const double* faceCentroid(const UnstructuredGrid& grid, int face_index);


/// \brief Traits of the face centroids of a grid.
///
/// This class exports two types: IteratorType, the type of the iterator
/// over the face centroids, and the ValueTpe, the type of the face centroid.
/// \tpatam G The type of the grid.
template<class G>
struct FaceCentroidTraits
{
};

template<>
struct FaceCentroidTraits<UnstructuredGrid>
{
    typedef const double* IteratorType;
    typedef const double* ValueType;
};

/// \brief Get an iterator over the face centroids positioned at the first cell.
FaceCentroidTraits<UnstructuredGrid>::IteratorType
beginFaceCentroids(const UnstructuredGrid& grid);

/// \brief Get a coordinate of a specific face centroid.
/// \param grid The grid.
/// \param face_index The index of the specific face.
/// \param coordinate The coordinate index.
FaceCentroidTraits<UnstructuredGrid>::ValueType
faceCentroid(const UnstructuredGrid& grid, int face_index);

/// \brief Get the normal of a face.
/// \param grid The grid that the face is part of.
/// \param face_index The index of the face in the grid.
const double* faceNormal(const UnstructuredGrid& grid, int face_index);

/// \brief Get the area of a face
/// \param grid The grid that the face is part of.
/// \param face_index The index of the face in the grid.
double faceArea(const UnstructuredGrid& grid, int face_index);

/// \brief Get Eclipse Cartesian tag of a face
/// \param grid The grid that the face is part of.
/// \param cell_face The face attached to a cell as obtained from cell2Faces()
/// \return 0, 1, 2, 3, 4, 5 for I-, I+, J-, J+, K-, K+
int faceTag(const UnstructuredGrid& grid, boost::iterator_range<const int*>::const_iterator cell_face);

/// \brief Maps the grid type to the associated type of the cell to faces mapping.
///
/// Provides a type named Type.
/// \tparam T The type of the grid.
template<class T>
struct Cell2FacesTraits
{
};

template<>
struct Cell2FacesTraits<UnstructuredGrid>
{
    typedef SparseTableView Type;
};

/// \brief Maps the grid type to the associated type of the face to vertices mapping.
///
/// Provides a type named Type.
/// \tparam T The type of the grid.
template<class T>
struct Face2VerticesTraits
{
};

template<>
struct Face2VerticesTraits<UnstructuredGrid>
{
    typedef SparseTableView Type;
};

/// \brief Get the cell to faces mapping of a grid.
Cell2FacesTraits<UnstructuredGrid>::Type 
cell2Faces(const UnstructuredGrid& grid);

/// \brief Get the face to vertices mapping of a grid.
Face2VerticesTraits<UnstructuredGrid>::Type 
face2Vertices(const UnstructuredGrid& grid);

/// \brief Get the coordinates of a vertex of the grid.
/// \param grid The grid the vertex is part of.
/// \param index The index identifying the vertex.
const double* vertexCoordinates(const UnstructuredGrid& grid, int index);

class FaceCellsProxy
{
public:
    FaceCellsProxy(const UnstructuredGrid& grid)
    : face_cells_(grid.face_cells)
    {}
    int operator()(int face_index, int local_index) const
    {
        return face_cells_[2*face_index+local_index];
    }
private:
    const int* face_cells_;
};

/// \brief Traits of the face to attached cell mappping of a grid.
///
/// Exports the type Type, the type of the mapping
/// \tparam T The type of the grid
template<class T>
struct FaceCellTraits
{};

template<>
struct FaceCellTraits<UnstructuredGrid>
{
    typedef FaceCellsProxy Type;
};

/// \brief Get the face to cell mapping of a grid.
FaceCellTraits<UnstructuredGrid>::Type faceCells(const UnstructuredGrid& grid);

/// \brief Increment an iterator over an array that reresents a dense row-major
///  matrix with dims columns
/// \param cc The iterator.
/// \param i The nzumber of rows to increment
/// \param dim The number of columns of the matrix.
template<class T>
T* increment(T* cc, int i, int dim)
{
    return cc+(i*dim);
}
/// \brief Increment an iterator over an array that reresents a dense row-major
///  matrix with dims columns
/// \param cc The iterator.
/// \param i The nzumber of rows to increment
template<class T>
T increment(const T& t, int i, int)
{
    return t+i;
}

/// \brief Get the i-th corrdinate of a centroid.
/// \param cc The array with the coordinates.
/// \param i The index of the coordinate.
/// \tparam T The type of the coordinate of the centroid.
template<class T>
double getCoordinate(T* cc, int i)
{
    return cc[i];
}

/// \brief Get the i-th corrdinate of an array.
/// \param t The iterator over the centroids
/// \brief i The index of the coordinate.
/// \tparam T The type of the iterator representing the centroid.
/// Its value_type has to provide an operator[] to access the coordinates.
template<class T>
double getCoordinate(T t, int i)
{
    return (*t)[i];
}

} // end namespace UGGridHelpers
} // end namespace OPM
#endif
