/*
  Copyright 2014 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2014 Statoil AS

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
#include <boost/range/iterator_range.hpp>

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
        
        IntRange(const int* start, const int* end)
            : BaseRowType(start, end)
        {}
    };
    /// \brief The type of the roww.
    typedef boost::iterator_range<const int*> row_type;

    /// \brief Creates a sparse table view
    /// \param data The array with data of the table.
    /// \param offset The offsets of the rows. Row i starts
    ///               at offset[i] and ends a offset[i+1]
    /// \param size   The number of entries/rows of the table
    SparseTableView(int* data, int *offset, std::size_t size)
        : data_(data), offset_(offset), size_(size)
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

/// \brief Get an iterator over the cell centroids positioned at the first cell.
///
/// The return type needs to be usable with the functions increment, and
/// getCoordinate.
const double* beginCellCentroids(const UnstructuredGrid& grid);

/// \brief Get a coordinate of a specific cell centroid.
/// \brief grid The grid.
/// \brief cell_index The index of the specific cell.
/// \breif coordinate The coordinate index.
double cellCentroidCoordinate(const UnstructuredGrid& grid, int cell_index,
                                 int coordinate);

/// \brief Get an iterator over the face centroids positioned at the first cell.
const double* beginFaceCentroids(const UnstructuredGrid& grid);

/// \brief Get a coordinate of a specific face centroid.
/// \brief grid The grid.
/// \brief face_index The index of the specific face.
/// \breif coordinate The coordinate index.
const double* faceCentroid(const UnstructuredGrid& grid, int face_index);

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

/// \brief Get the cell to faces mapping of a grid.
Cell2FacesTraits<UnstructuredGrid>::Type 
cell2Faces(const UnstructuredGrid& grid);

class FaceCellsProxy
{
public:
    FaceCellsProxy(const UnstructuredGrid& grid)
    : face_cells_(grid.face_cells)
    {}
    int operator()(int cell_index, int local_index)
    {
        return face_cells_[2*cell_index+local_index];
    }
private:
    const int* face_cells_;
};

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
