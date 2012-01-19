//===========================================================================
//
// File: SparseTable.hpp
//
// Created: Fri Apr 24 09:50:27 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPENRS_SPARSETABLE_HEADER
#define OPENRS_SPARSETABLE_HEADER

#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/range/iterator_range.hpp>
#include "ErrorMacros.hpp"

#include <ostream>

namespace Opm
{

    /// A SparseTable stores a table with rows of varying size
    /// as efficiently as possible.
    /// It is supposed to behave similarly to a vector of vectors.
    /// Its behaviour is similar to compressed row sparse matrices.
    template <typename T>
    class SparseTable
    {
    public:
        /// Default constructor. Yields an empty SparseTable.
        SparseTable()
        {
        }

        /// A constructor taking all the data for the table and row sizes.
        /// \param data_beg The start of the table data.
        /// \param data_end One-beyond-end of the table data.
        /// \param rowsize_beg The start of the row length data.
        /// \param rowsize_end One beyond the end of the row length data.
        template <typename DataIter, typename IntegerIter>
        SparseTable(DataIter data_beg, DataIter data_end,
                    IntegerIter rowsize_beg, IntegerIter rowsize_end)
            : data_(data_beg, data_end)
        {
	    setRowStartsFromSizes(rowsize_beg, rowsize_end);
        }


        /// Sets the table to contain the given data, organized into
	/// rows as indicated by the given row sizes.
        /// \param data_beg The start of the table data.
        /// \param data_end One-beyond-end of the table data.
        /// \param rowsize_beg The start of the row length data.
        /// \param rowsize_end One beyond the end of the row length data.
        template <typename DataIter, typename IntegerIter>
        void assign(DataIter data_beg, DataIter data_end,
                    IntegerIter rowsize_beg, IntegerIter rowsize_end)
        {
	    data_.assign(data_beg, data_end);
	    setRowStartsFromSizes(rowsize_beg, rowsize_end);
        }


        /// Request storage for table of given size.
        /// \param rowsize_beg Start of row size data.
        /// \param rowsize_end One beyond end of row size data.
        template <typename IntegerIter>
        void allocate(IntegerIter rowsize_beg, IntegerIter rowsize_end)
        {
            typedef typename std::vector<T>::size_type sz_t;

            sz_t ndata = std::accumulate(rowsize_beg, rowsize_end, sz_t(0));
            data_.resize(ndata);
            setRowStartsFromSizes(rowsize_beg, rowsize_end);
        }


        /// Appends a row to the table.
        template <typename DataIter>
        void appendRow(DataIter row_beg, DataIter row_end)
        {
            data_.insert(data_.end(), row_beg, row_end);
            if (row_start_.empty()) {
                row_start_.reserve(2);
                row_start_.push_back(0);
            }
            row_start_.push_back(data_.size());
        }

        /// True if the table contains no rows.
        bool empty() const
        {
            return row_start_.empty();
        }

        /// Returns the number of rows in the table.
        int size() const
        {
            return empty() ? 0 : row_start_.size() - 1;
        }

        /// Allocate storage for table of expected size
        void reserve(int exptd_nrows, int exptd_ndata)
        {
            row_start_.reserve(exptd_nrows + 1);
            data_.reserve(exptd_ndata);
        }

        /// Swap contents for other SparseTable<T>
        void swap(SparseTable<T>& other)
        {
            row_start_.swap(other.row_start_);
            data_.swap(other.data_);
        }

        /// Returns the number of data elements.
        int dataSize() const
        {
            return data_.size();
        }

        /// Returns the size of a table row.
        int rowSize(int row) const
        {
            ASSERT(row >= 0 && row < size());
            return row_start_[row + 1] - row_start_[row];
        }

        /// Makes the table empty().
        void clear()
        {
            data_.clear();
            row_start_.clear();
        }

        /// Defining the row type, returned by operator[].
        typedef boost::iterator_range<const T*> row_type;
        typedef boost::iterator_range<T*>       mutable_row_type;

        /// Returns a row of the table.
        row_type operator[](int row) const
        {
            ASSERT(row >= 0 && row < size());
            const T* start_ptr = data_.empty() ? 0 : &data_[0];
            return row_type(start_ptr + row_start_[row], start_ptr + row_start_[row + 1]);
        }

        /// Returns a mutable row of the table.
        mutable_row_type operator[](int row)
        {
            ASSERT(row >= 0 && row < size());
            T* start_ptr = data_.empty() ? 0 : &data_[0];
            return mutable_row_type(start_ptr + row_start_[row], start_ptr + row_start_[row + 1]);
        }

        /// Equality.
        bool operator==(const SparseTable& other) const
        {
            return data_ == other.data_ && row_start_ == other.row_start_;
        }

        template<class charT, class traits>
        void print(std::basic_ostream<charT, traits>& os) const
        {
            os << "Number of rows: " << size() << '\n';

            os << "Row starts = [";
            std::copy(row_start_.begin(), row_start_.end(),
                      std::ostream_iterator<int>(os, " "));
            os << "\b]\n";

            os << "Data values = [";
            std::copy(data_.begin(), data_.end(),
                      std::ostream_iterator<T>(os, " "));
            os << "\b]\n";
        }
        const T data(int i)const {
        	return data_[i];
        }

    private:
        std::vector<T> data_;
        // Like in the compressed row sparse matrix format,
        // row_start_.size() is equal to the number of rows + 1.
        std::vector<int> row_start_;

	template <class IntegerIter>
	void setRowStartsFromSizes(IntegerIter rowsize_beg, IntegerIter rowsize_end)
	{
            // Since we do not store the row sizes, but cumulative row sizes,
            // we have to create the cumulative ones.
            int num_rows = rowsize_end - rowsize_beg;
            if (num_rows < 1) {
                THROW("Must have at least one row. Got " << num_rows << " rows.");
            }
#ifndef NDEBUG
            if (*std::min_element(rowsize_beg, rowsize_end) < 0) {
                THROW("All row sizes must be at least 0.");
            }
#endif
            row_start_.resize(num_rows + 1);
            row_start_[0] = 0;
            std::partial_sum(rowsize_beg, rowsize_end, row_start_.begin() + 1);
            // Check that data_ and row_start_ match.
            if (int(data_.size()) != row_start_.back()) {
                THROW("End of row start indices different from data size.");
            }

	}
    };

} // namespace Opm


#endif // OPENRS_SPARSETABLE_HEADER
