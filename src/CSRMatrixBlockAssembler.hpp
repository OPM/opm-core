/*===========================================================================
//
// File: CSRMatrixBlockAssembler.hpp
//
// Created: 2011-10-03 12:40:56+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_CSRMATRIXBLOCKASSEMBLER_HPP_HEADER
#define OPM_CSRMATRIXBLOCKASSEMBLER_HPP_HEADER

#include <cstddef>
#include <cstdlib>

#include <algorithm>
#include <vector>

#include <sparse_sys.h>

#include "JacobianSystem.hpp"

namespace Opm {
    namespace ImplicitTransportDefault {

        template <>
        class MatrixZero <struct CSRMatrix> {
        public:
            static void
            zero(struct CSRMatrix& A) {
                csrmatrix_zero(&A);
            }
        };

        template <>
        class MatrixBlockAssembler<struct CSRMatrix> {
        public:
            template <class Block>
            void
            assembleBlock(::std::size_t ndof,
                          ::std::size_t i   ,
                          ::std::size_t j   ,
                          const Block&  b   ) {

                assert (ndof >  0);
                assert (ndof == ndof_);

                const ::std::size_t start = ia_[i*ndof + 0];
                const ::std::size_t off   =
                    csrmatrix_elm_index(i * ndof, j * ndof, &mat_) - start;

                for (::std::size_t row = 0; row < ndof; ++row) {
                    const ::std::size_t J = ia_[i*ndof + row] + off;

                    for (::std::size_t col = 0; col < ndof; ++col) {
                        sa_[J + col] += b[col*ndof + row];
                    }
                }
            }

            template <class Connections>
            void
            createBlockRow(::std::size_t      i  ,
                           const Connections& conn,
                           ::std::size_t      ndof) {

                assert (ndof >  0);
                assert (ndof == ndof_);
                assert (i    == (ia_.size() - 1) / ndof_);  (void) i;

                expandSortConn(conn, ndof);
                const int nconn = static_cast<int>(esconn_.size());

                for (::std::size_t dof = 0; dof < ndof; ++dof) {
                    ja_.insert(ja_.end(), esconn_.begin(), esconn_.end());
                    ia_.push_back(ia_.back() + nconn);
                }

                sa_.insert(sa_.end(), nconn * ndof, double(0.0));

                construct();
                setCSRSize();
            }

            void
            finalizeStructure() {
                setCSRSize();
            }

            void
            setSize(size_t ndof, size_t m, size_t n, size_t nnz = 0) {
                (void) n;

                clear();

                allocate(ndof, m, nnz);

                ia_.push_back(0);
                ndof_ = ndof;

                construct();
            }

            struct CSRMatrix&       matrix()       { return mat_; }
            const struct CSRMatrix& matrix() const { return mat_; }

        private:
            void
            allocate(::std::size_t ndof, ::std::size_t m, ::std::size_t nnz) {
                ia_.reserve(1 + ( m  * ndof));
                ja_.reserve(0 + (nnz * ndof));
                sa_.reserve(0 + (nnz * ndof));
            }

            void
            clear() {
                ia_.resize(0);
                ja_.resize(0);
                sa_.resize(0);
            }

            void
            construct() {
                mat_.ia = &ia_[0];
                mat_.ja = &ja_[0];
                mat_.sa = &sa_[0];
            }

            template <class Connections>
            void
            expandSortConn(const Connections& conn, ::std::size_t ndof) {
                sconn_.resize(0);
                sconn_.reserve(conn.size());

                for (typename Connections::const_iterator
                         c = conn.begin(), e = conn.end(); c != e; ++c) {
                    sconn_.push_back(static_cast<int>(*c));
                }

                ::std::sort(sconn_.begin(), sconn_.end());

                esconn_.resize(0);
                esconn_.reserve(ndof * sconn_.size());

                for (::std::vector<int>::iterator
                         c = sconn_.begin(), e = sconn_.end(); c != e; ++c) {
                    for (::std::size_t dof = 0; dof < ndof; ++dof) {
                        esconn_.push_back(static_cast<int>((*c)*ndof + dof));
                    }
                }
            }

            void
            setCSRSize() {
                mat_.m   = ia_.size() - 1;
                mat_.nnz = ja_.size()    ;
            }

            ::std::size_t         ndof_;

            ::std::vector<int>    ia_;
            ::std::vector<int>    ja_;
            ::std::vector<double> sa_;

            ::std::vector<int>    sconn_ ;
            ::std::vector<int>    esconn_;

            struct CSRMatrix      mat_;
        };
    }
}

#endif  /* OPM_CSRMATRIXBLOCKASSEMBLER_HPP_HEADER */
