/*===========================================================================
//
// File: JacobianSystem.hpp
//
// Created: 2011-09-30 19:23:31+0200
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

#ifndef OPM_JACOBIANSYSTEM_HPP_HEADER
#define OPM_JACOBIANSYSTEM_HPP_HEADER

#include <cassert>
#include <cmath>
#include <cstddef>

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>

namespace Opm {
    namespace ImplicitTransportDefault {
        template <class BaseVec>
        class VectorAdder {
        public:
            // y += x
            static void
            add(const BaseVec& x, BaseVec& y) {
                typedef typename BaseVec::value_type VT;

                ::std::transform(x.begin(), x.end(),
                                 y.begin(),
                                 y.begin(),
                                 ::std::plus<VT>());
            }
        };

        template <class BaseVec>
        class VectorNegater {
        public:
            // x *= -1
            static void
            negate(BaseVec& x) {
                typedef typename BaseVec::value_type VT;

                ::std::transform(x.begin(), x.end(),
                                 x.begin(),
                                 ::std::negate<VT>());
            }
        };

        template <class BaseVec>
        class VectorZero {
            static void
            zero(BaseVec& x) {
                typedef typename BaseVec::value_type VT;

                ::std::fill(x.begin(), x.end(), VT(0.0));
            }
        };

        template <class BaseVec>
        class VectorBlockAssembler {
        public:
            template <class Block>
            static void
            assemble(::std::size_t ndof,
                     ::std::size_t i   ,
                     const Block&  b   ,
                     BaseVec&      vec ) {

                for (::std::size_t d = 0; d < ndof; ++d) {
                    vec[i*ndof + d] += b[d];
                }
            }
        };

        template <class BaseVec>
        class VectorSizeSetter {
        public:
            VectorSizeSetter(BaseVec& v) : v_(v) {}

            void
            setSize(::std::size_t ndof, ::std::size_t m) {
                v_.resize(ndof * m);
            }

        private:
            BaseVec& v_;
        };

        template <class                  BaseVec                         ,
                  template <class> class VSzSetter = VectorSizeSetter    ,
                  template <class> class VAdd      = VectorAdder         ,
                  template <class> class VBlkAsm   = VectorBlockAssembler>
        class NewtonVectorCollection {
            enum { Residual = 0, Increment = 1, Solution = 2 };

        public:
            void
            setSize(::std::size_t ndof, ::std::size_t m) {
                typedef typename ::std::array<BaseVec, 3>::iterator VAI;

                for (VAI i = vcoll_.begin(), e = vcoll_.end(); i != e; ++i) {
                    VSzSetter<BaseVec>(*i).setSize(ndof, m);
                }

                ndof_ = ndof;
            }

            void
            addIncrement() {
                VAdd<BaseVec>::add(vcoll_[ Increment ], vcoll_[ Solution ]);
            }

            template <class Block>
            void
            assembleBlock(::std::size_t ndof,
                          ::std::size_t i   ,
                          const Block&  b   ) {

                assert (ndof_ >  0    );
                assert (ndof  == ndof_);

                VBlkAsm<BaseVec>::assemble(ndof, i, b, vcoll_[ Residual ]);
            }

            typedef BaseVec vector_type;

            const vector_type& increment() const { return vcoll_[ Increment ]; }
            const vector_type& residual () const { return vcoll_[ Residual  ]; }
            const vector_type& solution () const { return vcoll_[ Solution  ]; }

            // Write access for Newton solver purposes
            vector_type& writableIncrement()     { return vcoll_[ Increment ]; }
            vector_type& writableResidual ()     { return vcoll_[ Residual  ]; }
            vector_type& writableSolution ()     { return vcoll_[ Solution  ]; }

        private:
            ::std::size_t            ndof_ ;
            ::std::array<BaseVec, 3> vcoll_;
        };

        template <class Matrix>
        class MatrixBlockAssembler
        /* {
        public:
            template <class Block>
            void
            assembleBlock(size_t n, size_t i, size j, const Block& b);

            template <class Connections>
            void
            createBlockRow(size_t i, const Connections& conn, size_t n);

            void
            finalizeStructure();

            void
            setSize(size_t ndof, size_t m, size_t n, size_t nnz = 0);

            const Matrix&
            matrix();
        } */;

        template <class Matrix        ,
                  class NVecCollection>
        class JacobianSystem {
        public:
            JacobianSystem() {}

            typedef Matrix                               matrix_type;
            typedef MatrixBlockAssembler<Matrix>         assembler_type;
            typedef typename NVecCollection::vector_type vector_type;

            assembler_type&    matasm()       { return mba_         ; }
            NVecCollection&    vector()       { return sysvec_      ; }
            const matrix_type& matrix() const { return mba_.matrix(); }

            void
            setSize(::std::size_t ndof,
                    ::std::size_t m   ,
                    ::std::size_t nnz = 0) {

                mba_   .setSize(ndof, m, m, nnz);
                sysvec_.setSize(ndof, m        );
            }

        private:
            JacobianSystem           (const JacobianSystem&);
            JacobianSystem& operator=(const JacobianSystem&);

            MatrixBlockAssembler<Matrix> mba_   ; // Coefficient matrix
            NVecCollection               sysvec_; // Residual
        };
    }
}

#endif  /* OPM_JACOBIANSYSTEM_HPP_HEADER */
