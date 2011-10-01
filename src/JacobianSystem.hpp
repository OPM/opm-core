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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>

namespace Opm {
    namespace ImplicitTransportDefault {
        template <typename T>
        class MaxAbs : public std::binary_function<T, T, T> {
        public:
            T operator()(const T& x, const T& y) {
                return std::max(std::abs(x), std::abs(y));
            }
        };

        template <class Vector>
        class MaxNorm {
        public:
            typename Vector::value_type
            norm(const Vector& v) {
                typedef typename Vector::value_type VT;

                return std::accumulate(v.begin(), v.end(),
                                       VT(0), MaxAbs<VT>());
            }
        };

        template <typename T>
        class SumAbs : public std::binary_function<T, T, T> {
        public:
            T operator()(const T& x, const T& y) {
                return std::abs(x) + std::abs(y);
            }
        };

        template <class Vector>
        class TaxiCabNorm {
        public:
            typename Vector::value_type
            norm(const Vector& v) {
                typedef typename Vector::value_type VT;

                return std::accumulate(v.begin(), v.end(),
                                       VT(0), SumAbs<VT>());
            }
        };

        template <class Vector>
        class EuclidianNorm {
        public:
            typename Vector::value_type
            norm(const Vector& v) {
                typedef typename Vector::value_type VT;
                typedef typename Vector::iterator   VI;

                VT ret2 = 0;
                for (VI i = v.begin(), e = v.end(); i != e; ++i) {
                    VT  x = std::abs(*i);

                    ret2 += x * x;
                }

                return std::sqrt(ret2);
            }
        };

        template <class Vector,
                  template<class> class Norm = MaxNorm>
        class DefaultNewtonVector {
            typedef typename Vector::value_type VT;
        public:
            void
            resize(size_t m) { v_.resize(m); }

            void
            negate() {
                std::transform(v_.begin(),
                               v_.end  (),
                               v_.begin(),
                               std::negate<VT>());
            }

            typename Vector::value_type
            norm() const { return Norm<Vector>::norm(v_); }

            typename Vector::value_type&
            operator[](size_t i)       { return v_[i]; }

            const typename Vector::value_type&
            operator[](size_t i) const { return v_[i]; }

        private:
            Vector v_;
        };

        template <class Vector>
        class NewtonVectors {
            typedef typename std::array<Vector, 3> VC ;
            typedef typename VC ::iterator         VCI;

            enum { Residual = 0, Increment = 1, Solution = 2 };

        public:
            void
            setSize(size_t m) {
                for (VCI i = v_.begin(), e = v_.end(); i != e; ++i) {
                    i->resize(m);
                }
            }

            template <class Block>
            void
            assembleBlock(size_t n, size_t i, const Block& b) {
                for (size_t k = 0; k < n; ++k) {
                    residual()[i*n + k] += b[k];
                }
            }

            void
            addIncrement() {
                std::transform(v_[ Solution  ].begin(),
                               v_[ Solution  ].end  (),
                               v_[ Increment ].begin(),
                               v_[ Solution  ].begin(),
                               std::plus<typename Vector::value_type>());
            }

            Vector&       residual ()       { return v_[ Residual  ]; }
            const Vector& residual () const { return v_[ Residual  ]; }

            Vector&       increment()       { return v_[ Increment ]; }
            const Vector& increment() const { return v_[ Increment ]; }

            Vector&       solution ()       { return v_[ Solution  ]; }
            const Vector& solution () const { return v_[ Solution  ]; }

        private:
            VC v_;
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
            setSize(size_t m, size_t n, size_t nnz = 0);
        } */;

        template <class Matrix, class Vector>
        class JacobianSystem {
        public:
            typedef Matrix matrix_type;

            MatrixBlockAssembler<Matrix>& matrix() { return mba_ ; }
            NewtonVectors       <Vector>& vector() { return vecs_; }

        private:
            MatrixBlockAssembler<Matrix> mba_ ;
            NewtonVectors       <Vector> vecs_;
        };
    }
}

#endif  /* OPM_JACOBIANSYSTEM_HPP_HEADER */
