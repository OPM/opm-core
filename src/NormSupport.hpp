/*===========================================================================
//
// File: NormSupport.hpp
//
// Created: 2011-10-04 19:37:35+0200
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

#ifndef OPM_NORMSUPPORT_HPP_HEADER
#define OPM_NORMSUPPORT_HPP_HEADER

#include <algorithm>
#include <functional>
#include <numeric>

namespace Opm {
    namespace ImplicitTransportDefault {
        template <typename T>
        class MaxAbs : public ::std::binary_function <double, T, double> {
        public:
            double
            operator()(double x, const T& y) {
                return std::max(std::abs(x), std::abs(y));
            }

            static double
            postprocess(double nrm_inf) { return nrm_inf; }
        };

        template <typename T>
        class SumAbs : public ::std::binary_function <double, T, double> {
        public:
            double
            operator()(double x, const T& y) {
                return std::abs(x) + std::abs(y);
            }

            static double
            postprocess(double nrm_1) { return nrm_1; }
        };

        template <typename T>
        class Euclid : public ::std::binary_function <double, T, double> {
        public:
            double
            operator()(double x, const T& y) {
                const double ay = ::std::abs(y);

                return std::abs(x) + ay*ay;
            }

            static double
            postprocess(double nrm2) { return ::std::sqrt(nrm2); }
        };

        template <class Vector, template <typename> class NormImpl>
        class AccumulationNorm {
        public:
            static double
            norm(const Vector& v) {
                typedef typename Vector::value_type VT;

                double nrm = ::std::accumulate(v.begin(), v.end(), VT(0.0),
                                               NormImpl<VT>());

                return NormImpl<VT>::postprocess(nrm);
            }
        };
    }
}

#endif  /* OPM_NORMSUPPORT_HPP_HEADER */
