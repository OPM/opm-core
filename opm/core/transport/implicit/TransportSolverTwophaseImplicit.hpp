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

#ifndef OPM_TRANSPORTSOLVERTWOPHASEIMPLICIT_HEADER_INCLUDED
#define OPM_TRANSPORTSOLVERTWOPHASEIMPLICIT_HEADER_INCLUDED

#include <opm/core/transport/TransportSolverTwophaseInterface.hpp>
#include <opm/core/transport/implicit/SimpleFluid2pWrappingProps.hpp>
#include <opm/core/transport/implicit/SinglePointUpwindTwoPhase.hpp>
#include <opm/core/transport/implicit/ImplicitTransport.hpp>
#include <opm/core/transport/implicit/transport_source.h>
#include <opm/core/transport/implicit/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/implicit/NormSupport.hpp>
#include <opm/core/transport/implicit/ImplicitAssembly.hpp>
#include <opm/core/transport/implicit/ImplicitTransport.hpp>
#include <opm/core/transport/implicit/JacobianSystem.hpp>
#include <opm/core/transport/implicit/CSRMatrixBlockAssembler.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <memory>
#include <vector>



namespace Opm
{

    // Implicit transport solver using Newton-Raphson.
    class TransportSolverTwophaseImplicit : public TransportSolverTwophaseInterface
    {
    public:
        /// Construct solver.
        /// \param[in] grid      A 2d or 3d grid.
        /// \param[in] props     Rock and fluid properties.
        /// \param[in] porevol   Pore volumes
        /// \param[in] gravity   Gravity vector (null for no gravity).
        /// \param[in] half_trans Half-transmissibilities (one-sided)
        /// \param[in] maxit     Maximum number of non-linear iterations used.
        TransportSolverTwophaseImplicit(const UnstructuredGrid& grid,
                                        const Opm::IncompPropertiesInterface& props,
                                        const std::vector<double>& porevol,
                                        const double* gravity,
                                        const std::vector<double>& half_trans,
                                        const parameter::ParameterGroup& param);

        virtual ~TransportSolverTwophaseImplicit();

        /// Solve for saturation at next timestep.
        /// \param[in]      porevolume   Array of pore volumes.
        /// \param[in]      source       Transport source term. For interpretation see Opm::computeTransportSource().
        /// \param[in]      dt           Time step.
        /// \param[in, out] state        Reservoir state. Calling solve() will read state.faceflux() and
        ///                              read and write state.saturation().
        virtual void solve(const double* porevolume,
                           const double* source,
                           const double dt,
                           TwophaseState& state);

    private:
        // Disallow copying and assignment.
        TransportSolverTwophaseImplicit(const TransportSolverTwophaseImplicit&);
        TransportSolverTwophaseImplicit& operator=(const TransportSolverTwophaseImplicit&);

        // Defining types for the underlying transport solver.
        typedef Opm::SimpleFluid2pWrappingProps TwophaseFluid;
        typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;
        typedef ImplicitTransportDefault::NewtonVectorCollection< ::std::vector<double> >      NVecColl;
        typedef ImplicitTransportDefault::JacobianSystem        < struct CSRMatrix, NVecColl > JacSys;
        template <class Vector>
        class MaxNorm {
        public:
            static double
            norm(const Vector& v) {
                return ImplicitTransportDefault::AccumulationNorm <Vector, ImplicitTransportDefault::MaxAbs>::norm(v);
            }
        };
        typedef Opm::ImplicitTransport<TransportModel,
                                       JacSys        ,
                                       MaxNorm       ,
                                       ImplicitTransportDefault::VectorNegater ,
                                       ImplicitTransportDefault::VectorZero    ,
                                       ImplicitTransportDefault::MatrixZero    ,
                                       ImplicitTransportDefault::VectorAssign  > TransportSolver;

        // Data members.
        Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver linsolver_;
        Opm::SimpleFluid2pWrappingProps fluid_;
        SinglePointUpwindTwoPhase<Opm::SimpleFluid2pWrappingProps> model_;
        TransportSolver tsolver_;
        const UnstructuredGrid& grid_;
        Opm::ImplicitTransportDetails::NRControl ctrl_;
        const Opm::IncompPropertiesInterface& props_;
        TransportSource* tsrc_;
        double initial_porevolume_cell0_;
    };

} // namespace Opm

#endif // OPM_TRANSPORTSOLVERTWOPHASEIMPLICIT_HEADER_INCLUDED
