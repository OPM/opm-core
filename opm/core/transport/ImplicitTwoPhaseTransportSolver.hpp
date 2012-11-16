/*===========================================================================
//
// File: ImpliciteTwoPhaseTransportSolver.hpp
//
// Author: hnil <hnil@sintef.no>
//
// Created: 9 Nov 2012
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

#ifndef IMPLICITETWOPHASETRANSPORTSOLVER_HPP
#define IMPLICITTWOPHASETRANSPORTSOLVER_HPP
#include <vector>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/transport/TwoPhaseTransportSolver.hpp>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/transport/SimpleFluid2pWrappingProps.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>
#include <boost/scoped_ptr.hpp>

#include <opm/core/fluid/RockCompressibility.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/simulator/WellState.hpp>
namespace Opm{
    // implicite transprot solver
    class ImplicitTwoPhaseTransportSolver : public TwoPhaseTransportSolver
    {
    public:
        /// Construct solver.
        /// \param[in] grid      A 2d or 3d grid.
        /// \param[in] props     Rock and fluid properties.
        /// \param[in] tol       Tolerance used in the solver.
        /// \param[in] maxit     Maximum number of non-linear iterations used.
        ImplicitTwoPhaseTransportSolver(
                const Opm::WellsManager& wells,
                const Opm::RockCompressibility& rock_comp,
                const ImplicitTransportDetails::NRControl& ctrl,
                SinglePointUpwindTwoPhase<Opm::SimpleFluid2pWrappingProps>& model,
                const UnstructuredGrid& grid,
                const Opm::IncompPropertiesInterface& props,
                const parameter::ParameterGroup& param);

        ~ImplicitTwoPhaseTransportSolver(){
           destroy_transport_source(tsrc_);
        }

        /// Solve for saturation at next timestep.
        /// \param[in] darcyflux         Array of signed face fluxes.
        /// \param[in] porevolume        Array of pore volumes.
        /// \param[in] source            Transport source term.
        /// \param[in] dt                Time step.
        /// \param[in, out] saturation   Phase saturations.
        void solve(const double* porevolume,
                   const double* source,
                   const double dt,
                   Opm::TwophaseState& state,
                   Opm::WellState& well_state);

    private:
        ImplicitTwoPhaseTransportSolver(const ImplicitTwoPhaseTransportSolver&);
        ImplicitTwoPhaseTransportSolver& operator=(const ImplicitTwoPhaseTransportSolver&);
        typedef Opm::SimpleFluid2pWrappingProps TwophaseFluid;
        typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;

        //using namespace ImplicitTransportDefault;

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
        //Opm::LinearSolverFactory
        Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver linsolver_;
        TransportSolver tsolver_;
        const UnstructuredGrid& grid_;
        const Opm::ImplicitTransportDetails::NRControl& ctrl_;
        const Opm::IncompPropertiesInterface& props_;
        const Opm::RockCompressibility& rock_comp_;
        const Opm::WellsManager& wells_;

        TransportSource* tsrc_;

    };
}
#endif // IMPLICITETWOPHASETRANSPORTSOLVER_HPP
