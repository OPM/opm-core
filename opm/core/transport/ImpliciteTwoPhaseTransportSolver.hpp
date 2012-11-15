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
#include <vector>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#ifndef IMPLICITETWOPHASETRANSPORTSOLVER_HPP
#define IMPLICITETWOPHASETRANSPORTSOLVER_HPP

// implicite transprot solver
class ImpliciteTwoPhaseTransportSolver : public TwoPhaseTransportSolver
{
public:
    /// Construct solver.
    /// \param[in] grid      A 2d or 3d grid.
    /// \param[in] props     Rock and fluid properties.
    /// \param[in] tol       Tolerance used in the solver.
    /// \param[in] maxit     Maximum number of non-linear iterations used.
    TransportModelTwophase(const UnstructuredGrid& grid,
                           const Opm::IncompPropertiesInterface& props,
                           Opm::parameter::ParameterGroup& param);


    /// Solve for saturation at next timestep.
    /// \param[in] darcyflux         Array of signed face fluxes.
    /// \param[in] porevolume        Array of pore volumes.
    /// \param[in] source            Transport source term.
    /// \param[in] dt                Time step.
    /// \param[in, out] saturation   Phase saturations.
    void solve(const double* darcyflux,
               const double* porevolume,
               const double* source,
               const double dt,
               std::vector<double>& saturation);
private:
    typedef SimpleFluid2pWrappingProps TwophaseFluid;
    typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;

    using namespace Opm::ImplicitTransportDefault;

    typedef NewtonVectorCollection< ::std::vector<double> >      NVecColl;
    typedef JacobianSystem        < struct CSRMatrix, NVecColl > JacSys;

    template <class Vector>
    class MaxNorm {
    public:
        static double
        norm(const Vector& v) {
            return AccumulationNorm <Vector, MaxAbs>::norm(v);
        }
    };

    typedef Opm::ImplicitTransport<TransportModel,
                                   JacSys        ,
                                   MaxNorm       ,
                                   VectorNegater ,
                                   VectorZero    ,
                                   MatrixZero    ,
                                   VectorAssign  > TransportSolver;
    TransportSolver tsolver_,
    UnstructuredGrid* grid_,
    Opm::ImplicitTransportDetails::NRControl ctrl_,
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props_,
    boost::scoped_ptr<Opm::RockCompressibility> rock_comp_(rock_comp),
    boost::scoped_ptr<Opm::WellsManager> wells_(wells),

};

#endif // IMPLICITETWOPHASETRANSPORTSOLVER_HPP
