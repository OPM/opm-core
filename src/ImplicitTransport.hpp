/*===========================================================================
//
// File: ImplicitTransport.hpp
//
// Created: 2011-09-29 10:38:42+0200
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

#ifndef OPM_IMPLICITTRANSPORT_HPP_HEADER
#define OPM_IMPLICITTRANSPORT_HPP_HEADER

#include "ImplicitAssembly.hpp"

namespace Opm {
    namespace ImplicitTransportDetails {
        struct NRControl {
            NRControl()
                : max_it(1),
                  atol(1.0e-6),
                  rtol(5.0e-7),
                  dxtol(1.0e-8)
            {}

            int    max_it;
            double atol  ;
            double rtol  ;
            double dxtol ;
        };

        struct NRReport {
            int    nit;
            int    flag;
            double norm_res;
            double norm_dx;
        };
    }

    template <class Model                  ,
              class JacobianSystem         ,
              template <class> class VNorm ,
              template <class> class VNeg  ,
              template <class> class VZero ,
              template <class> class MZero >
    class ImplicitTransport {
    public:
        ImplicitTransport(Model& model)
            : model_(model),
              asm_  (model)
        {}

        template <class Grid          ,
                  class SourceTerms   ,
                  class ReservoirState,
                  class LinearSolver  >
        void solve(const Grid&                                g       ,
                   const SourceTerms*                         src     ,
                   const double                               dt      ,
                   const ImplicitTransportDetails::NRControl& ctrl    ,
                   ReservoirState&                            state   ,
                   LinearSolver&                              linsolve,
                   ImplicitTransportDetails::NRReport&        rpt     ) {

            typedef typename JacobianSystem::vector_type vector_type;
            typedef typename JacobianSystem::matrix_type matrix_type;

            asm_.createSystem(g, sys_);

            model_.initStep(state, g, sys_);
            model_.initIteration(state, g, sys_);

            MZero<matrix_type>::zero(sys_.writableMatrix());
            VZero<vector_type>::zero(sys_.vector().writableResidual());

            asm_.assemble(state, g, src, dt, sys_);

            const double nrm_res0 =
                VNorm<vector_type>::norm(sys_.vector().residual());

            rpt.norm_res = nrm_res0;
            rpt.norm_dx  = -1.0;
            rpt.nit      = 0;

            bool done = rpt.norm_res < ctrl.atol;

            while (! done) {
                linsolve.solve(sys_.matrix(),
                               sys_.vector().residual(),
                               sys_.vector().writableIncrement());

                VNeg<vector_type>::negate(sys_.vector().writableIncrement());

                model_.finishIteration(state, g, sys_.vector());

                rpt.norm_dx =
                    VNorm<vector_type>::norm(sys_.vector().increment());

                std::cerr << rpt.norm_dx << '\n';

                sys_.vector().addIncrement();
                model_.initIteration(state, g, sys_);

                MZero<matrix_type>::zero(sys_.writableMatrix());
                VZero<vector_type>::zero(sys_.vector().writableResidual());

                asm_.assemble(state, g, src, dt, sys_);
                rpt.norm_res =
                    VNorm<vector_type>::norm(sys_.vector().residual());

                rpt.nit++;

                std::cerr <<  "Iteration " << rpt.nit << " norm :" << rpt.norm_res << '\n';

                done = (rpt.norm_res < ctrl.atol)            ||
                       (rpt.norm_res < ctrl.rtol * nrm_res0) ||
                       (rpt.norm_dx  < ctrl.dxtol)           ||
                       (rpt.nit == ctrl.max_it);
            }

            model_.finishStep(g, sys_.vector().solution(), state);

            if      (rpt.norm_res < ctrl.atol)            { rpt.flag =  1; }
            else if (rpt.norm_res < ctrl.rtol * nrm_res0) { rpt.flag =  2; }
            else if (rpt.norm_dx  < ctrl.dxtol)           { rpt.flag =  3; }
            else                                          { rpt.flag = -1; }
        }

    private:
        ImplicitTransport           (const ImplicitTransport&);
        ImplicitTransport& operator=(const ImplicitTransport&);

#if 0
        using Model::initStep;
        using Model::initIteration;
        using Model::finishIteration;
        using Model::finishStep;
#endif

        Model&                  model_;
        ImplicitAssembly<Model> asm_;
        JacobianSystem          sys_;
    };
}
#endif  /* OPM_IMPLICITTRANSPORT_HPP_HEADER */
