/*===========================================================================
//
// File: ImplicitTransport.hpp
//
// Created: 2011-09-29 10:38:42+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
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

    template <class Model         ,
              class JacobianSystem>
    class ImplicitTransport : private Model {
    public:
        template <class Grid          ,
                  class SourceTerms   ,
                  class ReservoirState>
        void solve(const Grid&                                g    ,
                   const SourceTerms&                         src  ,
                   const double                               dt   ,
                   const ImplicitTransportDetails::NRControl& ctrl ,
                   ReservoirState&                            state,
                   ImplicitTransportDetails::NRReport&        rpt  ) {

            asm_.createSystem(g, sys_);

            this->initStep(state, g, sys_);
            this->initIteration(state, g, sys_);

            asm_.assemble(state, g, src, dt, sys_);

            const double nrm_res0 = sys_.vector().residual().norm();
            rpt.norm_res          = nrm_res0;
            rpt.norm_dx           = -1.0;
            rpt.nit               = 0;

            bool done = rpt.norm_res < ctrl.atol;

            while (! done) {
                sys_.solve ();
                sys_.vector().increment().negate();

                this->finishIteration(state, g, sys_.vector());

                nrm_dx = sys_.vector().increment().norm();

                sys_.vector().addIncrement();
                this->initIteration(state, g, sys_);

                asm_.assemble(state, g, src, dt, sys_);
                rpt.norm_res = sys_.vector().residual().norm();

                rpt.nit++;

                done = (rpt.norm_res < ctrl.atol)            ||
                       (rpt.norm_res < ctrl.rtol * nrm_res0) ||
                       (rpt.nit == ctrl.max_it);
            }

            this->finisStep(g, sys_.vector().solution(), state);

            if      (rpt.norm_res < ctrl.atol)            { rpt.flag =  1; }
            else if (rpt.norm_res < ctrl.rtol * nrm_res0) { rpt.flag =  2; }
            else                                          { rpt.flag = -1; }
        }

    private:
        using Model::initStep;
        using Model::initIteration;
        using Model::finishIteration;
        using Model::finishStep;

        ImplicitAssembly<Model> asm_;
        JacobianSystem          sys_;
    };
}
#endif  /* OPM_IMPLICITTRANSPORT_HPP_HEADER */
