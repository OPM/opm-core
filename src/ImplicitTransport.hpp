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
#include <boost/lambda/lambda.hpp>
namespace Opm {
    namespace ImplicitTransportDetails {
        struct NRControl {
            NRControl()
                : max_it(1),
                  atol(1.0e-6),
                  rtol(5.0e-7),
                  dxtol(1.0e-8),
                  max_it_ls(5)
            {}

            int    max_it;
            double atol  ;
            double rtol  ;
            double dxtol ;
            int    max_it_ls;
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
            typedef typename JacobianSystem::vector_collection_type vector_collection_type;
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

            bool done = false;//rpt.norm_res < ctrl.atol;

            while (! done) {
                linsolve.solve(sys_.matrix(),
                               sys_.vector().residual(),
                               sys_.vector().writableIncrement());
                std::ofstream myfile2;
                myfile2.open ("jacobi.txt");
                Dune::printSparseMatrix(myfile2, sys_.matrix(), "", "", 26, 15);
                myfile2.close();

                Dune::writeMatrixToMatlab(sys_.matrix(),"matrix_matlab"); // write jacobi
                std::ofstream myfile;
                myfile.open ("residual.txt");
                Dune::printvector(myfile,sys_.vector().residual(),"redidual","", 1, 10, 15); //write residual
                myfile.close();

                VNeg<vector_type>::negate(sys_.vector().writableIncrement());

                //model_.finishIteration(state, g, sys_.vector());

                rpt.norm_dx =
                    VNorm<vector_type>::norm(sys_.vector().increment());

                int lin_it=0;
                double residual=VNorm<vector_type>::norm(sys_.vector().residual());
                bool finnished=rpt.norm_res<ctrl.atol;//residual < rpt.norm_res;
                double alpha=2.0;
                // store old solution and increasement before line search
                vector_type dx_old(sys_.vector().increment());
                vector_type x_old(sys_.vector().solution());
                while(! finnished){
                	alpha/=2.0;
                	sys_.vector().writableIncrement()=dx_old;
                	sys_.vector().writableIncrement()*=alpha;
                	sys_.vector().writableSolution()=x_old;
                	/*
                	 	should be used if vector_type is std::vector<double>
                		std::vector<double> operator*=(std::vector<double>& vec,const double a){
                	    for_each(vec.begin(),vec.end(), boost::lambda::_1*=a );
                	    return vec;
                	*/
                	sys_.vector().addIncrement();
                    model_.initIteration(state, g, sys_);
                    MZero<matrix_type>::zero(sys_.writableMatrix());
                    VZero<vector_type>::zero(sys_.vector().writableResidual());
                    asm_.assemble(state, g, src, dt, sys_);
                	residual=VNorm<vector_type>::norm(sys_.vector().residual());
                	lin_it +=1;
                	finnished=(residual < rpt.norm_res) || (lin_it> ctrl.max_it_ls);
                	//std::cerr <<  "Line search iteration " << std::scientific  << lin_it << " norm :" << residual <<  " alpha " << alpha << '\n';
                }
                rpt.norm_res =
                    VNorm<vector_type>::norm(sys_.vector().residual());

                rpt.nit++;

                std::cout <<  "Iteration " << std::scientific  << rpt.nit
                		<< " norm :" << rpt.norm_res <<  " alpha " << alpha << std::endl;

                done = (rpt.norm_res < ctrl.atol)            ||
                       (rpt.norm_res < ctrl.rtol * nrm_res0) ||
                       (rpt.norm_dx  < ctrl.dxtol)           ||
                       (lin_it       > ctrl.max_it_ls)       ||
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
