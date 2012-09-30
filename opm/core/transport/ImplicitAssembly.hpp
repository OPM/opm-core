/*===========================================================================
//
// File: ImplicitAssembly.hpp
//
// Created: 2011-09-28 10:00:46+0200
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

#ifndef OPM_IMPLICITASSEMBLY_HPP_HEADER
#define OPM_IMPLICITASSEMBLY_HPP_HEADER

#include <cstddef>
#include <algorithm>
#include <vector>

namespace Opm {
    template <class Model>
    class ImplicitAssembly {
        enum { DofPerCell = Model::DofPerCell };

    public:
        ImplicitAssembly(Model& model)
            : model_(model),
              nconn_(-1)   ,
              asm_buffer_()
        {}

        template <class Grid          ,
                  class JacobianSystem>
        void
        createSystem(const Grid&     g  ,//const SourceTerms* tsrc,
                     JacobianSystem& sys) const {

            ::std::size_t m   = g.number_of_cells;
            ::std::size_t nnz = g.number_of_cells
            		          + countConnections(g);
//            		          + tsrc->pf;

            sys.setSize(DofPerCell, m, nnz);

            for (int c = 0; c < g.number_of_cells; ++c) {
                this->createRowStructure(g, c, sys);
            }
            /* is done by modifiying grid
            // add elements for periodic boundary
            for (int pf = 0; pf<tsrc->pf; ++pf){
            	int nconn = 1;
            	::std::vector<int> connections(1);
            	connections[0]=(tsrc->periodic_cells[2*pf]);
            	sys.matasm().createBlockRow(tsrc->periodic_cells[2*pf+1], connections, DofPerCell);
            	connections[0]=(tsrc->periodic_cells[2*pf+1]);
            	sys.matasm().createBlockRow(tsrc->periodic_cells[2*pf], connections, DofPerCell);
            }
            */
            sys.matasm().finalizeStructure();
        }

        template <class ReservoirState,
                  class Grid          ,
                  class SourceTerms   ,
                  class JacobianSystem>
        void
        assemble(const ReservoirState& state,
                 const Grid&           g    ,
                 const SourceTerms*    src  ,
                 const double          dt   ,
                 JacobianSystem&       sys  ) {

            for (int c = 0; c < g.number_of_cells; ++c) {
                this->computeCellContrib(state, g, dt, c);
                this->assembleCellContrib(g, c, sys);
            }

            if (src != 0) {
                this->assembleSourceContrib(g, src, dt, sys);
            }
        }

    private:
        template <class Grid>
        int
        countConnections(const Grid& g, int c) const {
            int i, f, c1, c2, n;

            for (i = g.cell_facepos[c + 0], n = 0;
                 i < g.cell_facepos[c + 1]; ++i) {
                f  = g.cell_faces[i];
                c1 = g.face_cells[2*f + 0];
                c2 = g.face_cells[2*f + 1];

                n += (c1 >= 0) && (c2 >= 0);
            }

            return n;
        }

        template <class Grid>
        int
        countConnections(const Grid& g) const {
            int n = 0;

            for (int c = 0; c < g.number_of_cells; ++c) {
                n += countConnections(g, c);
            }

            return n;
        }

        template <class Grid, class System>
        void
        createRowStructure(const Grid& g  ,
                           const int   c  ,
                           System&     sys) const {

            int nconn = countConnections(g, c);

            ::std::vector<int> connections;
            connections.reserve  (nconn + 1);
            connections.push_back(c);

            for (int i = g.cell_facepos[c + 0];
                 i     < g.cell_facepos[c + 1]; ++i) {
                int f  = g.cell_faces[i];
                int c1 = g.face_cells[2*f + 0];
                int c2 = g.face_cells[2*f + 1];

                if ((c1 >= 0) && (c2 >= 0)) {
                    connections.push_back((c1 == c) ? c2 : c1);
                }
            }
            sys.matasm().createBlockRow(c, connections, DofPerCell);

        }

        template <class ReservoirState, class Grid>
        int
        computeCellContrib(const ReservoirState& state,
                           const Grid&           g    ,
                           const double          dt   ,
                           const int             c    ) {
            const int ndof  = DofPerCell;
            const int ndof2 = ndof * ndof;

            nconn_ = countConnections(g, c);

            asm_buffer_.resize((2*nconn_ + 1)*ndof2 + (nconn_ + 2)*ndof);
            std::fill(asm_buffer_.begin(), asm_buffer_.end(), 0.0);

            double* F  = &asm_buffer_[(2*nconn_ + 1) * ndof2];
            double* J1 = &asm_buffer_[(0*nconn_ + 1) * ndof2];
            double* J2 = J1         + (1*nconn_ + 0) * ndof2 ;

            model_.initResidual(c, F);
            F += ndof;

            for (int i = g.cell_facepos[c + 0];
                 i     < g.cell_facepos[c + 1]; ++i) {
                int f  = g.cell_faces[i];
                int c1 = g.face_cells[2*f + 0];
                int c2 = g.face_cells[2*f + 1];

                if ((c1 >= 0) && (c2 >= 0)) {
                    model_.fluxConnection(state, g, dt, c, f, J1, J2, F);
                    J1 += ndof2;  J2 += ndof2;   F += ndof;
                }
            }

            model_.accumulation(g, c, &asm_buffer_[0], F);

            return 1;
        }

        template <class Grid, class System>
        void
        assembleCellContrib(const Grid& g  ,
                            const int   c  ,
                            System&     sys) const {
            const int ndof  = DofPerCell;
            const int ndof2 = ndof * ndof;

            typedef std::vector<int>::size_type sz_t;

            const double* J1 = &asm_buffer_[0];
            const double* J2 = J1 + ((1*nconn_ + 1) * ndof2);

            // Assemble contributions from accumulation term
            sys.matasm().assembleBlock(ndof, c, c, J1);  J1 += ndof2;

            // Assemble connection contributions.
            for (int i = g.cell_facepos[c + 0];
                 i     < g.cell_facepos[c + 1]; ++i) {
                int f  = g.cell_faces[i];
                int c1 = g.face_cells[2*f + 0];
                int c2 = g.face_cells[2*f + 1];

                c2 = (c1 == c) ? c2 : c1;

                if (c2 >= 0) {
                    sys.matasm().assembleBlock(ndof, c, c , J1);
                    sys.matasm().assembleBlock(ndof, c, c2, J2);

                    J1 += ndof2;
                    J2 += ndof2;
                }
            }
            // Assemble residual
            const double* F = &asm_buffer_[(2*nconn_ + 1) * ndof2];
            for (int conn = 0; conn < nconn_ + 2; ++conn, F += ndof) {
                sys.vector().assembleBlock(ndof, c, F);
            }
        }

        template <class Grid, class SourceTerms, class System>
        void
        assembleSourceContrib(const Grid&        g,
                              const SourceTerms* src,
                              const double       dt,
                              System&            sys) {
            const int ndof  = DofPerCell;
            const int ndof2 = ndof * ndof;

            for (int i = 0; i < src->nsrc; ++i) {
                std::fill_n(asm_buffer_.begin(), ndof2 + ndof, 0.0);

                double *J = &asm_buffer_[0];
                double *F = J + ndof2;

                model_.sourceTerms(g, src, i, dt, J, F);

                const int c = src->cell[i];

                sys.matasm().assembleBlock(ndof, c, c, J);
                sys.vector().assembleBlock(ndof, c,    F);
            }
        }

        Model&              model_     ;
        int                 nconn_     ;
        std::vector<double> asm_buffer_;
    };
}
#endif  /* OPM_IMPLICITASSEMBLY_HPP_HEADER */
