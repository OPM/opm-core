/*===========================================================================
//
// File: SinglePointUpwindTwoPhase.hpp
//
// Created: 2011-09-28 14:21:34+0200
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

#ifndef OPM_SINGLEPOINTUPWINDTWOPHASE_HPP_HEADER
#define OPM_SINGLEPOINTUPWINDTWOPHASE_HPP_HEADER

#include <cassert>

#include <algorithm>
#include <array>
#include <vector>

namespace Opm {
    namespace spu_2p {
        class ModelParameterStorage {
        public:
            ModelParameterStorage(int nc, int totconn)
                : drho_(0.0), mob_(0), dmob_(0),
                  porevol_(0), dg_(0), ds_(0),
                  data_()
            {
                size_t alloc_sz;

                alloc_sz  = 2 * nc;      // mob_
                alloc_sz += 2 * nc;      // dmob_
                alloc_sz += 1 * nc;      // porevol_
                alloc_sz += 1 * totconn; // dg_
                alloc_sz += 1 * nc;      // ds_

                data_.resize(alloc_sz);

                mob_     = &data_[0];
                dmob_    = mob_     + (2 * nc     );
                porevol_ = dmob_    + (2 * nc     );
                dg_      = porevol_ + (1 * nc     );
                ds_      = dg_      + (1 * totconn);
            }

            double&       drho   ()            { return drho_            ; }
            double        drho   ()      const { return drho_            ; }

            double*       mob    (int c)       { return mob_  + (2*c + 0); }
            const double* mob    (int c) const { return mob_  + (2*c + 0); }

            double*       dmob   (int c)       { return dmob_ + (2*c + 0); }
            const double* dmob   (int c) const { return dmob_ + (2*c + 0); }

            double*       porevol()            { return porevol_         ; }
            double        porevol(int c) const { return porevol_[c]      ; }

            double&       dg(int i)            { return dg_[i]           ; }
            double        dg(int i)      const { return dg_[i]           ; }

            double&       ds(int c)            { return ds_[c]           ; }
            double        ds(int c)      const { return ds_[c]           ; }

        private:
            double  drho_   ;
            double *mob_    ;
            double *dmob_   ;
            double *porevol_;
            double *dg_     ;
            double *ds_     ;

            std::vector<double> data_;
        };
    }


    template <class TwophaseFluid>
    class SinglePointUpwindTwoPhase : private TwophaseFluid {
    public:
        template <class Grid>
        SinglePointUpwindTwoPhase(const Grid&                g        ,
                                  const std::vector<double>& porevol  ,
                                  const double*              grav  = 0,
                                  const double*              trans = 0)
            : TwophaseFluid()                                   ,
              gravity_     ((grav != 0) && (trans != 0))        ,
              f2hf_        (2 * g.number_of_faces, -1)          ,
              store_       (g.number_of_cells,
                            g.cell_facepos[ g.number_of_cells ]),
        {
            if (gravity_) {
                store_.drho() = this->density(0) - this->density(1);

                this->computeStaticGravity(g, grav, trans);
            }

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                for (; i < g.cell_facepos[c + 1]; ++i) {
                    const int f = g.cell_faces[i];
                    const int p = 1 - (g.face_cells[2*f + 0] == c);

                    f2hf_[2*f + p] = i;
                }
            }

            std::copy(porevol.begin(), porevol.end(), data_.porevol());
        }

        // -----------------------------------------------------------------
        // System assembly innards
        // -----------------------------------------------------------------
        void
        initResidual(const int c, double* F) const {
            (void) c;       // Suppress 'unused' warning
            *F = 0.0;
        }

        template <class ReservoirState,
                  class Grid          >
        void
        fluxConnection(const ReservoirState& state,
                       const Grid&           g    ,
                       const double          dt   ,
                       const int             c    ,
                       const int             f    ,
                       double*               J1   ,
                       double*               J2   ,
                       double*               F    ) const {

            const int *n = g.face_cells + (2 * f);
            double dflux = state.faceflux[f];
            double gflux = gravityFlux(f);

            int    pix[2];
            double m[2], dm[2];
            upwindMobility(dflux, gflux, n, pix, m, dm);

            assert (! ((m[0] < 0) || (m[1] < 0)));

            double mt = m[0] + m[1];
            assert (mt > 0);

            double sgn  = 2.0*(c1 == c) - 1.0;
            dflux      *= sgn;
            gflux      *= sgn;

            double       f1 = m[0] / mt;
            const double v1 = dflux + m[1]*gflux;

            // Assemble residual contributions
            *F += dt * f1 * v1;

            // Assemble Jacobian (J1 <-> c, J2 <-> other)
            double *J[2];
            if (n[0] == c) { J[0] = J1; J[1] = J2; }
            else           { J[0] = J2; J[1] = J1; }

            // dF/dm_1 \cdot dm_1/ds
            *J[ pix[0] ] += dt * (1 - f1) / mt * v1    * dm[0];

            /* dF/dm_2 \cdot dm_2/ds */
            *J[ pix[1] ] -= dt * f1       / mt * v1    * dm[1];
            *J[ pix[1] ] += dt * f1            * gflux * dm[1];
        }

        template <class Grid>
        void
        accumulation(const Grid& g,
                     const int   c,
                     double*     J,
                     double*     F) const {
            (void) g;

            const double pv = data_.porevol(c);

            *J += pv;
            *F += pv * data_.ds(c);
        }

        template <class Grid       ,
                  class SourceTerms>
        void
        sourceTerms(const Grid&        g  ,
                    const SourceTerms& src,
                    const int          i  ,
                    const double       dt ,
                    double*            J  ,
                    double*            F  ) const {

            double dflux = -src.flux[i]; // .flux[] is rate of *inflow*

            if (dflux < 0) {
                // src -> cell, affects residual only.
                *F += dt * dflux * src.saturation[2*i + 0];
            } else {
                // cell -> src
                const int     c  = src.cell[i];
                const double* m  = data_.mob (c);
                const double* dm = data_.dmob(c);

                const double  mt = m[0] + m[1];

                assert (! ((m[0] < 0) || (m[1] < 0)));
                assert (mt > 0);

                const double f  = m[0] / mt;
                const double df = ((1 - f)*dm[0] - f*dm[1]) / mt;

                *F += dt * dflux *  f;
                *J += dt * dflux * df;
            }
        }

        // -----------------------------------------------------------------
        // Newton control
        // -----------------------------------------------------------------
        template <class ReservoirState,
                  class Grid          ,
                  class JacobianSystem>
        void
        initStep(const ReservoirState& state,
                 const Grid&           g    ,
                 JacobianSystem&       sys  ) {
            (void) state;  (void) g; // Suppress 'unused'

            sys.vector().solution().fill(0.0);
        }

        template <class ReservoirState,
                  class Grid          ,
                  class JacobianSystem>
        void
        initIteration(const ReservoirState& state,
                      const Grid&           g    ,
                      JacobianSystem&       sys  ) {
            std::array<double, 2*2> dmob;

            const double *s = &state.saturation[0*2 + 0];

            for (int c = 0; c < g.number_of_cells; ++c, s += 0) {
                this->mobility(c, s, store_.mob(c), dmob);

                store_.dmob(c)[0] =  dmob[0*2 + 0];
                store_.dmob(c)[1] = -dmob[1*2 + 1];
            }
        }

        template <class ReservoirState,
                  class Grid          ,
                  class NewtonIterate >
        void
        finishIteration(const ReservoirState& state,
                        const Grid&           g    ,
                        NewtonIterate&        it   ) {
            // Nothing to do at end of iteration in this model.
            (void) state;  (void) g;  (void) it;
        }

        template <class Grid          ,
                  class SolutionVector,
                  class ReservoirState>
        void
        finishStep(const Grid&           g    ,
                   const SolutionVector& x    ,
                   ReservoirState&       state) {

            double *s = &state.saturation[0*2 + 0];

            for (int c = 0; c < g.number_of_cells; ++c, s += 2) {
                s[0] += x[c]    ;
                s[1]  = 1 - s[0];
            }
        }

    private:
        void
        upwindMobility(const double dflux,
                       const double gflux,
                       const int*   n    ,
                       int*         pix  ,
                       double*      m    ,
                       double*      dm   ) const {
            bool equal_sign = ( (! (dflux < 0)) && (! (gflux < 0)) ) ||
                              ( (! (dflux > 0)) && (! (gflux > 0)) );

            if (equal_sign) {

                if (! (dflux < 0) && ! (gflux < 0)) { pix[0] = 0; }
                else                                { pix[0] = 1; }

                m[0] = data_.mob(n[ pix[0] ]) [ 0 ];

                if (! (dflux - m[0]*gflux < 0))     { pix[1] = 0; }
                else                                { pix[1] = 1; }

                m[1] = data_.mob(n[ pix[1] ]) [ 1 ];

            } else {

                if (! (dflux < 0) && ! (gflux > 0)) { pix[1] = 0; }
                else                                { pix[1] = 1; }

                m[1] = data_.mob(n[ pix[1] ]) [ 1 ];

                if (dflux + m[1]*gflux > 0)         { pix[0] = 0; }
                else                                { pix[0] = 1; }

                m[0] = data_.mob(n[ pix[0] ]) [ 0 ];
            }

            dm[0] = data_.dmob(n[ pix[0] ]) [ 0 ];
            dm[1] = data_.dmob(n[ pix[1] ]) [ 1 ];
        }

        template <class Grid>
        void
        computeStaticGravity(const Grid&   g    ,
                             const double* grav ,
                             const double* trans) {
            const int d = g.dimensions;

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                const double* cc = g.cell_centroids + (c * d);

                for (; i < g.cell_facepos[c + 1]; ++i) {
                    const int     f  = g.cell_faces[i];
                    const double* fc = g.face_centroids + (f * d);

                    double dg = 0.0;
                    for (int j = 0; j < d; ++j) {
                        dg += grav[j] * (fc[j] - cc[j]);
                    }

                    data_.dg(i) = trans[f] * dg;
                }
            }
        }

        double
        gravityFlux(const int f) const {
            double gflux;

            if (gravity_) {
                int i1 = f2hf_[2*f + 0];
                int i2 = f2hf_[2*f + 1];

                assert ((i1 >= 0) && (i2 >= 0));

                gflux  = data_.dg(i1) - data_.dg(i2);
                gflux *= data_.drho();
            } else {
                gflux = 0.0;
            }

            return gflux;
        }

        bool                          gravity_;
        std::vector<int>              f2hf_   ;
        spu_2p::ModelParameterStorage store_  ;
    };
}
#endif  /* OPM_SINGLEPOINTUPWINDTWOPHASE_HPP_HEADER */
