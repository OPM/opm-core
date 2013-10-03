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

/**
 * \file
 * Numerical model and support classes needed to model transport of two
 * incompressible fluid phases.  Intended for the ImplicitTransport system.
 */

#ifndef OPM_SINGLEPOINTUPWINDTWOPHASE_HPP_HEADER
#define OPM_SINGLEPOINTUPWINDTWOPHASE_HPP_HEADER

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <vector>
#include <iostream>

namespace Opm {
    namespace spu_2p {
        /**
         * Internal class to manage the direct and derived quantities needed to
         * formulate the fluid transport system.
         *
         * Note: This class elides off-diagonal elements of the phase mobility
         * Jacobian, \f$(\partial_{s_\beta} \lambda_\alpha)_{\alpha\beta}\f$.
         * These elements are assumed to be strictly equal to zero.  In other
         * words, the relative permeability of phase \f$\alpha\f$ is assumed to
         * depend only on the saturation of phase \f$\alpha\f$.  This convention
         * allows storing only the two diagonals of the mobility Jacobian per
         * grid cell.
         *
         * The static gravity term is the scalar value
         * \f[
         * \Delta G_i = \mathsf{T}_f\, \vec{g}\cdot(\Bar{x}_f - \Bar{x}_c)
         * \f]
         * in which @c i is the half face index corresponding to the cell-face
         * pair <CODE>(f,c)</CODE> and \f$\mathsf{T}_f\f$ is the absolute
         * (bacground) two-point transmissibility of face @c f.
         *
         * The fluid transport problem is formulated in terms of saturation
         * changes, \f$\Delta s\f$, per cell.  These changes are the primary
         * degrees of freedom in this model.
         *
         * Capillary pressures are defined by the fluid model, but usually
         * correspond to \f$p_w - p_n\f$ (e.g., \f$p_\mathit{oil} -
         * p_\mathit{water}\f$).
         */
        class ModelParameterStorage {
        public:
            /**
             * Constructor.
             *
             * @param[in] nc      Total number of grid cells.
             * @param[in] totconn Total number of connections, accumulated per
             *                    cell (``half faces'').
             */
            ModelParameterStorage(int nc, int totconn)
                : drho_(0.0), mob_(0), dmob_(0),
                  porevol_(0), dg_(0), ds_(0), pc_(0), dpc_(0), trans_(0),
                  data_()
            {
                std::size_t alloc_sz;

                alloc_sz  = 2 * nc;      // mob_
                alloc_sz += 2 * nc;      // dmob_
                alloc_sz += 1 * nc;      // porevol_
                alloc_sz += 1 * totconn; // dg_
                alloc_sz += 1 * nc;      // ds_
                alloc_sz += 1 * nc;      // pc_
                alloc_sz += 1 * nc;      // dpc_
                alloc_sz += 1 * totconn; // dtrans
                data_.resize(alloc_sz);

                mob_     = &data_[0];
                dmob_    = mob_     + (2 * nc     );
                porevol_ = dmob_    + (2 * nc     );
                dg_      = porevol_ + (1 * nc     );
                ds_      = dg_      + (1 * totconn);
                pc_      = ds_      + (1 * nc     );
                dpc_     = pc_      + (1 * nc     );
                trans_   = dpc_      + (1 * nc     );
            }

            /**
             * Modifiable density difference.
             * @return Reference to modifiable internal representation of fluid
             * phase density difference.
             */
            double&       drho   ()            { return drho_            ; }
            /**
             * Read-only density difference.
             * @return Read-only value of current fluid phase difference value.
             */
            double        drho   ()      const { return drho_            ; }

            /**
             * Phase mobility in cell.
             * @param[in] c Cell.
             * @return Read-write reference to two consecutive phase mobilities
             * in cell @c c.
             */
            double*       mob    (int c)       { return mob_  + (2*c + 0); }
            /**
             * Phase mobility in cell.
             * @param[in] c Cell.
             * @return Read-only reference to two consecutive phase mobilities
             * in cell @c c.
             */
            const double* mob    (int c) const { return mob_  + (2*c + 0); }

            /**
             * Diagonal elements of phase mobility derivative (Jacobian).
             *
             * @param[in] c Cell.
             * @return Read-write reference to diagonal elements of phase
             * mobility Jacobian in cell @c c.
             */
            double*       dmob   (int c)       { return dmob_ + (2*c + 0); }
            /**
             * Diagonal elements of phase mobility derivative (Jacobian).
             *
             * @param[in] c Cell.
             * @return Read-only reference to two consecutive diagonal elements
             * of phase mobility Jacobian in cell @c c.
             */
            const double* dmob   (int c) const { return dmob_ + (2*c + 0); }

            /**
             * Retrieve pore volumes for all cells.
             * @return Modifiable vector of pore volumes for all cells.
             */
            double*       porevol()            { return porevol_         ; }
            /**
             * Pore volume of single cell.
             * @param[in] c Cell.
             * @return Pore volume of cell @c c.
             */
            double        porevol(int c) const { return porevol_[c]      ; }

            /**
             * Static gravity term associated to single half face.
             *
             * @param[in] i Half face index corresponding to particular
             *              cell-face pair.
             * @return Read-write reference to static gravity term of single
             * half face.
             */
            double&       dg(int i)            { return dg_[i]           ; }
            /**
             * Static gravity term associated to single half face.
             * @param[in] i Half face index corresponding to particular
             *              cell-face pair.
             * @return Read-only reference to static gravity term of single half
             * face.
             */
            double        dg(int i)      const { return dg_[i]           ; }

            /**
             * Saturation change in particular cell.
             *
             * @param[in] c
             * @return Read-write reference to saturation change (scalar) in
             * cell @c c.
             */
            double&       ds(int c)            { return ds_[c]           ; }
            /**
             * Saturation change in particular cell.
             *
             * @param[in] c
             * @return Read-only reference to saturation change (scalar) in cell
             * @c c.
             */
            double        ds(int c)      const { return ds_[c]           ; }

            /**
             * Capillary pressure in particular cell.
             *
             * @param[in] c Cell.
             * @return Read-write reference to capillary pressure in cell @c c.
             */
            double&       pc(int c)            { return pc_[c]           ; }
            /**
             * Capillary pressure in particular cell.
             *
             * @param[in] c Cell
             * @return Read-only reference to capillary pressure in cell @c c.
             */
            double        pc(int c)      const { return pc_[c]           ; }

            /**
             * Derivative of capillary pressure with respect to saturation.
             *
             * @param[in] c Cell
             * @return Read-write reference to capillary pressure derivative
             * with respect to primary saturation in cell @c c.
             */
            double&       dpc(int c)           { return dpc_[c]          ; }
            /**
             * Derivative of capillary pressure with respect to saturation.
             *
             * @param[in] c Cell
             * @return Read-only reference to capillary pressure derivative with
             * respect to primary saturation in cell @c c.
             */
            double        dpc(int c)     const { return dpc_[c]          ; }

            /**
             * Background (absolute) face transmissibility of particular face.
             *
             * @param[in] f Face
             * @return Read-write reference to background face transmissibility
             * of face @c f.
             */
            double&       trans(int f)         { return trans_[f]        ; }
            /**
             * Background (absolute) face transmissibility of particular face.
             *
             * @param[in] f Face
             * @return Read-only reference to bacground face transmissibility of
             * face @c f.
             */
            double        trans(int f)   const { return trans_[f]        ; }

        private:
            double  drho_   ;  /**< Fluid phase density difference */
            double *mob_    ;  /**< Fluid phase mobility in all cells */
            double *dmob_   ;  /**< Derivative of phase mobility in all cells */
            double *porevol_;  /**< Pore volume in all cells */
            double *dg_     ;  /**< Static gravity term on all half faces */
            double *ds_     ;  /**< Saturation change in all cells */
            double *pc_     ;  /**< Capillary pressure in all cells */
            double *dpc_    ;  /**< Derivative of cap. pressure in all cells */
            double *trans_  ;  /**< Absolute transmissibility on all faces */

            /**
             * Data storage from which individual quantities are managed.
             */
            std::vector<double> data_;
        };
    }


    template <class TwophaseFluid>
    class SinglePointUpwindTwoPhase {
    public:
        template <class Grid>
        SinglePointUpwindTwoPhase(const TwophaseFluid&       fluid    ,
                                  const Grid&                g        ,
                                  const std::vector<double>& porevol  ,
                                  const double*              grav  = 0,
				  const bool                 guess_previous = true)
            : fluid_  (fluid)                              ,
              gravity_(grav)        ,
              f2hf_   (2 * g.number_of_faces, -1)          ,
              store_  (g.number_of_cells,
                       g.cell_facepos[ g.number_of_cells ]),
	      init_step_use_previous_sol_(guess_previous),
	      sat_tol_(1e-5)
        {

            if (gravity_) {
                store_.drho() = fluid_.density(0) - fluid_.density(1);
            }

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                for (; i < g.cell_facepos[c + 1]; ++i) {
                    const int f = g.cell_faces[i];
                    const int p = 1 - (g.face_cells[2*f + 0] == c);
                    f2hf_[2*f + p] = i;
                }
            }

            std::copy(porevol.begin(), porevol.end(), store_.porevol());
        }

        void
        makefhfQPeriodic(const std::vector<int>& p_faces ,
                         const std::vector<int>& hf_faces,
                         const std::vector<int>& nb_faces)
        {
	    if (p_faces.empty()) {
		return;
	    }
            assert (p_faces.size()  == hf_faces.size());
            assert (hf_faces.size() == nb_faces.size());

            std::vector<int> nbhf(hf_faces.size());

            for (std::vector<int>::size_type i = 0; i < p_faces.size(); ++i) {
                const int nbf = nb_faces[i];

                assert (2*std::vector<int>::size_type(nbf) + 1 < f2hf_.size());
                assert ((f2hf_[2*nbf + 0] < 0) ^ (f2hf_[2*nbf + 1] < 0));

                const int p = (f2hf_[2*nbf + 0] < 0) ? 1 : 0;  // "Self"
                nbhf[ i ]   =  f2hf_[2*nbf + p];
            }

            for (std::vector<int>::size_type i = 0; i < p_faces.size(); ++i) {
                const int f  = p_faces [i];
                const int hf = hf_faces[i];

                assert (0 <= f);
                assert (0 <= hf);
                assert (2*std::vector<int>::size_type(f) + 1 < f2hf_.size());

                assert ((f2hf_[2*f + 0] <  0 ) ^ (f2hf_[2*f + 1] <  0 ));
                assert ((f2hf_[2*f + 0] == hf) ^ (f2hf_[2*f + 1] == hf));

                const int p = (f2hf_[2*f + 0] == hf) ? 1 : 0;  // "Other"

                f2hf_[2*f + p] = nbhf[ i ];
            }
        }

        // -----------------------------------------------------------------
        // System assembly innards
        // -----------------------------------------------------------------

        enum { DofPerCell = 1 };

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
            double dflux = state.faceflux()[f];
            double gflux = gravityFlux(f);
            double pcflux,dpcflux[2];
            capFlux(f,n, pcflux, dpcflux);
            gflux += pcflux;

            int    pix[2];
            double m[2], dm[2];
            upwindMobility(dflux, gflux, n, pix, m, dm);

            assert (! ((m[0] < 0) || (m[1] < 0)));

            double mt = m[0] + m[1];
            assert (mt > 0);

            double sgn  = 2.0*(n[0] == c) - 1.0;
            dflux      *= sgn;
            gflux      *= sgn;


            double       f1 = m[0] / mt;
            const double v1 = dflux + m[1]*gflux;

            // Assemble residual contributions
            *F += dt * f1 * v1;

            // Assemble Jacobian (J1 <-> c, J2 <-> other)
            double *J[2];
            if (n[0] == c) {
                J[0] = J1; J[1] = J2;
                // sign is positive
                J1[0*2 + 0] += sgn*dt * f1            * dpcflux[0] * m[1];
                J2[0*2 + 0] += sgn*dt * f1            * dpcflux[1] * m[1];
            } else {
                J[0] = J2; J[1] = J1;
                // sign is negative
                J1[0*2 + 0] += sgn*dt * f1            * dpcflux[1] * m[1];
                J2[0*2 + 0] += sgn*dt * f1            * dpcflux[0] * m[1];
            }

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

            const double pv = store_.porevol(c);

            *J += pv;
            *F += pv * store_.ds(c);
        }

        template <class Grid       ,
                  class SourceTerms>
        void
        sourceTerms(const Grid&        g  ,
                    const SourceTerms* src,
                    const int          i  ,
                    const double       dt ,
                    double*            J  ,
                    double*            F  ) const {

            (void) g;

            double dflux = -src->flux[i]; // ->flux[] is rate of *inflow*

            if (dflux < 0) {
                // src -> cell, affects residual only.
                *F += dt * dflux * src->saturation[2*i + 0];
            } else {
                // cell -> src
                const int     c  = src->cell[i];
                const double* m  = store_.mob (c);
                const double* dm = store_.dmob(c);

                const double  mt = m[0] + m[1];

                assert (! ((m[0] < 0) || (m[1] < 0)));
                assert (mt > 0);

                const double f  = m[0] / mt;
                const double df = ((1 - f)*dm[0] - f*dm[1]) / mt;

                *F += dt * dflux *  f;
                *J += dt * dflux * df;
            }
        }
        template <class Grid>
        void
        initGravityTrans(const Grid&  g    ,
                         const std::vector<double> &  htrans) {

            assert (htrans.size() ==
                    static_cast<std::vector<double>::size_type>(g.cell_facepos[ g.number_of_cells ]));

            for (int f = 0; f < g.number_of_faces; ++f) {
                store_.trans(f) = 0.0;
            }

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                for (; i < g.cell_facepos[c + 1]; ++i) {
                    int f = g.cell_faces[i];

                    assert (htrans[i] > 0.0);

                    store_.trans(f) += 1.0 / htrans[i];
                }
            }

            for (int f = 0; f < g.number_of_faces; ++f) {
                store_.trans(f) = 1.0 / store_.trans(f);
            }

            if (gravity_) {
                this->computeStaticGravity(g);
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

            (void) state;       // Suppress 'unused' warning.

            typename JacobianSystem::vector_type& x =
                sys.vector().writableSolution();

            assert (x.size() == (::std::size_t) (g.number_of_cells));

	    if (init_step_use_previous_sol_) {
		std::fill(x.begin(), x.end(), 0.0);
            } else {
		const std::vector<double>& s = state.saturation();
		for (int c = 0, nc = g.number_of_cells; c < nc; ++c) {
		    // Impose s=0.5 at next time level as an NR initial value.
		    x[c] = 0.5 - s[2*c + 0];
		}
	    }
        }

        template <class ReservoirState,
                  class Grid          ,
                  class JacobianSystem>
        bool
        initIteration(const ReservoirState& state,
                      const Grid&           g    ,
                      JacobianSystem&       sys) {

            double s[2],  mob[2],  dmob[2 * 2], pc, dpc;

            const typename JacobianSystem::vector_type& x =
                sys.vector().solution();
            const ::std::vector<double>& sat = state.saturation();

            bool in_range = true;
            for (int c = 0; c < g.number_of_cells; ++c) {
                store_.ds(c) = x[c]; // Store sat-change for accumulation().

                s[0] = sat[c*2 + 0] + x[c];

                double s_min = fluid_.s_min(c);
                double s_max = fluid_.s_max(c);

                if ( s[0] < (s_min - sat_tol_) || s[0] > (s_max + sat_tol_) ) {
                    // if (s[0] < s_min){
		    // 	std::cout << "Warning: s out of range, s-s_min = " << s_min-s[0] << std::endl;
                    // }
                    // if (s[0] > s_max){
		    // 	std::cout << "Warning: s out of range, s-s_max = " << s[0]-s_max << std::endl;
                    // }
                    in_range = false; //line search fails
                }
                s[0] = std::max(s_min, s[0]);
                s[0] = std::min(s_max, s[0]);
                s[1] = 1 - s[0];

                fluid_.mobility(c, s, mob, dmob);
                fluid_.pc(c, s, pc, dpc);

                store_.mob (c)[0] =  mob [0];
                store_.mob (c)[1] =  mob [1];
                store_.dmob(c)[0] =  dmob[0*2 + 0];
                store_.dmob(c)[1] = -dmob[1*2 + 1];
                store_.pc(c)      = pc;
                store_.dpc(c)     = dpc;
            }
            if (!in_range) {
#ifdef VERBOSE
                std::cout << "Warning: initIteration() - s was clamped in some cells.\n";
#endif
            }
            return in_range;
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

            double *s = &state.saturation()[0*2 + 0];

            for (int c = 0; c < g.number_of_cells; ++c, s += 2) {
                s[0] += x[c];
                double s_min = fluid_.s_min(c);
                double s_max = fluid_.s_max(c);

#if 0
                assert(s[0] >= s_min - sat_tol_);
                assert(s[0] <= s_max + sat_tol_);
#endif

                s[0] = std::max(s_min, s[0]);
                s[0] = std::min(s_max, s[0]);
                s[1]  = 1.0 - s[0];
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

                m[0] = store_.mob(n[ pix[0] ]) [ 0 ];

                if (! (dflux - m[0]*gflux < 0))     { pix[1] = 0; }
                else                                { pix[1] = 1; }

                m[1] = store_.mob(n[ pix[1] ]) [ 1 ];

            } else {

                if (! (dflux < 0) && ! (gflux > 0)) { pix[1] = 0; }
                else                                { pix[1] = 1; }

                m[1] = store_.mob(n[ pix[1] ]) [ 1 ];

                if (dflux + m[1]*gflux > 0)         { pix[0] = 0; }
                else                                { pix[0] = 1; }

                m[0] = store_.mob(n[ pix[0] ]) [ 0 ];
            }

            dm[0] = store_.dmob(n[ pix[0] ]) [ 0 ];
            dm[1] = store_.dmob(n[ pix[1] ]) [ 1 ];
        }

        template <class Grid>
        void
        computeStaticGravity(const Grid& g) {

            const int d = g.dimensions;

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                const double* cc = g.cell_centroids + (c * d);

                for (; i < g.cell_facepos[c + 1]; ++i) {
                    const int     f  = g.cell_faces[i];
                    const double* fc = g.face_centroids + (f * d);

                    double dg = 0.0;
                    for (int j = 0; j < d; ++j) {
                        dg += gravity_[j] * (fc[j] - cc[j]);
                    }

                    store_.dg(i) = store_.trans(f) * dg;
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

                gflux  = store_.dg(i1) - store_.dg(i2);
                gflux *= store_.drho();
            } else {
                gflux = 0.0;
            }

            return gflux;
        }
        void
        capFlux(const int f,const int* n,double& pcflux, double* dpcflux) const {
            //double capflux;
            int i1 = n[0];
            int i2 = n[1];
            assert ((i1 >= 0) && (i2 >= 0));
            //double sgn=-1.0;
            pcflux  = store_.trans(f)*(store_.pc(i2) - store_.pc(i1));
            dpcflux[0]  = -store_.trans(f)*store_.dpc(i1);
            dpcflux[1]  = store_.trans(f)*store_.dpc(i2);
        }

        TwophaseFluid                 fluid_  ;
        const double*                 gravity_;
        std::vector<int>              f2hf_   ;
        spu_2p::ModelParameterStorage store_  ;
	bool init_step_use_previous_sol_;
	double sat_tol_;
    };
}
#endif  /* OPM_SINGLEPOINTUPWINDTWOPHASE_HPP_HEADER */
