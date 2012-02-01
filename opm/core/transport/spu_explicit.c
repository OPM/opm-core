/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */


#include <assert.h>
#include <stdio.h>

#include <opm/core/grid.h>
#include <opm/core/transport/spu_explicit.h>


/* Twophase mobility-weighted upwind */
void
spu_explicit(struct UnstructuredGrid *g, double *s0, double *s, double *mob,
             double *dflux, double *gflux, double *src,
             double dt)
{
    int i, f, c1, c2;
    int nc = g->number_of_cells;
    int nf = g->number_of_faces;

    double m1, m2, flux;


    /* Contribution form sources */
    for (i=0; i<nc; ++i) {
        s[i] = s0[i];

        /* Injection */
        if (src[i] > 0.0) {
            /* Assume sat==1.0 in source, and f(1.0)=1.0; */
            s[i] += dt*src[i];
        }

        /* Production */
        else {
            m1 = mob[2*i];
            m2 = mob[2*i+1];
            s[i] += dt*src[i]* m1/(m1 + m2);
        }
    }

    /* Contribution from internal faces */
    for (f=0; f<nf; ++f) {
        c1 = g->face_cells[2*f+0];
        c2 = g->face_cells[2*f+1];
        if ((c1 !=-1) && (c2 !=-1)) {

            if ((dflux[f]>0.0 && gflux[f]>0.0) ||
                (dflux[f]<0.0 && gflux[f]<0.0) ) {
                /* Water mobility */
                if (dflux[f]>0) {
                    m1 = mob[2*c1];
                }
                else {
                    m1 = mob[2*c2];
                }
                /* Oil mobility */
                if (dflux[f] - m1*gflux[f]>0) {
                    m2 = mob[2*c1+1];
                }
                else {
                    m2 = mob[2*c2+1];
                }
            }

            else {
                /* Oil mobility */
                if (dflux[f]>0) {
                    m2 = mob[2*c1+1];
                }
                else {
                    m2 = mob[2*c2+1];
                }
                /* Water mobility */
                if (dflux[f]+m2*gflux[f]>0) {
                    m1 = mob[2*c1];
                }
                else {
                    m1 = mob[2*c2];
                }
            }

            /* Water flux */            
            assert(m1+m2>0.0);
            flux = m1/(m1+m2)*(dflux[f] + m2*gflux[f]);
            s[c1] -= flux*dt;
            s[c2] += flux*dt;
        }
    }
}   
