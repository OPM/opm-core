/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "call_umfpack.h"
#include "spu_implicit.h"




/* Assume uniformly spaced table. */
double interpolate(int n, double h, double x0, double *tab, double x)
{
    int           i;
    double        a;

    assert(h > 0);
    assert((x-x0) < h*INT_MAX);
    assert((x-x0) > h*INT_MIN);

    if ( x < x0  ) { 
        return tab[0];
    }
    
    i = ((x-x0)/h);

    assert(i>=0);
    
    if (i+1 > n-1) { 
        return tab[n-1];
    }
    
    a = (x-x0 - i*h) / h;
    
    return (1-a) * tab[i] + a * tab[i+1];
}

/* Assume uniformly spaced table. */
double differentiate(int n, double h, double x0, double *tab, double x)
{
    int           i;
    double        a;

    assert(h > 0);
    assert((x-x0) < h*INT_MAX);
    assert((x-x0) > h*INT_MIN);
    assert(n>1);

    if ( x < x0  ) { 
        return (tab[1]-tab[0])/h;
    }
    
    i = ((x-x0)/h);

    assert(i>=0);
    
    if (i+1 > n-1) { 
        return (tab[n-1]-tab[n-2])/h;
    }
    
    a = (x-x0 - i*h) / h;
    
    return (tab[i+1]-tab[i])/h;
}

void compute_mobilities(int n, double *s, double *mob, double *dmob, int ntab, double h, double x0, double *tab)
{
    double *tabw=tab;
    double *tabo=tab+ntab;
    
    int i;
    for (i=0; i<n; ++i) {
        *mob++  = interpolate  (ntab, h, x0, tabw, *s);
        *mob++  = interpolate  (ntab, h, x0, tabo, *s);
        *dmob++ = differentiate(ntab, h, x0, tabw, *s);
        *dmob++ = differentiate(ntab, h, x0, tabo, *s++);
    }
}


typedef struct Sparse {
    int     m;
    int     n;
    int    *ia;
    int    *ja;
    double *sa;
} sparse_t;

#if 0

/*
 * To compute the derivative of the (mobility-weighted) upwind flux in
 * terms of the mobilities and the first derivatives of mobilities, we
 * use
 *
 *             mw
 *    F  =  -------
 *          mo + mw
 *
 *
 *    dF      dt      mo       1
 *    --- =   -- · -------  -----  (dflux + mo*gflux)
 *    dmw     pv   mw + mo  mw + mo
 * 
 *            dt  fo
 *        =   --  --  (dflux + mo*gflux)
 *            pv  mt
 *   
 *   
 *    dF          dt  fw
 *    --- = -1 *  --  --  (dflux - mw*gflux)
 *    dmo         pv  mt
 *
 *
 * and the chain rule  
 *
 *    dF    dF  dmw    dF  dmo
 *    --  = --- ---  + --- ---.
 *    dS    dmw dS     dmo dS
 *
 *
 * If the gravity term is zero, we can use the simpler rule 
 *
 *    dF    mo * dmw - mw * dmo
 *    --  = -------------------
 *    dS       [ mw + mo]²
 *
 *
 *
*/ 
double 
spu_implicit(grid_t *g, double *s0, double *s, double *mob, double *dmob,
             double *dflux, double *gflux, double *src, double dt)
{
    int    i, k, f, c1, c2, c;
    int    nc   = g->number_of_cells;
    int    nf   = g->number_of_faces;
    int    nhf  = g->cell_facepos[nc]; /* Too much */

    double m1,  m2;
    double dm1, dm2;
    double mt2;
    
    double *b;
    double *x;
    double infnorm;
    
    int           *pja;
    double        *psa;
    double        *d;
    double        *p1, *p2;
    double         flux;
    double         sgn;
    double         darcyflux;
    double         gravityflux;
    sparse_t      *S = malloc(sizeof *S);
    if (S) {
        S->ia = malloc((nc+1) * sizeof *S->ia);
        S->ja = malloc( nhf   * sizeof *S->ja);
        S->sa = malloc( nhf   * sizeof *S->sa);
    }
    else {
        assert(0);
    }
    S->n = nc;

    b = malloc(nc * sizeof *b);
    x = malloc(nc * sizeof *x);

    pja = S->ja;
    psa = S->sa;
    
    /* Assemble system */
    S->ia[0] = 0;    
    for (i=0; i<nc; ++i) {

        /* Store position of diagonal element */
        d      = psa++;
        *pja++ = i;
        *d     = 0.0;

        /* Accumulation term */
        b[i]   = (s[i] - s0[i]);        
        *d    += 1.0;
        
        /* Flux terms follows*/
        for (k=g->cell_facepos[i]; k<g->cell_facepos[i+1]; ++k) {
            f   = g->cell_faces[k];
            c1  = g->face_cells[2*f+0];
            c2  = g->face_cells[2*f+1];

            /* Skip all boundary terms (for now). */
            if (c1 == -1 || c2 == -1) { continue; } 

            /* Set cell index of other cell, set correct sign of fluxes. */
            c            = (i==c1 ? c2 : c1);
            sgn          = (i==c1)*2.0 - 1.0;
            darcyflux    = sgn * dt*dflux[f];
            gravityflux  = sgn * dt*gflux[f];

            /* ====================================================== */
            /* If darcy flux and gravity flux have same sign...       */
            if ( (darcyflux>0.0 && gravityflux>0.0) ||
                 (darcyflux<0.0 && gravityflux<0.0) ) {

                /* Set water mobility */
                if (darcyflux>0) {
                    /* Positive water phase flux */
                    m1  =  mob[2*i];
                    dm1 = dmob[2*i];
                    p1  = d;
                }
                else {
                    /* Negative water phase flux */
                    m1  =  mob[2*c];
                    dm1 = dmob[2*c];
                    p1  = psa;
                }

                /* Set oil mobility */
                if (darcyflux - m1*gravityflux>0) {
                    /* Positive oil phase flux */
                    m2  =  mob[2*i+1];
                    dm2 = dmob[2*i+1];
                    p2  = d;
                }
                else {
                    /* Negative oil phase flux */
                    m2  =  mob[2*c+1];
                    dm2 = dmob[2*c+1];
                    p2  = psa;
                }
            }
            
            /* ====================================================== */
            /* If Darcy flux and gravity flux have opposite sign...   */
            else {
                /* Set oil mobility */
                if (darcyflux>0) {
                    /* Positive oil phase flux */
                    m2  =  mob[2*i+1];
                    dm2 = dmob[2*i+1];
                    p2  = d;
                }
                else {
                    /* Negative oil phase flux */
                    m2  =  mob[2*c+1];
                    dm2 = dmob[2*c+1];
                    p2  = psa;
                }
                /* Set water mobility */
                if (darcyflux+m2*gravityflux>0) {
                    /* Positive water phase flux */
                    m1  =  mob[2*i];
                    dm1 = dmob[2*i];
                    p1  = d;
                }
                else {
                    /* Negative water phase flux */
                    m1  =  mob[2*c];
                    dm1 = dmob[2*c];
                    p1  = psa;
                }
            }
            
            b[i]  += m1/(m1+m2)*(darcyflux + m2*gravityflux);
            if (p1 == psa || p2 == psa) {
                *psa++ = 0.0;
                *pja++ = c;
            }
            
            /* Add contribitions from dFw/dmw·dmw/dsw and dFw/dmo·dmo/dsw */
            mt2  = (m1+m2)*(m1+m2);

            /* dFw/dmw·dmw/dsw */
            *p1 +=  m2/mt2*(darcyflux + m2*gravityflux)*dm1;

            /* dFw/dmo·dmo/dsw */
            *p2 += -m1/mt2*(darcyflux - m1*gravityflux)*dm2;
            
            
        }
        /* Injection */
        if (src[i] > 0.0) {
            /* Assume sat==1.0 in source, and f(1.0)=1.0; */
            m1    = 1.0;
            m2    = 0.0;
            b[i] -= dt*src[i] * m1/(m1+m2);
        }
        
        /* Production */
        else {
            m1  = mob [2*i+0];
            m2  = mob [2*i+1];
            dm1 = dmob[2*i+0];
            dm2 = dmob[2*i+1];
            
            *d   -= dt*src[i] *(m2*dm1-m1*dm2)/(m1+m2)/(m1+m2);
            b[i] -= dt*src[i] *m1/(m1+m2);
            
        }
        S->ia[i+1] = pja - S->ja;
    }

    d = malloc(nc * sizeof *d);
    
    
    infnorm = 0.0;
    for(i=0; i<nc; ++i) {infnorm = (infnorm > fabs(b[i]) ? infnorm : fabs(b[i]));}
    /* fprintf(stderr, "norm %e\n", infnorm); */
    
    
#if 0
    for(i=0; i<nc; ++i) {
        for(f=0; f<nc; ++f) {d[f] = 0.0;}
        for(k=S->ia[i]; k<S->ia[i+1]; ++k) {
            d[S->ja[k]] = S->sa[k];
        }
        for(f=0; f<nc; ++f) {if (d[f]==0.0) fprintf(stderr, "    .   ");else fprintf(stderr, "%+6.3e ", d[f]);}
        fprintf(stderr, "\n");
    }
#endif
#if 0
    fprintf(stderr, "b: ");
    for(i=0; i<nc; ++i) {fprintf(stderr, "%+6.3e ", b[i]);}
    fprintf(stderr, "\n");

#endif
#if 1
    /* Solve system */
    callMWUMFPACK(S->n, S->ia, S->ja, S->sa, b, x);


    /* Return x. */
    for(i=0; i<nc; ++i) {
        s[i] = s[i] - x[i];
        /* fprintf(stderr, "%e %e\n", s[i], x[i]); */
    }
    /* fprintf(stderr, "\n"); */
        
#endif    
    


    free(b); 
    free(x);
    free(S->ia);
    free(S->ja);
    free(S->sa);
    free(S);

    return infnorm;
}   
#endif

double 
spu_implicit2(grid_t *g, double *s0, double *s, double *mob, double *dmob,
              double *dflux, double *gflux, double *src, double dt)
{
    int    nc   = g->number_of_cells;
    int    nf   = g->number_of_faces;
    int    nhf  = g->cell_facepos[nc]; /* Too much */
    double infnorm;
    double *b; 
    double *x;
    int i;

    /* Allocate space for linear system etc.*/
    sparse_t      *S = malloc(sizeof *S);
    if (S) {
        S->ia = malloc((nc+1) * sizeof *S->ia);
        S->ja = malloc( nhf   * sizeof *S->ja);
        S->sa = malloc( nhf   * sizeof *S->sa);
        S->n  = nc;
    }
    else {
        assert(0);
    }

    b = malloc(nc * sizeof *b);
    x = malloc(nc * sizeof *x);


    spu_implicit_assemble(g, s0, s, mob, dmob, dflux, gflux, src, dt, S, b);

    /* Compute inf-norm of residual */
    infnorm = 0.0;
    for(i=0; i<nc; ++i) {
        infnorm = (infnorm > fabs(b[i]) ? infnorm : fabs(b[i]));
    }

    /* Solve system */
    callMWUMFPACK(S->n, S->ia, S->ja, S->sa, b, x);


    /* Return x. */
    for(i=0; i<nc; ++i) {
        s[i] = s[i] - x[i];
    }
        

    free(b); 
    free(x);
    free(S->ia);
    free(S->ja);
    free(S->sa);
    free(S);

    return infnorm;    
}


void 
spu_implicit_assemble(grid_t *g, double *s0, double *s, double *mob, double *dmob,
                      double *dflux, double *gflux, double *src, double dt, sparse_t *S, double *b)
{
    int     i, k, f, c1, c2, c;
    int     nc   = g->number_of_cells;
    /* int    nf   = g->number_of_faces; */
    /* int    nhf  = g->cell_facepos[nc]; /\* Too much *\/ */

    double  m1,  m2;
    double  dm1, dm2;
    double  mt2;
    
    int    *pja;
    double *psa;
    double *d;
    double *p1, *p2;
    double  flux;
    double  sgn;
    double  darcyflux;
    double  gravityflux;

    pja = S->ja;
    psa = S->sa;
    
    /* Assemble system */
    S->ia[0] = 0;    
    for (i=0; i<nc; ++i) {

        /* Store position of diagonal element */
        d      = psa++;
        *pja++ = i;
        *d     = 0.0;

        /* Accumulation term */
        b[i]   = (s[i] - s0[i]);        
        *d    += 1.0;
        
        /* Flux terms follows*/
        for (k=g->cell_facepos[i]; k<g->cell_facepos[i+1]; ++k) {
            f   = g->cell_faces[k];
            c1  = g->face_cells[2*f+0];
            c2  = g->face_cells[2*f+1];

            /* Skip all boundary terms (for now). */
            if (c1 == -1 || c2 == -1) { continue; } 

            /* Set cell index of other cell, set correct sign of fluxes. */
            c            = (i==c1 ? c2 : c1);
            sgn          = (i==c1)*2.0 - 1.0;
            darcyflux    = sgn * dt*dflux[f];
            gravityflux  = sgn * dt*gflux[f];

            /* ====================================================== */
            /* If darcy flux and gravity flux have same sign...       */
            if ( (darcyflux>0.0 && gravityflux>0.0) ||
                 (darcyflux<0.0 && gravityflux<0.0) ) {

                /* Set water mobility */
                if (darcyflux>0) {
                    /* Positive water phase flux */
                    m1  =  mob[2*i];
                    dm1 = dmob[2*i];
                    p1  = d;
                }
                else {
                    /* Negative water phase flux */
                    m1  =  mob[2*c];
                    dm1 = dmob[2*c];
                    p1  = psa;
                }

                /* Set oil mobility */
                if (darcyflux - m1*gravityflux>0) {
                    /* Positive oil phase flux */
                    m2  =  mob[2*i+1];
                    dm2 = dmob[2*i+1];
                    p2  = d;
                }
                else {
                    /* Negative oil phase flux */
                    m2  =  mob[2*c+1];
                    dm2 = dmob[2*c+1];
                    p2  = psa;
                }
            }
            
            /* ====================================================== */
            /* If Darcy flux and gravity flux have opposite sign...   */
            else {
                /* Set oil mobility */
                if (darcyflux>0) {
                    /* Positive oil phase flux */
                    m2  =  mob[2*i+1];
                    dm2 = dmob[2*i+1];
                    p2  = d;
                }
                else {
                    /* Negative oil phase flux */
                    m2  =  mob[2*c+1];
                    dm2 = dmob[2*c+1];
                    p2  = psa;
                }
                /* Set water mobility */
                if (darcyflux+m2*gravityflux>0) {
                    /* Positive water phase flux */
                    m1  =  mob[2*i];
                    dm1 = dmob[2*i];
                    p1  = d;
                }
                else {
                    /* Negative water phase flux */
                    m1  =  mob[2*c];
                    dm1 = dmob[2*c];
                    p1  = psa;
                }
            }
            
            b[i]  += m1/(m1+m2)*(darcyflux + m2*gravityflux);
            if (p1 == psa || p2 == psa) {
                *psa++ = 0.0;
                *pja++ = c;
            }
            
            /* Add contribitions from dFw/dmw·dmw/dsw and dFw/dmo·dmo/dsw */
            mt2  = (m1+m2)*(m1+m2);

            /* dFw/dmw·dmw/dsw */
            *p1 +=  m2/mt2*(darcyflux + m2*gravityflux)*dm1;

            /* dFw/dmo·dmo/dsw */
            *p2 += -m1/mt2*(darcyflux - m1*gravityflux)*dm2;
            
            
        }
        /* Injection */
        if (src[i] > 0.0) {
            /* Assume sat==1.0 in source, and f(1.0)=1.0; */
            m1    = 1.0;
            m2    = 0.0;
            b[i] -= dt*src[i] * m1/(m1+m2);
        }
        
        /* Production */
        else {
            m1  = mob [2*i+0];
            m2  = mob [2*i+1];
            dm1 = dmob[2*i+0];
            dm2 = dmob[2*i+1];
            
            *d   -= dt*src[i] *(m2*dm1-m1*dm2)/(m1+m2)/(m1+m2);
            b[i] -= dt*src[i] *m1/(m1+m2);
            
        }
        S->ia[i+1] = pja - S->ja;
    }
}   
