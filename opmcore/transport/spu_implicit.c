/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "grid.h"
#include "spu_implicit.h"




/* Assume uniformly spaced table. */
static double
interpolate(int n, double h, double x0, double *tab, double x)
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
static double 
differentiate(int n, double h, double x0, double *tab, double x)
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

static void 
compute_mobilities(int n, double *s, double *mob, double *dmob, int ntab, double h, double x0, double *tab)
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
spu_implicit(grid_t *g, double *s0, double *s, double h, double x0, int ntab, double *tab,
              double *dflux, double *gflux, double *src, double dt, 
              void (*linear_solver)(int, int*, int*, double *, double *, double *))
{
    int    nc   = g->number_of_cells;
    int    nhf  = g->cell_facepos[nc]; /* Too much */
    double infnorm;
    double *b; 
    double *x;
    char *work;
    double *mob, *dmob;
    int i;
    int it;

    /* Allocate space for linear system etc.*/
    sparse_t      *S = malloc(sizeof *S);
    if (S) {
        S->ia = malloc((nc+1) * sizeof *S->ia);
        S->ja = malloc( nhf   * sizeof *S->ja);
        S->sa = malloc( nhf   * sizeof *S->sa);
        S->n  = nc;
        S->m  = nc;
    }
    else {
        assert(0);
    }

    b = malloc(nc * sizeof *b);
    x = malloc(nc * sizeof *x);    

    work  = malloc(6*sizeof (double));
    mob   = malloc(g->number_of_cells *2* sizeof *mob);
    dmob  = malloc(g->number_of_cells *2* sizeof *dmob);
    
    infnorm = 1.0;
    it      = 0;
    while (infnorm > 1e-9 && it++ < 20) {
        compute_mobilities(g->number_of_cells, s, mob, dmob, ntab, h, x0, tab);
        spu_implicit_assemble(g, s0, s, mob, dmob, dflux, gflux, src, dt, S, b, work);
        
        /* Compute inf-norm of residual */
        infnorm = 0.0;
        for(i=0; i<nc; ++i) {
            infnorm = (infnorm > fabs(b[i]) ? infnorm : fabs(b[i]));
        }
        fprintf(stderr, "  Max norm of residual: %e\n", infnorm);
        
        /* Solve system */
        (*linear_solver)(S->m, S->ia, S->ja, S->sa, b, x);
        
        /* Return x. */
        for(i=0; i<nc; ++i) {
            s[i] = s[i] - x[i];
            s[i] = s[i]>1.0 ? 1.0 : (s[i] < 0.0 ? 0.0 : s[i]);
        }
    }

    free(work);
    free(mob);
    free(dmob);

    free(b); 
    free(x);
    free(S->ia);
    free(S->ja);
    free(S->sa);
    free(S);

    return infnorm;    
}


static void
phase_upwind_mobility(double darcyflux, double gravityflux,  int i, int c, 
                      double *mob, double *dmob, double *m, double *dm, int *cix)
{
    /* ====================================================== */
    /* If darcy flux and gravity flux have same sign...       */
    if ( (darcyflux>0.0 && gravityflux>0.0) ||
         (darcyflux<0.0 && gravityflux<0.0) ) {

        /* Positive water phase flux */
        if (darcyflux>0) {
            m [0]  = mob [2*i+0];
            dm[0]  = dmob[2*i+0];
            cix[0] = i;
        }
        
        /* Negative water phase flux */
        else {
            m [0]  = mob [2*c+0];
            dm[0]  = dmob[2*c+0];
            cix[0] = c;
        }

        /* Positive oil phase flux */
        if (darcyflux - m[0]*gravityflux>0) {
            m [1]  =  mob[2*i+1];
            dm[1]  = dmob[2*i+1];
            cix[1] = i;
        }

        /* Negative oil phase flux */
        else {
            m [1]  =  mob[2*c+1];
            dm[1]  = dmob[2*c+1];
            cix[1] = c;
        }
    }
            
    /* ====================================================== */
    /* If Darcy flux and gravity flux have opposite sign...   */
    else {
        
        /* Positive oil phase flux */
        if (darcyflux>0) {
            m [1]  =  mob[2*i+1];
            dm[1]  = dmob[2*i+1];
            cix[1] = i;
        }

        /* Negative oil phase flux */
        else {
            m [1]  =  mob[2*c+1];
            dm[1]  = dmob[2*c+1];
            cix[1] = c;
        }

        /* Positive water phase flux */
        if (darcyflux+m[1]*gravityflux>0) {
            m [0]  =  mob[2*i+0];
            dm[0]  = dmob[2*i+0];
            cix[0] = i;
        }

        /* Negative water phase flux */
        else {
            m [0]  =  mob[2*c+0];
            dm[0]  = dmob[2*c+0];
            cix[0] = c;
        }
    }
}

void 
spu_implicit_assemble(grid_t *g, double *s0, double *s, double *mob, double *dmob,
                      double *dflux, double *gflux, double *src, double dt, sparse_t *S, 
                      double *b, char *work)
{
    int     i, k, f, c1, c2, c;
    int     nc   = g->number_of_cells;

    double *m   = (double*)work;
    double *dm  =  m + 2;
    int    *cix = (int*)(dm + 2);
    double  m1,  m2, dm1, dm2, mt2;
    
    int    *pja;
    double *psa;
    double *d;
    double  sgn;
    double  df, gf;

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
            c   = (i==c1 ? c2 : c1);
            sgn = (i==c1)*2.0 - 1.0;
            df  = sgn * dt*dflux[f];
            gf  = sgn * dt*gflux[f];

            phase_upwind_mobility(df, gf,  i, c, mob, dmob, m, dm, cix);

            /* Ensure we do not divide by zero. */
            if (m[0] + m[1] > 0.0) {

                b[i]  += m[0]/(m[0]+m[1])*(df + m[1]*gf);

                mt2  = (m[0]+m[1])*(m[0]+m[1]);
                *psa = 0.0;
                *pja = c;
                
                /* dFw/dmw·dmw/dsw */
                if (cix[0] == c ) { *psa +=  m[1]/mt2*(df + m[1]*gf)*dm[0]; }
                else              { *d   +=  m[1]/mt2*(df + m[1]*gf)*dm[0]; }
                
                /* dFw/dmo·dmo/dsw */
                if (cix[1] == c) { *psa += -m[0]/mt2*(df - m[0]*gf)*dm[1];  }
                else             { *d   += -m[0]/mt2*(df - m[0]*gf)*dm[1];  }
                
                
                if (cix[0] == c || cix[1] == c) {
                    ++psa;
                    ++pja;
                }
            }            
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
