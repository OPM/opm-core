/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <opm/core/grid.h>
#include <opm/core/transport/reorder/twophase.h>
#include <opm/core/transport/reorder/nlsolvers.h>

static struct Parameters get_parameters(struct vdata *vd, 
                                        const struct cdata *cd, 
                                        int cell);
static double fluxfun(double s);

static double 
fluxfun(double s)
{
    return s*s/(s*s + (1-s)*(1-s));
}

struct Parameters
{
    double s0;
    double influx;  /* sum_j min(v_ij, 0)*f(s_j) */
    double outflux; /* sum_j max(v_ij, 0)        */
    double dtpv;    /* dt/pv(i)                  */
};
#include <stdio.h>
static double 
G(double s, void *data)
{
    struct Parameters *p = data;
    fprintf(stderr, "s:%f dtpv:%f s0:%f, influx:%f outflux:%f G:%f\n", 
            s, 
            p->dtpv,
            p->s0, 
            p->influx,
            p->outflux,
            s - p->s0 +  p->dtpv*(p->outflux*fluxfun(s) + p->influx));
    /* G(s) = s - s0 + dt/pv*(influx - outflux*f(s) ) */
    return s - p->s0 +  p->dtpv*(p->outflux*fluxfun(s) + p->influx);
}

static struct Parameters 
get_parameters(struct vdata *vd, const struct cdata *cd, int cell)
{
    struct UnstructuredGrid *g  = cd->grid;
    struct Parameters        p;

    p.s0      = vd->saturation[cell];
    p.influx  = cd->source[cell] >  0 ? -cd->source[cell] : 0.0;
    p.outflux = cd->source[cell] <= 0 ? -cd->source[cell] : 0.0;
    p.dtpv    = cd->dt/cd->porevolume[cell];
    
    fprintf(stderr, "src:%f pv:%f\n",
            cd->source[cell],
            cd->porevolume[cell]);

    int i;
    vd->saturation[cell] = 0;
    for (i=g->cell_facepos[cell]; i<g->cell_facepos[cell+1]; ++i) {
        int    f = g->cell_faces[i];
        double flux;
        int    other;

        /* Compute cell flux*/
        if (cell == g->face_cells[2*f]) {
            flux  = cd->darcyflux[f];
            other = g->face_cells[2*f+1];
        }
        else {
            flux  =-cd->darcyflux[f];
            other = g->face_cells[2*f];
        }
       
        fprintf(stderr, "flux:%f\n", flux);
        if (other != -1) {
            if (flux < 0.0) {
                p.influx  += flux*vd->fractionalflow[other];
            }
            else {
                p.outflux += flux;
            }
        }
    }
    return p;
}


static enum Method {RIDDER, REGULAFALSI, BISECTION} method = REGULAFALSI;

void solve(void *vdata, const void *cdata, int cell)
{
    const struct cdata *cd = cdata;
    struct vdata       *vd = vdata;

    struct Parameters   p  = get_parameters(vdata, cdata, cell);

    int    iterations;
    int    maxit = 20;
    double TOL = 1e-9;

    switch (method) {
    default:
    case BISECTION:
	vd->saturation[cell] = bisection  ( &G, &p, TOL, maxit, &iterations);
        break;

    case RIDDER:
	vd->saturation[cell] = ridder     ( &G, &p, TOL, maxit, &iterations); 
        break;

    case REGULAFALSI:
	vd->saturation[cell] = regulafalsi( &G, &p, TOL, maxit, &iterations); 
        break;
    }
    
    vd->fractionalflow[cell]   = fluxfun(vd->saturation[cell]);
}



/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */

