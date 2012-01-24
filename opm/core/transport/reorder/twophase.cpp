/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdlib.h>

#include <opm/core/grid.h>
#include <opm/core/transport/reorder/twophase.hpp>
#include <opm/core/transport/reorder/nlsolvers.h>
#include <opm/core/transport/reorder/fluid.h>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>


/* Parameters used in solution of single-cell boundary-value problem */
struct Parameters
{
    double s0;
    double influx;  /* sum_j min(v_ij, 0)*f(s_j) */
    double outflux; /* sum_j max(v_ij, 0)        */
    double dtpv;    /* dt/pv(i)                  */
};


static struct Parameters get_parameters(struct SolverData *d, int cell);
static double residual(double s, void *data);
/* static double fluxfun(double s, int cell); */

void
destroy_solverdata(struct SolverData *d)
{
    if (d!=NULL)
    {
        free(d->fractionalflow);
    }
    free(d);
}

struct SolverData *
init_solverdata(struct UnstructuredGrid *grid,
		const Opm::IncompPropertiesInterface* props,
		const double *darcyflux,
                const double *porevolume,
		const double *source,
                const double dt,
		double *saturation)
{
    int i;
    struct SolverData *d = (struct SolverData*) malloc(sizeof *d);

    if(d!=NULL)
    {
        d->grid       = grid;
	d->props      = props;
        d->darcyflux  = darcyflux;
        d->porevolume = porevolume;
        d->source     = source;
        d->dt         = dt;

        d->saturation     = saturation;
        d->fractionalflow = (double*) malloc(grid->number_of_cells *
					     sizeof *d->fractionalflow);
        if (d->fractionalflow == NULL)
        {
            destroy_solverdata(d);
            d = NULL;
        }
        for(i=0; i<grid->number_of_cells; ++i)
        {
            d->fractionalflow[i] = 0.0;
        }
    }
    return d;
}

/* Solver for single-cell bvp calls root-finder in nlsolvers.c */
void solvecell(void *data, struct NonlinearSolverCtrl *ctrl, int cell)
{
    struct SolverData   *d   = (struct SolverData*) data;
    struct Parameters   prm = get_parameters(d, cell);
    
    d->saturation[cell] = find_zero(residual, &prm, ctrl);
    d->fractionalflow[cell] = fluxfun(d->saturation[cell], -999);
}


/* ====================== Internals =================================*/


/* static double  */
/* fluxfun(double s, int region) */
/* { */
/*     return s*s/(s*s + (1-s)*(1-s)); */
/* } */

/* Residual function r(s) for a single-cell bvp */
/*
 *     r(s) = s - s0 + dt/pv*(influx - outflux*f(s) )
 */
/* influx is water influx, outflux is total outflux */
static double
residual(double s, void *data)
{
    struct Parameters *p = (struct Parameters*) data;
    return s - p->s0 +  p->dtpv*(p->outflux*fluxfun(s, -999) + p->influx);
}

static struct Parameters
get_parameters(struct SolverData *d, int cell)
{
    int i;
    struct UnstructuredGrid *g  = d->grid;
    struct Parameters        p;
    double flux;
    int f, other;

    p.s0      = d->saturation[cell];
    p.influx  = d->source[cell] >  0 ? -d->source[cell] : 0.0;
    p.outflux = d->source[cell] <= 0 ? -d->source[cell] : 0.0;
    p.dtpv    = d->dt/d->porevolume[cell];

    d->saturation[cell] = 0;
    for (i=g->cell_facepos[cell]; i<g->cell_facepos[cell+1]; ++i) {
        f = g->cell_faces[i];

        /* Compute cell flux*/
        if (cell == g->face_cells[2*f]) {
            flux  = d->darcyflux[f];
            other = g->face_cells[2*f+1];
        }
        else {
            flux  =-d->darcyflux[f];
            other = g->face_cells[2*f];
        }

        if (other != -1) {
            if (flux < 0.0) {
                p.influx  += flux*d->fractionalflow[other];
            }
            else {
                p.outflux += flux;
            }
        }
    }
    return p;
}



/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
