/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/transport/reorder/tarjan.h>
#include <opm/core/transport/reorder/twophase.h>

#include <opm/core/transport/reorder/twophasetransport.h>

void twophasetransport(
    const double *porevolume,
    const double *source,
    double dt,
    struct UnstructuredGrid *grid, 
    const double *darcyflux,
    double *saturation) 
{
    int    i;

    int    ncomponents;
    int    *sequence;
    int    *components;

    struct cdata cd;
    struct vdata vd;

    /* Compute sequence of single-cell problems */
    sequence   = malloc(  grid->number_of_cells    * sizeof *sequence);
    components = malloc(( grid->number_of_cells+1) * sizeof *components);

    compute_sequence(grid, darcyflux, sequence, components, &ncomponents);
    assert(ncomponents == grid->number_of_cells);


    vd.saturation     = saturation;
    vd.fractionalflow = malloc(grid->number_of_cells * 
                               sizeof *vd.fractionalflow);
    for(i=0; i<grid->number_of_cells; ++i) 
    {
        vd.fractionalflow[i] = 0.0;
    }

    cd.grid       = grid;
    cd.darcyflux  = darcyflux;
    cd.source     = source;
    cd.porevolume = porevolume;
    cd.dt         = dt;

    /* Assume all strong components are single-cell domains. */
    for(i=0; i<grid->number_of_cells; ++i) 
    { 
        solve(&vd, &cd, sequence[i]);
    }

    free(vd.fractionalflow);
    free(sequence);
    free(components);
}
