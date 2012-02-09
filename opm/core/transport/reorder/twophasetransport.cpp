/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/transport/reorder/tarjan.h>
#include <opm/core/transport/reorder/nlsolvers.h>
#include <opm/core/transport/reorder/twophase.hpp>
#include <opm/core/transport/reorder/twophasetransport.hpp>
#include <opm/core/transport/reorder/TransportModelTwophase.hpp>


void twophasetransport(
    const double            *porevolume,
    const double            *source,
    double                   dt,
    struct UnstructuredGrid *grid,
    const Opm::IncompPropertiesInterface* props,
    const double            *darcyflux,
    double                  *saturation)
{
    int    i;

    int    ncomponents;
    int    *sequence;
    int    *components;

    struct SolverData *data = init_solverdata(grid, props, darcyflux,
                                              porevolume, source, dt, saturation);


    struct NonlinearSolverCtrl ctrl;


    /* Compute sequence of single-cell problems */
    sequence   = (int*) malloc(  grid->number_of_cells    * sizeof *sequence);
    components = (int*) malloc(( grid->number_of_cells+1) * sizeof *components);

    compute_sequence(grid, darcyflux, sequence, components, &ncomponents);
    assert(ncomponents == grid->number_of_cells);



    ctrl.maxiterations = 20;
    ctrl.nltolerance   = 1e-9;
    ctrl.method        = NonlinearSolverCtrl::REGULAFALSI;
    ctrl.initialguess  = 0.5;
    ctrl.min_bracket   = 0.0;
    ctrl.max_bracket   = 1.0;

    /* Assume all strong components are single-cell domains. */
    for(i=0; i<grid->number_of_cells; ++i)
    {
#ifdef MATLAB_MEX_FILE
        if (interrupt_signal)
        {
            mexPrintf("Reorder loop interrupted by user: %d of %d "
                      "cells finished.\n", i, grid->number_of_cells);
            break;
        }
#endif
        solvecell(data, &ctrl, sequence[i]);
        /*print("iterations:%d residual:%20.16f\n",
         *  ctrl.iterations, ctrl.residual);*/
    }

    destroy_solverdata(data);

    free(sequence);
    free(components);
}


void Opm::reorderTransportTwophase(const double            *porevolume,
				   const double            *source,
				   const double             dt,
				   const UnstructuredGrid  *grid,
				   const IncompPropertiesInterface* props,
				   const double            *darcyflux,
				   double                  *saturation)
{
    TransportModelTwophase tmodel(grid, props, darcyflux,
				  porevolume, source, dt, saturation);

    // Compute sequence of single-cell problems
    std::vector<int> sequence(grid->number_of_cells);
    std::vector<int> components(grid->number_of_cells + 1);
    int ncomponents;
    compute_sequence(grid, darcyflux, &sequence[0], &components[0], &ncomponents);
    assert(ncomponents == grid->number_of_cells);


    // Assume all strong components are single-cell domains.
    for (int i = 0; i < grid->number_of_cells; ++i)
    {
#ifdef MATLAB_MEX_FILE
        if (interrupt_signal)
        {
            mexPrintf("Reorder loop interrupted by user: %d of %d "
                      "cells finished.\n", i, grid->number_of_cells);
            break;
        }
#endif
        tmodel.solveSingleCell(sequence[i]);
        // print("iterations:%d residual:%20.16f\n",
        //  ctrl.iterations, ctrl.residual);
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
