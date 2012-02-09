/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <opm/core/transport/reorder/twophasetransport.hpp>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/transport/reorder/TransportModelTwophase.hpp>
#include <opm/core/grid.h>

#include <cassert>



void Opm::reorderTransportTwophase(const double            *porevolume,
				   const double            *source,
				   const double             dt,
				   const UnstructuredGrid  *grid,
				   const IncompPropertiesInterface* props,
				   const double            *darcyflux,
				   double                  *saturation)
{
    // Set up transport model.
    TransportModelTwophase tmodel(grid, props, darcyflux,
				  porevolume, source, dt, saturation);

    // Compute sequence of single-cell problems
    std::vector<int> sequence(grid->number_of_cells);
    std::vector<int> components(grid->number_of_cells + 1);
    int ncomponents;
    compute_sequence(grid, darcyflux, &sequence[0], &components[0], &ncomponents);

    // Assume all strong components are single-cell domains.
    assert(ncomponents == grid->number_of_cells);
    for (int i = 0; i < grid->number_of_cells; ++i) {
#ifdef MATLAB_MEX_FILE
        if (interrupt_signal) {
            mexPrintf("Reorder loop interrupted by user: %d of %d "
                      "cells finished.\n", i, grid->number_of_cells);
            break;
        }
#endif
        tmodel.solveSingleCell(sequence[i]);
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
