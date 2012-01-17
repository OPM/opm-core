/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef TWOPHASETRANSPORT_H_INCLUDED
#define TWOPHASETRANSPORT_H_INCLUDED

struct UnstructuredGrid;
int twophasetransport(
    const double *porevolume,
    const double *source,
    double dt,
    struct UnstructuredGrid *grid, 
    const double *darcyflux,
    double *saturation);

#endif /* TWOPHASETRANSPORT_H_INCLUDED */

