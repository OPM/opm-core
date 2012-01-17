/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef TWOPHASE_H_INCLUDED
#define TWOPHASE_H_INCLUDED


struct vdata {
    double *saturation;      /* one per cell */
    double *fractionalflow;  /* one per cell */
};

struct UnstructuredGrid;
struct cdata  {
    struct UnstructuredGrid *grid;
    const double *darcyflux;   /* one flux per face  in cdata::grid*/
    const double *source;      /* one source per cell */
    const double *porevolume;  /* one volume per cell */
    double  dt;
};

void solve(void *vdata, const void *cdata, int cell);
#endif /* TWOPHASE_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
