#ifndef HYBSYS_GLOBAL_H_INCLUDED
#define HYBSYS_GLOBAL_H_INCLUDED

#include "grid.h"
#include "well.h"
#include "sparse_sys.h"

struct CSRMatrix *
hybsys_define_globconn(grid_t *G, well_t *W);


void
hybsys_global_assemble_cell(int nconn, int *l2g,
                            const double     *S,
                            const double     *r,
                            struct CSRMatrix *A,
                            double           *b);

void
hybsys_global_assemble_well_sym(int ngconn_tot,
                                int ngconn, const int *gconn,
                                int nwconn, const int *wconn,
                                const double     *r2w,
                                const double     *w2w,
                                const double     *r,
                                struct CSRMatrix *A,
                                double           *b);

#endif  /* HYBSYS_GLOBAL_H_INCLUDED */
