/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#ifndef COARSE_SYS_H_INCLUDED
#define COARSE_SYS_H_INCLUDED

#include "grid.h"

/* ---------------------------------------------------------------------- */

struct coarse_sys {
    int *blkdof_pos;         /* Start pointers to each block's dofs */
    int *basis_pos;          /* Start pointers to each block's bf's */
    int *cell_ip_pos;        /* Start pointers to each block's IP */

    int    *blkdof;          /* Each block's dofs */
    double *basis;           /* All basis functions */
    double *cell_ip;         /* Fine-scale IP contributions */
    double *Binv;            /* Coarse-scale inverse IP per block */
};


/* ---------------------------------------------------------------------- */

struct CSRMatrix;
typedef void (*LocalSolver)(struct CSRMatrix *A,
                            double           *b,
                            double           *x);

struct coarse_sys *
coarse_sys_construct(grid_t *g, const int   *p,
                     struct coarse_topology *ct,
                     const double           *perm,
                     const double           *src,
                     const double           *totmob,
                     LocalSolver             linsolve);

void
coarse_sys_destroy(struct coarse_sys *sys);

void
coarse_sys_compute_cell_ip(int                nc,
                           int                max_nconn,
                           int                nb,
                           const int         *pconn,
                           const double      *Binv,
                           const int         *b2c_pos,
                           const int         *b2c,
                           struct coarse_sys *sys);

void
coarse_sys_compute_Binv(int                nb,
                        int                max_bcells,
                        const double      *totmob,
                        const int         *b2c_pos,
                        const int         *b2c,
                        struct coarse_sys *sys,
                        double            *work);

#endif  /* COARSE_SYS_H_INCLUDED */
