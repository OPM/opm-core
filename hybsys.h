/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#ifndef HYBSYS_H_INCLUDED
#define HYBSYS_H_INCLUDED

struct hybsys {
    double *L;                  /* C2' * inv(B) * C1 - P */
    double *q;                  /* g - F2*G */
    double *F1;                 /* C1' * inv(B)          */
    double *F2;                 /* C2' * inv(B)          */
    double *r;                  /* system rhs in single cell */
    double *S;                  /* system matrix in single cell */
    double *one;                /* ones(max_nconn, 1) */
};


struct hybsys_well {
    double *F1;
    double *F2;
    double *r;

    double *w2r;
    double *r2w;
    double *w2w;

    double *data;
};


struct hybsys *
hybsys_allocate_symm(int max_nconn, int nc, int nconn_tot);

struct hybsys *
hybsys_allocate_unsymm(int max_nconn, int nc, int nconn_tot);

struct hybsys_well *
hybsys_well_allocate_symm(int max_nconn, int nc, int *cwpos);

struct hybsys_well *
hybsys_well_allocate_unsymm(int max_nconn, int nc, int *cwpos);

void
hybsys_free(struct hybsys *sys);

void
hybsys_well_free(struct hybsys_well *wsys);

void
hybsys_init(int max_nconn, struct hybsys *sys);

/*
 * Schur complement reduction (per grid cell) of block matrix
 *
 *    [  B   C   D  ]
 *    [  C'  0   0  ]
 *    [  D'  0   0  ]
 *
 */
void
hybsys_schur_comp_symm(int nc, const int *pconn,
                       const double *Binv, struct hybsys *sys);

/*
 * Schur complement reduction (per grid cell) of block matrix
 *
 *    [    B     C   D  ]
 *    [  (C-V)'  P   0  ]
 *    [    D'    0   0  ]
 *
 */
void
hybsys_schur_comp_unsymm(int nc, const int *pconn,
                         const double *Binv, const double *BIV,
                         const double *P, struct hybsys *sys);

/*
 * Schur complement reduction (per grid cell) of block matrix
 *
 *    [   B   C   D  ]
 *    [  C2'  P   0  ]
 *    [   D'  0   0  ]
 *
 */
void
hybsys_schur_comp_gen(int nc, const int *pconn,
                      const double *Binv, const double *C2,
                      const double *P, struct hybsys *sys);

void
hybsys_well_schur_comp_symm(int nc, const int *cwpos,
                            double             *WI,
                            struct hybsys      *sys,
                            struct hybsys_well *wsys);

void
hybsys_cellcontrib_symm(int c, int nconn, int p1, int p2,
                        const double *gpress, const double *src,
                        const double *Binv, struct hybsys *sys);

void
hybsys_cellcontrib_unsymm(int c, int nconn, int p1, int p2,
                          const double *gpress, const double *src,
                          const double *Binv, struct hybsys *sys);

void
hybsys_well_cellcontrib_symm(int c, int ngconn, int p1,
                             const int *cwpos, const int *cwells,
                             const double *WI, const double *wdp,
                             struct hybsys *sys, struct hybsys_well *wsys);

void
hybsys_compute_press_flux(int nc, const int *pconn, const int *conn,
                          const double *gpress, const double *src,
                          const double *Binv, const struct hybsys *sys,
                          const double *pi, double *press, double *flux,
                          double *work);

#endif  /* HYBSYS_H_INCLUDED */
