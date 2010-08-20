/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#ifndef HYBSYS_H_INCLUDED
#define HYBSYS_H_INCLUDED

#ifndef MAT_SIZE_T
#define MAT_SIZE_T int
#endif

struct hybsys {
    double *L;                  /* C2' * inv(B) * C1 - P */
    double *F1;                 /* C1' * inv(B)          */
    double *F2;                 /* C2' * inv(B)          */
    double *r;                  /* system rhs in single cell */
    double *S;                  /* system matrix in single cell */
    double *one;                /* ones(max_nconn, 1) */
};



struct Sparse
{
    int        m;
    int        n;
    MAT_SIZE_T *ia;
    MAT_SIZE_T *ja;
    double     *sa;
};



struct hybsys *
hybsys_allocate_symm(int max_nconn, int nc, int nconn_tot);

struct hybsys *
hybsys_allocate_unsymm(int max_nconn, int nc, int nconn_tot);

void
hybsys_free(struct hybsys *sys);

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
hybsys_cellcontrib_symm(int c, int nconn, int p1, int p2,
                        const double *gpress, const double *src,
                        const double *Binv, struct hybsys *sys);

void
hybsys_cellcontrib_unsymm(int c, int nconn, int p1, int p2,
                          const double *gpress, const double *src,
                          const double *Binv, struct hybsys *sys);

void
hybsys_compute_press_flux(int nc, const int *pconn, const int *conn,
                          const double *gpress, const double *src,
                          const double *Binv, const struct hybsys *sys,
                          const double *pi, double *press, double *flux,
                          double *work);

void
hybsys_assemble(int nc, int nf,
                const int *pconn, const int *conn,
                const double *S, const double *R,
                struct Sparse *A, double **b);

#endif  /* HYBSYS_H_INCLUDED */
