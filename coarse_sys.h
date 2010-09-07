#ifndef COARSE_SYS_H_INCLUDED
#define COARSE_SYS_H_INCLUDED

/* ---------------------------------------------------------------------- */

struct coarse_sys {
    int *dof_pos;              /* Start pointers to each block's dofs */
    int *bf_pos;               /* Start pointers to each block's bf's */
    int *ip_pos;               /* Start pointers to each block's IP */

    int    *blkdof;            /* Each block's dofs */
    double *basis;             /* All basis functions */
    double *cell_ip;           /* Fine-scale IP contributions */
    double *Binv;              /* Coarse-scale inverse IP per block */
};


/* ---------------------------------------------------------------------- */

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
