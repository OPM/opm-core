#ifndef COARSE_SYS_H_INCLUDED
#define COARSE_SYS_H_INCLUDED

/* ---------------------------------------------------------------------- */

struct coarse_sys {
    int *num_bf;                /* Number of bf's per block */
    int *bf_pos;                /* Start pointers to each block's bf's */
    int *ip_pos;                /* Start pointers to each block's IP */

    double *basis;              /* All basis functions */
    double *cell_ip;            /* Fine-scale IP contributions */
    double *Binv;               /* Coarse-scale inverse IP per block */
};

#endif  /* COARSE_SYS_H_INCLUDED */
