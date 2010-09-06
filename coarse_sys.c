#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "blas_lapack.h"
#include "coarse_conn.h"
#include "coarse_sys.h"


/* ---------------------------------------------------------------------- */
static int
max_element(int n, int *a)
/* ---------------------------------------------------------------------- */
{
    int i, ret;

    assert ((n > 0) && (a != NULL));

    ret = a[0];
    for (i = 1; i < n; i++) {
        ret = (a[i] > ret) ? a[i] : ret;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
void
coarse_sys_compute_cell_ip(int                     nc,
                           int                     max_nconn,
                           const int              *pconn,
                           const double           *Binv,
                           const int              *b2c_pos,
                           const int              *b2c,
                           struct coarse_topology *ct,
                           struct coarse_sys      *sys)
/* ---------------------------------------------------------------------- */
{
    int i, i1, i2, b, c, n, bf, *pconn2;
    int max_nbf, loc_nc;

    size_t p, nbf_pairs, bf_off, bf_sz;

    MAT_SIZE_T mm, nn, kk, nrhs, ld1, ld2, ld3, info;

    double a1, a2;
    double *work, *BI, *Psi, *IP;

    max_nbf = max_element(ct->nblocks, sys->num_bf);

    pconn2  = malloc((nc + 1) * sizeof *pconn2);
    work    = malloc(((max_nconn * max_nconn) + /* BI */
                      (max_nconn * max_nbf  ) + /* Psi */
                      (max_nbf   * max_nbf  ))  /* Ip */
                     * sizeof *work);

    if ((pconn2 != NULL) && (work != NULL)) {
        BI  = work + 0                      ;
        Psi = BI   + (max_nconn * max_nconn);
        IP  = Psi  + (max_nconn * max_nbf  );

        pconn2[0] = 0;

        for (i = 1; i <= nc; i++) {
            n         = pconn[i] - pconn[i - 1];
            pconn2[i] = pconn2[i - 1] + (n * n);
        }

        for (b = 0; b < ct->nblocks; b++) {
            loc_nc = b2c_pos[b + 1] - b2c_pos[b];
            bf_off = 0;

            assert ((sys->bf_pos[b + 1] -
                     sys->bf_pos[b + 0]) % sys->num_bf[b] == 0);

            bf_sz  = (sys->bf_pos[b + 1] - sys->bf_pos[b]) / sys->num_bf[b];

            nbf_pairs = sys->num_bf[b] * (sys->num_bf[b] + 1) / 2;

            for (i = 0; i < loc_nc; i++) {
                c = b2c[b2c_pos[b] + i];
                n = pconn[c + 1] - pconn[c];

                /* Linearise (collect) BF values for cell */
                p = sys->bf_pos[b] + bf_off;
                for (bf = 0; bf < sys->num_bf[b]; bf++, p += bf_sz) {
                    memcpy(Psi + bf*n, &sys->basis[p], n * sizeof *Psi);
                }

                /* Extract cell's inv(B) values... */
                memcpy(BI, &Binv[pconn2[c]], n * n * sizeof *BI);

                /* ...and (Cholesky) factor it... */
                nn  = n;
                ld1 = n;
                dpotrf_("Upper Triangular", &nn, BI, &ld1, &info);

                /* ...and solve BI*X = Psi (=> Psi (=X) <- B*Psi) */
                nrhs = sys->num_bf[b];
                ld2  = n;
                dpotrs_("Upper Triangular", &nn, &nrhs, BI, &ld1,
                        Psi, &ld2, &info);

                /* Finally, compute IP = Psi'*X = Psi'*B*Psi... */
                mm  = nn = sys->num_bf[b];
                kk  = n;
                ld1 = bf_sz;   ld2 = n;    ld3 = sys->num_bf[b];
                a1  = 1.0;     a2  = 0.0;
                dgemm_("Transpose", "No Transpose", &mm, &nn, &kk,
                       &a1, &sys->basis[sys->bf_pos[b] + bf_off], &ld1,
                       Psi, &ld2, &a2, IP, &ld3);

                /* ...and fill results into ip-contrib for this cell... */
                p = sys->ip_pos[b] + i*nbf_pairs;
                for (i2 = 0; i2 < sys->num_bf[b]; i2++) {
                    for (i1 = 0; i1 <= i2; i1++, p++) {
                        sys->cell_ip[p] = IP[i1 + i2*n];
                    }
                }

                /* ...and prepare for next cell. */
                bf_off += n;
            }
        }
    }

    free(work);  free(pconn2);
}


/* ---------------------------------------------------------------------- */
void
coarse_sys_compute_Binv(int                     max_bcells,
                        const double           *totmob,
                        const int              *b2c_pos,
                        const int              *b2c,
                        struct coarse_topology *ct,
                        struct coarse_sys      *sys,
                        double                 *work)
/* ---------------------------------------------------------------------- */
{
    int    b, i, i1, i2, nbf_pairs, loc_nc, nbf;

    double a1, a2;
    double *Lti, *B;

    size_t p1, p2;

    MAT_SIZE_T mm, nn, ld, incx, incy, info;

    Lti = work + 0;
    B   = work + max_bcells;

    incx = incy = 1;
    p2   = 0;
    for (b = 0; b < ct->nblocks; b++) {
        loc_nc = b2c_pos[b + 1] - b2c_pos[b];

        for (i = 0; i < loc_nc; i++) {
            Lti[i] = 1.0 / totmob[b2c[b2c_pos[b] + i]];
        }

        /* Form coarse inner-product matrix for block 'b' as (inverse)
         * mobility weighted sum of cell contributions */
        nbf       = sys->num_bf[b];
        nbf_pairs = nbf * (nbf + 1) / 2;

        mm = ld = nbf_pairs;
        nn = loc_nc;
        a1 = 1.0;
        a2 = 0.0;
        dgemv_("No Transpose", &mm, &nn, &a1,
               &sys->cell_ip[sys->ip_pos[b]], &ld,
               Lti, &incx, &a2, B, &incy);

        /* Factor (packed) SPD inner-product matrix... */
        dpptrf_("Upper Triangular", &mm, B, &info);

        /* ...and invert it... */
        dpptri_("Upper Triangular", &mm, B, &info);

        /* ...and write result to permanent storage suitable for
         * reduction functions hybsys_schur_comp*() */
        p1 = 0;
        for (i2 = 0; i2 < nbf; i2++) {
            for (i1 = 0; i1 <= i2; i1++, p1++) {
                sys->Binv[p2 + i1 + i2*nbf] = B[p1]; /* col i1 */
                sys->Binv[p2 + i2 + i1*nbf] = B[p1]; /* row i1 */
            }
        }

        p2 += nbf * nbf;
    }
}
