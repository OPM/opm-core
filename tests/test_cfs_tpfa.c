#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <suitesparse/umfpack.h>

#include <opmcore/grid.h>
#include <opmcore/well.h>

#include <opmcore/linalg/sparse_sys.h>

#include <opmcore/pressure/flow_bc.h>

#include <opmcore/pressure/tpfa/cfs_tpfa.h>
#include <opmcore/pressure/tpfa/compr_quant.h>
#include <opmcore/pressure/tpfa/trans_tpfa.h>



struct CSCMatrix {
    UF_long  n;
    UF_long  nnz;

    UF_long *p;
    UF_long *i;
    double  *x;
};


/* ---------------------------------------------------------------------- */
static void
csc_deallocate(struct CSCMatrix *csc)
/* ---------------------------------------------------------------------- */
{
    if (csc != NULL) {
        if (csc->x != NULL) { free(csc->x); }
        if (csc->i != NULL) { free(csc->i); }
        if (csc->p != NULL) { free(csc->p); }

        free(csc);
    }
}


/* ---------------------------------------------------------------------- */
static struct CSCMatrix *
csc_allocate(UF_long n, UF_long nnz)
/* ---------------------------------------------------------------------- */
{
    struct CSCMatrix *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->p = malloc((n + 1) * sizeof *new->p);
        new->i = malloc(nnz     * sizeof *new->i);
        new->x = malloc(nnz     * sizeof *new->x);

        if ((new->p == NULL) || (new->i == NULL) || (new->x == NULL)) {
            csc_deallocate(new);
            new = NULL;
        } else {
            new->n   = n;
            new->nnz = nnz;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
csr_to_csc(const int        *ia,
           const int        *ja,
           const double     *sa,
           struct CSCMatrix *csc)
/* ---------------------------------------------------------------------- */
{
    UF_long i, nz;

    /* Clear garbage, prepare for counting */
    for (i = 0; i <= csc->n; i++) { csc->p[i] = 0; }

    /* Count column connections */
    for (nz = 0; nz < csc->nnz; nz++) {
        csc->p[ ja[nz] + 1 ] += 1;
    }

    /* Define column start pointers */
    for (i = 1; i <= csc->n; i++) {
        csc->p[0] += csc->p[i];
        csc->p[i]  = csc->p[0] - csc->p[i];
    }

    assert (csc->p[0] == csc->nnz);

    /* Fill matrix whilst defining column end pointers */
    for (i = nz = 0; i < csc->n; i++) {
        for (; nz < ia[i + 1]; nz++) {
            csc->i[ csc->p[ ja[nz] + 1 ] ] = i;      /* Insertion sort */
            csc->x[ csc->p[ ja[nz] + 1 ] ] = sa[nz]; /* Insert mat elem */

            csc->p        [ ja[nz] + 1 ]  += 1;      /* Advance col ptr */
        }
    }

    assert (csc->p[csc->n] == csc->nnz);

    csc->p[0] = 0;
}


/* ---------------------------------------------------------------------- */
static void
solve_umfpack(struct CSCMatrix *csc, const double *b, double *x)
/* ---------------------------------------------------------------------- */
{
    void *Symbolic, *Numeric;
    double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];

    umfpack_dl_defaults(Control);

    umfpack_dl_symbolic(csc->n, csc->n, csc->p, csc->i, csc->x,
                        &Symbolic, Control, Info);
    umfpack_dl_numeric (csc->p, csc->i, csc->x,
                        Symbolic, &Numeric, Control, Info);

    umfpack_dl_free_symbolic(&Symbolic);

    umfpack_dl_solve(UMFPACK_A, csc->p, csc->i, csc->x, x, b,
                     Numeric, Control, Info);

    umfpack_dl_free_numeric(&Numeric);
}


/*---------------------------------------------------------------------------*/
void
call_UMFPACK(struct CSRMatrix *A, double *b, double *x)
/*---------------------------------------------------------------------------*/
{
    struct CSCMatrix *csc;

    csc = csc_allocate(A->m, A->ia[A->m]);

    if (csc != NULL) {
        csr_to_csc(A->ia, A->ja, A->sa, csc);

        solve_umfpack(csc, b, x);
    }

    csc_deallocate(csc);
}

static void
deallocate_cart_grid(grid_t *G)
{
    if (G != NULL) {
        free(G->node_coordinates);

        free(G->face_nodes);
        free(G->face_nodepos);
        free(G->face_cells);
        free(G->face_centroids);
        free(G->face_normals);
        free(G->face_areas);

        free(G->cell_faces);
        free(G->cell_facepos);
        free(G->cell_centroids);
        free(G->cell_volumes);
    }

    free(G);
}

grid_t *
cart_grid(int nx, int ny, int nz)
{
    /* The following applies to any grid derived from base_grid_t. */
    int i,j,k;
    int Nx = nx+1;
    int Ny = ny+1;
    int Nz = nz+1;

    grid_t *G = malloc(1 * sizeof *G);

    G->dimensions = 3;

    int nxf = Nx*ny*nz;
    int nyf = nx*Ny*nz;
    int nzf = nx*ny*Nz;

    G->number_of_cells  = nx*ny*nz;
    G->number_of_faces  = nxf+nyf+nzf;
    G->number_of_nodes  = Nx*Ny*Nz;

    G->node_coordinates = malloc(G->number_of_nodes * 3 * sizeof *(G->node_coordinates));

    G->face_nodes       = malloc(G->number_of_faces * 4 * sizeof *(G->face_nodes));
    G->face_nodepos     = malloc((G->number_of_faces+1) * sizeof *(G->face_nodepos));
    G->face_cells       = malloc(G->number_of_faces * 2 * sizeof *(G->face_cells));
    G->face_centroids   = malloc(G->number_of_faces * 3 * sizeof *(G->face_centroids));
    G->face_normals     = malloc(G->number_of_faces * 3 * sizeof *(G->face_normals));
    G->face_areas       = malloc(G->number_of_faces * 1 * sizeof *(G->face_areas));

    G->cell_faces       = malloc(G->number_of_cells * 6 * sizeof *(G->cell_faces));
    G->cell_facepos     = malloc((G->number_of_cells+1) * sizeof *(G->cell_facepos));
    G->cell_centroids   = malloc(G->number_of_cells * 3 * sizeof *(G->cell_centroids));
    G->cell_volumes     = malloc(G->number_of_cells * 1 * sizeof *(G->cell_volumes));


    int    *cfaces     = G->cell_faces;
    int    *cfacepos   = G->cell_facepos;
    double *ccentroids = G->cell_centroids;
    double *cvolumes   = G->cell_volumes;
    for (k=0; k<nz; ++k)  {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *cfaces++ = i+  Nx*(j+  ny* k   );
                *cfaces++ = i+1+Nx*(j+  ny* k   );
                *cfaces++ = i+  nx*(j+  Ny* k   )  +nxf;
                *cfaces++ = i+  nx*(j+1+Ny* k   )  +nxf;
                *cfaces++ = i+  nx*(j+  ny* k   )  +nxf+nyf;
                *cfaces++ = i+  nx*(j+  ny*(k+1))  +nxf+nyf;

                cfacepos[1] = cfacepos[0]+6;
                ++cfacepos;

                *ccentroids++ = i+0.5;
                *ccentroids++ = j+0.5;
                *ccentroids++ = k+0.5;

                *cvolumes++ = 1;
            }
        }
    }


    int    *fnodes     = G->face_nodes;
    int    *fnodepos   = G->face_nodepos;
    int    *fcells     = G->face_cells;
    double *fnormals   = G->face_normals;
    double *fcentroids = G->face_centroids;
    double *fareas     = G->face_areas;

    /* Faces with x-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx+1; ++i) {
                *fnodes++ = i+Nx*(j   + Ny * k   );
                *fnodes++ = i+Nx*(j+1 + Ny * k   );
                *fnodes++ = i+Nx*(j+1 + Ny *(k+1));
                *fnodes++ = i+Nx*(j   + Ny *(k+1));
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (i==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (i == nx) {
                    *fcells++ =  i-1+nx*(j+ny*k);
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i-1 + nx*(j+ny*k);
                    *fcells++ =  i   + nx*(j+ny*k);
                }

                *fnormals++ = 1;
                *fnormals++ = 0;
                *fnormals++ = 0;

                *fcentroids++ = i;
                *fcentroids++ = j+0.5;
                *fcentroids++ = k+0.5;

                *fareas++ = 1;
            }
        }
    }
    /* Faces with y-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx; ++i) {
                *fnodes++ = i+    Nx*(j + Ny * k   );
                *fnodes++ = i   + Nx*(j + Ny *(k+1));
                *fnodes++ = i+1 + Nx*(j + Ny *(k+1));
                *fnodes++ = i+1 + Nx*(j + Ny * k   );
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (j==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (j == ny) {
                    *fcells++ =  i+nx*(j-1+ny*k);
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i+nx*(j-1+ny*k);
                    *fcells++ =  i+nx*(j+ny*k);
                }

                *fnormals++ = 0;
                *fnormals++ = 1;
                *fnormals++ = 0;

                *fcentroids++ = i+0.5;
                *fcentroids++ = j;
                *fcentroids++ = k+0.5;

                *fareas++ = 1;
            }
        }
    }
    /* Faces with z-normal */
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *fnodes++ = i+    Nx*(j   + Ny * k);
                *fnodes++ = i+1 + Nx*(j   + Ny * k);
                *fnodes++ = i+1 + Nx*(j+1 + Ny * k);
                *fnodes++ = i+    Nx*(j+1 + Ny * k);
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (k==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (k == nz) {
                    *fcells++ =  i+nx*(j+ny*(k-1));
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i+nx*(j+ny*(k-1));
                    *fcells++ =  i+nx*(j+ny*k);
                }

                *fnormals++ = 0;
                *fnormals++ = 0;
                *fnormals++ = 1;

                *fcentroids++ = i+0.5;
                *fcentroids++ = j+0.5;
                *fcentroids++ = k;

                *fareas++ = 1;
            }
        }
    }

    double *coord = G->node_coordinates;
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx+1; ++i) {
                *coord++ = i;
                *coord++ = j;
                *coord++ = k;
            }
        }
    }

    return G;
}


static void
deallocate_cq(struct compr_quantities *cq)
{
    if (cq != NULL) {
        free(cq->phasemobf);
        free(cq->Af);
        free(cq->Ac);
        free(cq->voldiscr);
        free(cq->totcompr);
    }

    free(cq);
}

static struct compr_quantities *
allocate_cq(size_t nc, size_t nf, int np)
{
    struct compr_quantities *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->totcompr  = malloc(nc           * sizeof *new->totcompr );
        new->voldiscr  = malloc(nc           * sizeof *new->voldiscr );
        new->Ac        = malloc(nc * np * np * sizeof *new->Ac       );
        new->Af        = malloc(nf * np * np * sizeof *new->Af       );
        new->phasemobf = malloc(nf * np      * sizeof *new->phasemobf);

        if ((new->totcompr == NULL) || (new->voldiscr == NULL) ||
            (new->Ac == NULL) || (new->Af == NULL) ||
            (new->phasemobf == NULL)) {
            deallocate_cq(new);
            new = NULL;
        } else {
            new->nphases = np;
        }
    }

    return new;
}

static void
vector_assign(size_t n, double value, double *v)
{
    size_t i;

    for (i = 0; i < n; i++) {
        v[i] = value;
    }
}

static void
vector_ones(size_t n, double *v)
{
    vector_assign(n, 1.0, v);
}

static void
set_incompressible(size_t nc, size_t nf, struct compr_quantities *cq)
{
    size_t i, j, np;

    np = cq->nphases;

    vector_zero(nc          , cq->totcompr );
    vector_zero(nc          , cq->voldiscr );
    vector_zero(nc * np * np, cq->Ac);
    vector_zero(nf * np * np, cq->Af);
    vector_ones(nf * np     , cq->phasemobf);

    for (i = 0; i < nc; i++) {
        for (j = 0; j < np; j++) {
            cq->Ac[i*np*np + j*(np + 1)] = 1.0;
        }
    }

    for (i = 0; i < nf; i++) {
        for (j = 0; j < np; j++) {
            cq->Af[i*np*np + j*(np + 1)] = 1.0;
        }
    }
}


static void
set_homoperm(size_t nc, size_t dim, double *perm)
{
    size_t i, j;

    vector_zero(nc * dim * dim, perm);

    for (i = 0; i < nc; i++) {
        for (j = 0; j < dim; j++) {
            perm[i*dim*dim + j*(dim + 1)] = 1.0;
        }
    }
}


static void
destroy_wells(well_t *W)
{
    if (W != NULL) {
        free(W->well_cells);
        free(W->well_connpos);
    }

    free(W);
}


static well_t *
create_wells(grid_t *G)
{
    well_t *new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->well_connpos = malloc((2 + 1) * sizeof *new->well_connpos);
        new->well_cells   = malloc(2       * sizeof *new->well_cells  );

        if ((new->well_connpos == NULL) ||
            (new->well_cells   == NULL)) {
            destroy_wells(new);
            new = NULL;
        } else {
            new->number_of_wells = 2;

            new->well_connpos[0] = 0;
            new->well_connpos[1] = 1;
            new->well_connpos[2] = 2;

            new->well_cells[0] = 0;
            new->well_cells[1] = G->number_of_cells - 1;
        }
    }

    return new;
}


static void
destroy_well_control(well_control_t *ctrl)
{
    if (ctrl != NULL) {
        free(ctrl->target);
        free(ctrl->ctrl);
        free(ctrl->type);
    }

    free(ctrl);
}


static well_control_t *
create_well_control(void)
{
    well_control_t *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->type   = malloc(2 * sizeof *new->type  );
        new->ctrl   = malloc(2 * sizeof *new->ctrl  );
        new->target = malloc(2 * sizeof *new->target);

        if ((new->type == NULL) || (new->ctrl == NULL) ||
            (new->target == NULL)) {
            destroy_well_control(new);
            new = NULL;
        } else {
            new->type[0] = INJECTOR;  new->type[1] = PRODUCER;
            new->ctrl[0] = new->ctrl[1] = BHP;
            new->target[0] = 2;
            new->target[1] = 1;
        }
    }

    return new;
}


static void
destroy_completion_data(struct completion_data *cd)
{
    if (cd != NULL) {
        free(cd->WI);
        free(cd->gpot);
        free(cd->A);
        free(cd->phasemob);
    }

    free(cd);
}


static struct completion_data *
create_completion_data(well_t *W, int np)
{
    size_t i, j, totconn;
    struct completion_data *new;

    totconn = W->well_connpos[ W->number_of_wells ];
    new     = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->WI       = malloc(totconn           * sizeof *new->WI);
        new->gpot     = malloc(totconn * np      * sizeof *new->gpot);
        new->A        = malloc(totconn * np * np * sizeof *new->A);
        new->phasemob = malloc(totconn * np      * sizeof *new->phasemob);

        if ((new->WI      == NULL) ||
            (new->gpot    == NULL) ||
            (new->A       == NULL) ||
            (new->phasemob == NULL)) {
            destroy_completion_data(new);
            new = NULL;
        } else {
            vector_ones(totconn          , new->WI);
            vector_zero(totconn * np     , new->gpot);
            vector_zero(totconn * np * np, new->A);
            vector_ones(totconn * np     , new->phasemob);

            for (i = 0; i < totconn; i++) {
                for (j = 0; j < (size_t)np; j++) {
                    new->A[i*np*np + j*(np + 1)] = 1.0;
                }
            }
        }
    }

    return new;
}


int
main(void)
{
    int       i, nphases, dim;
    flowbc_t *bc;
    double    dt;
    double   *htrans, *trans, *perm, *totmob, *src;
    double   *gravcap_f, *porevol, *cpress0;
    double   *cpress, *fpress, *fflux, *wpress, *wflux;

    grid_t                 *G;
    well_t                 *W;
    well_control_t         *wctrl;
    struct completion_data *wdata;

    struct compr_quantities *cq;

    struct cfs_tpfa_data *h;

    nphases = 2;

    G = cart_grid(5, 5, 1);
    W = create_wells(G);
    wctrl = create_well_control();
    wdata = create_completion_data(W, nphases);

    dim    = G->dimensions;

    htrans = malloc(G->cell_facepos[ G->number_of_cells ] * sizeof *htrans);
    trans  = malloc(G->number_of_faces                    * sizeof *trans );
    perm   = malloc(dim * dim * G->number_of_cells        * sizeof *perm  );
    totmob = malloc(G->number_of_cells                    * sizeof *totmob);
    src    = malloc(G->number_of_cells                    * sizeof *src   );

    gravcap_f = malloc(G->number_of_faces * nphases * sizeof *gravcap_f);
    porevol   = malloc(G->number_of_cells           * sizeof *porevol  );
    cpress0   = malloc(G->number_of_cells           * sizeof *cpress0  );
    cpress    = malloc(G->number_of_cells           * sizeof *cpress   );
    fpress    = malloc(G->number_of_faces           * sizeof *fpress   );
    fflux     = malloc(G->number_of_faces           * sizeof *fflux    );

    wpress = malloc(W->number_of_wells * sizeof *wpress);
    wflux  = malloc(W->well_connpos[ W->number_of_wells ] * sizeof *wflux);

    bc = allocate_flowbc(G->number_of_faces);

    set_homoperm(G->number_of_cells, G->dimensions, perm);

    cq = allocate_cq(G->number_of_cells, G->number_of_faces, nphases);
    set_incompressible(G->number_of_cells, G->number_of_faces, cq);

    h = cfs_tpfa_construct(G, W, nphases);

    vector_assign(G->number_of_cells, nphases, totmob);

    tpfa_htrans_compute(G, perm, htrans);
    tpfa_eff_trans_compute(G, totmob, htrans, trans);

    vector_zero(G->number_of_cells, src);

    dt = 1;
    cfs_tpfa_assemble(G, dt, W, bc, src, cq, trans, gravcap_f, wctrl, wdata,
                      cpress0, porevol,
                      h);

    call_UMFPACK(h->A, h->b, h->x);

    cfs_tpfa_press_flux(G, bc, W, nphases, trans, cq->phasemobf, gravcap_f,
                        wdata, h, cpress, fflux, wpress, wflux);

    cfs_tpfa_fpress(G, bc, nphases, htrans, cq->phasemobf, gravcap_f, h,
                    cpress, fflux, fpress);

    for (i = 0; i < G->number_of_cells; i++) {
        fprintf(stderr, "press(%02d) = %g;\n", i + 1, cpress[i]);
    }

    cfs_tpfa_destroy(h);
    deallocate_cq(cq);
    deallocate_flowbc(bc);

    free(wflux); free(wpress);
    free(fflux); free(fpress); free(cpress);
    free(cpress0); free(porevol); free(gravcap_f);
    free(src); free(totmob); free(perm); free(trans); free(htrans);

    deallocate_cart_grid(G);

    return 0;
}
