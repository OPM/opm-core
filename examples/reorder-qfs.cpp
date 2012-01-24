/* Copyright 2011 (c) Jostein R. Natvig <jostein.natvig@gmail.com> */




#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/grid.h>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/linalg/call_umfpack.h>
#include <opm/core/utility/cart_grid.h>

#include <opm/core/transport/reorder/twophasetransport.hpp>


struct Rock {
    double *perm;
    double *poro;
};

static void
destroy_rock(struct Rock *rock)
{
    if (rock!=NULL)
    {
        free(rock->perm);
        free(rock->poro);
    }
    free(rock);
}

static struct Rock*
init_rock(int nc, int d)
{
    int i,j;
    struct Rock *rock = (struct Rock*) malloc(sizeof *rock);

    if (rock!=NULL)
    {
        rock->perm = (double*) malloc(nc*d*d*sizeof *rock->perm);
        rock->poro = (double*) malloc(nc * sizeof *rock->poro);

        if ((rock->perm==NULL) || (rock->poro==NULL))
        {
            destroy_rock(rock);
            rock = NULL;
        }
        else
        {
            vector_zero(nc*d*d, rock->perm);
            for (i=0; i<nc; ++i)
            {
                rock->poro[i] = 1.0;
                for (j=0; j<d; ++j)
                {
                    rock->perm[j*(d+1)+i*d*d] = 1.0;
                }
            }
        }
    }

    return rock; 
}


static void
qfs(grid_t *g, struct Rock *rock, double *src, 
    double *pressure, double *faceflux)
{
    size_t  alloc_sz, totconn;
    double *htrans, *trans, *gpress;
    struct ifs_tpfa_data *h;

    fprintf(stderr, "QFS: Allocate space\n");
    totconn = g->cell_facepos[ g->number_of_cells ];

    alloc_sz  = totconn;            /* htrans */
    alloc_sz += g->number_of_faces; /* trans */
    alloc_sz += totconn;            /* gpress */

    htrans = (double*) malloc(alloc_sz * sizeof *htrans);
    trans  = htrans + totconn;
    gpress = trans  + g->number_of_faces;

    h = ifs_tpfa_construct(g);

    fprintf(stderr, "QFS: Compute transmissibilities\n");
    tpfa_htrans_compute(g, rock->perm, htrans);
    tpfa_trans_compute (g, htrans    , trans );

    vector_zero(totconn, gpress);

    assert (g->number_of_cells > 1);

    fprintf(stderr, "QFS: Assemble pressure system\n");
    ifs_tpfa_assemble(g, trans, src, gpress, h);

    fprintf(stderr, "QFS: Solve pressure system\n");
    call_UMFPACK(h->A, h->b, h->x);

    fprintf(stderr, "QFS: Recover pressure and flux\n");
    ifs_tpfa_press_flux(g, trans, h, pressure, faceflux);

    fprintf(stderr, "QFS: Release memory resources\n");
    ifs_tpfa_destroy(h);
    free(htrans);

    fprintf(stderr, "QFS: Done\n");
}

static void
compute_porevolume(struct UnstructuredGrid *g, double *poro, 
                   double *pv)
{
    int i;
    for (i=0; i<g->number_of_cells; ++i)
    {
        pv[i] = g->cell_volumes[i]*poro[i];
    }
}

int main(void)
{
    struct UnstructuredGrid *g = create_cart_grid_2d(10,10);
    struct Rock *rock = init_rock(g->number_of_cells, g->dimensions);
    double *sat   = (double*) malloc(g->number_of_cells *sizeof *sat);
    double *press = (double*) malloc(g->number_of_cells *sizeof *press);
    double *flux  = (double*) malloc(g->number_of_faces *sizeof *flux);
    double *pv    = (double*) malloc(g->number_of_cells *sizeof *pv);
    double *src   = (double*) malloc(g->number_of_cells *sizeof *src);

    vector_zero(g->number_of_cells, sat);
    vector_zero(g->number_of_cells, src);
    src[0]=1;src[g->number_of_cells-1]=-1;

    qfs(g, rock, src, press, flux);
    
    compute_porevolume(g, rock->poro, pv);

    twophasetransport(pv, src, 10, g, flux, NULL, sat);

    vector_write(g->number_of_cells, sat, "saturation.txt");

    free(pv);
    free(src);
    free(sat);
    free(press);
    free(flux);
    destroy_rock(rock);
    return 0;
}
