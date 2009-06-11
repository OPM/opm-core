#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>

#include "sparsetable.h"
#include "grdecl.h"
#include "uniquepoints.h"
#include "mxgrdecl.h"
#include "matalloc.h"
#include "facetopology.h"

/* No checking of input arguments in this code! */
#define min(i,j) ((i)<(j) ? (i) : (j))
#define max(i,j) ((i)>(j) ? (i) : (j))



void interpolate_pillar(double *coord, double *x, double *y, double *z)
{
  double a = (*z-coord[2])/(coord[5]-coord[2]);
  *x       = coord[0] + a*(coord[3]-coord[0]);
  *y       = coord[1] + a*(coord[4]-coord[1]);
}


/*-------------------------------------------------------*/
static void igetvectors(int dims[3], int i, int j, int *field, int *v[])
{

  int im = max(1,       i  ) - 1;
  int ip = min(dims[0], i+1) - 1;
  int jm = max(1,       j  ) - 1;
  int jp = min(dims[1], j+1) - 1;

  v[0] = field + dims[2]*(im + dims[0]* jm);
  v[1] = field + dims[2]*(im + dims[0]* jp);
  v[2] = field + dims[2]*(ip + dims[0]* jm);
  v[3] = field + dims[2]*(ip + dims[0]* jp);
}




/* Gateway routine for Matlab mex function.              */
/*-------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i,j;

  /* Set up data passed from Matlab */
  struct Grdecl g;
  mxInitGrdecl(&g, prhs);


  /* Code below assumes k index runs fastests, ie. that dimensions of
     table is permuted to (dims[2], dims[0], dims[1]) */
  /* ---------------------------------------------------------------*/


  /* Set up space for return values */


  /* Unstructured storage of unique zcorn values for each pillar. */
  /* ztab->data may need extra space temporarily due to simple boundary treatement  */
  int            sz         = 8*(g.dims[0]+1)*(g.dims[1]+1)*g.dims[2]; /* Big for convenient bc */
  int            numpillars = (g.dims[0]+1)*(g.dims[1]+1);
  sparse_table_t *ztab      = malloc_sparse_table(numpillars, sz, sizeof(double));
  int    *plist             = malloc( 2*g.dims[0]*2*g.dims[1]*(2*g.dims[2]+2)*sizeof(int));

  finduniquepoints(&g, plist, ztab);


  /*======================================================*/

  fprintf(stderr, "process face geomtery\n");
  /* Process face geometry and cell-face topology */
  /* internal constant i pillar pairs */
  int n = 2 + 2*g.dims[2];
  int k;

  const int BIGNUM = 100;
  /* Unstructured storage of face nodes */
  sparse_table_t *ftab = malloc_sparse_table(BIGNUM*sz*sizeof(int),
					     BIGNUM*sz*sizeof(int),
					     sizeof(int));



  int *pts[4];
  int *neighbors     = malloc(2* n*n* sizeof(*neighbors));
  int *intersections = malloc(4* n*n* sizeof(*intersections));

  int *work          = calloc(2* n,   sizeof(*work));
  int  d[3] = {2*g.dims[0], 2*g.dims[1], 2+2*g.dims[2]};

  int startface;
  int numpts         = ztab->ptr[numpillars];
  int startpts       = numpts;
  int isodd;
  ftab->ptr[0]=0;

  for (j=0; j<g.dims[1]; ++j) {
    /* ------------------ */
    /* Add boundary faces */
    /* ------------------ */
    igetvectors(d, 0, 2*j+1, plist, pts);
    startface = ftab->position;
    findconnections(n, pts, &numpts, intersections, neighbors, work, ftab);


    /* odd */
    for(k=2*startface+1; k<2*ftab->position; k+=2){
      if (neighbors[k] != -1){
	neighbors[k] = g.dims[0]*(j + g.dims[1]*neighbors[k]);
      }
    }

    /* even */
    for(k=2*startface; k<2*ftab->position; k+=2){
      neighbors[k] = -1;
    }


    /* -------------- */
    /* internal faces */
    /* -------------- */
    for (i=1; i<g.dims[0]; ++i){

      /* Vectors of point numbers */
      igetvectors(d, 2*i, 2*j+1, plist, pts);
      startface = ftab->position;
      findconnections(n, pts, &numpts, intersections, neighbors, work, ftab);



      /* Compute cell numbers from cell k-index (stored in neighbors) */
      /* For even k, neighbors[k] refers to stack (i-1,j),
	 for odd  k, neighbors[k] refers to stack (i,j) of cells */
      isodd = 1;
      for(k=2*startface; k<2*ftab->position; ++k){
	isodd = !isodd;
	if (neighbors[k] != -1){
	  neighbors[k] = i-1+isodd + g.dims[0]*(j + g.dims[1]*neighbors[k]);
	}
      }

      /* compute intersections */

    }

    /* ------------------ */
    /* Add boundary face */
    /* ------------------ */
    startface = ftab->position;
    igetvectors(d, 2*g.dims[0], 2*j+1, plist, pts);
    findconnections(n, pts, &numpts, intersections, neighbors, work, ftab);

    /* even indices */
    for(k=2*startface; k<2*ftab->position; k+=2){
      if (neighbors[k] != -1){
	neighbors[k] = g.dims[0]-1 + g.dims[0]*(j + g.dims[1]*neighbors[k]);
      }
    }
    /* odd indices */
    for(k=2*startface+1; k<2*ftab->position; k+=2){
      neighbors[k] = -1;
    }
  }



  /*                                                       */
  /*                                                       */
  /*                                                       */
  /*                 D E B U G    CODE                     */
  /*                                                       */
  /*                                                       */
  /*                                                       */

#if 1
      fprintf(stderr, "\nfaces\nnumfaces %d\n", ftab->position);
      for (i=0; i<ftab->position; ++i){
	for (k=ftab->ptr[i]; k<ftab->ptr[i+1]; ++k){
	  fprintf(stderr, "%d ", ((int*)ftab->data)[k]);
	}
	fprintf(stderr, "\n");
      }


      fprintf(stderr, "\nneighbors\n");
      int *iptr = neighbors;
      for(k=0; k<ftab->position; ++k){
	fprintf(stderr, " (%d %d)\n",  iptr[0], iptr[1]);
	++iptr;
	++iptr;
      }

      fprintf(stderr, "\nline intersections:\n");
      iptr = intersections;
      int numintersections = numpts - startpts;
      for(k=0; k<numintersections; ++k){
	fprintf(stderr, " (%d %d %d %d)\n",  iptr[0], iptr[1], iptr[2], iptr[3]);
	iptr+=4;
      }
#endif



#if 0
  int xstride = 2*g.dims[2];
  int ystride = 4*g.dims[2]*g.dims[0];
  int zstride = 1;

  /* internal constant j pillar pairs */
  for (i=0; i<g.dims[0]; ++i){
    for (j=0; j<g.dims[1]-1; ++j){
      int p1 = 2*i + ystride* (2*j+1);
      int p2 = p1  + xstride;
      process_pair(p1, p2, g.dims, g.zcorn, ystride, zstride);
    }
  }
#endif


  free_sparse_table(ztab);
  free_sparse_table(ftab);

  free (intersections);
  free (neighbors);
  free (plist);
  /* Free whatever was allocated in initGrdecl. */
  freeGrdecl(&g);

}
