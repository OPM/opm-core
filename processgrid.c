#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>


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
/*   fprintf(stderr, "%d %d %d %d\n", im,ip, jm,jp); */
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

  /* zlist may need extra space temporarily due to simple boundary treatement  */

  int    sz         = 8*(g.dims[0]+1)*(g.dims[1]+1)*g.dims[2];
  int    numpillars = (g.dims[0]+1)*(g.dims[1]+1);
  double *zlist     = malloc(sz*sizeof(*zlist));
  int    *zptr      = malloc((numpillars+1)*sizeof(*zptr));
  int    *plist     = malloc( 2*g.dims[0]*2*g.dims[1]*(2*g.dims[2]+2)*sizeof(int));


  fprintf(stderr, "Allocate %d ints for plist\n", 2*g.dims[0]*2*g.dims[1]*(2*g.dims[2]+2));

  finduniquepoints(&g, zlist, zptr, plist);


#if 0
  for (i=0; i<numpillars  ; ++i){
    fprintf(stderr, "pillar %d\n", i);
    for (k=zptr[i]; k<zptr[i+1]; ++k){
      fprintf(stderr, "%f\n", zlist[k]);
    }

    mexPrintf("\n");
  }

  for (i=0; i<8*g.n; ++i) mexPrintf("%d\n", plist[i]);
#endif


  /*======================================================*/

  fprintf(stderr, "process face geomtery\n");
  /* Process face geometry and cell-face topology */
  /* internal constant i pillar pairs */
  int n = 2 + 2*g.dims[2];
  int k;

  const int BIGNUM = 100;
  int    *faces      = malloc(BIGNUM*sz*sizeof(int));
  int    *fptr       = malloc(BIGNUM*sz*sizeof(int));



  int *pts[4];
  int *neighbors     = malloc(2* n*n* sizeof(*neighbors));
  int *intersections = malloc(4* n*n* sizeof(*intersections));
  int *work          = calloc(2* n,   sizeof(*work));
  int  d[3] = {2*g.dims[0], 2*g.dims[1], 2+2*g.dims[2]};

  int numfaces       = 0;
  int numpts         = zptr[numpillars];
  int startpts = numpts;

  fptr[0]=0;

  for (j=0; j<g.dims[1]; ++j) {
    /* Add boundary faces */
    igetvectors(d, 0, 2*j+1, plist, pts);
    findconnections(n, pts, &numpts, intersections,
		    neighbors+2*numfaces, faces, fptr, &numfaces,
		    work);

    /* internal faces */
    for (i=1; i<g.dims[0]; ++i){

      /* Vectors of point numbers */
      igetvectors(d, 2*i, 2*j+1, plist, pts);
      int startfaces = numfaces;
      findconnections(n, pts, &numpts, intersections,
		      neighbors+2*numfaces, faces, fptr, &numfaces,
		      work);



      /* Compute cell numbers from cell k-index (stored in neighbors) */
      for (k=startfaces; k<numfaces; ++k){
	if (neighbors[2*k] != -1){
	  neighbors[2*k] = i-1 + g.dims[0]*(j + g.dims[1]*neighbors[2*k]);
	}
	if (neighbors[2*k+1] != -1){

	  neighbors[2*k+1] = i + g.dims[0]*(j + g.dims[1]*neighbors[2*k+1]);
	}
      }
      /* compute intersections */

    }

    /* Add boundary face */
    igetvectors(d, 2*g.dims[0], 2*j+1, plist, pts);
    findconnections(n, pts, &numpts, intersections,
		    neighbors+2*numfaces, faces, fptr, &numfaces,
		    work);
  }


#if 1
      fprintf(stderr, "\nfaces\nnumfaces %d\n", numfaces);
      for (i=0; i<numfaces; ++i){
	for (k=fptr[i]; k<fptr[i+1]; ++k){
	  fprintf(stderr, "%d ", faces[k]);
	}
	fprintf(stderr, "\n");
      }


      fprintf(stderr, "\nneighbors\n");
      int *iptr = neighbors;
      for(k=0; k<numfaces; ++k){
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

  free (intersections);
  free (faces);
  free (fptr);
  free (neighbors);
  free (zptr);
  free (plist);
  free (zlist);
  /* Free whatever was allocated in initGrdecl. */
  freeGrdecl(&g);

}
