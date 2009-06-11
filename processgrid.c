#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>


#include "grdecl.h"
#include "uniquepoints.h"
#include "mxgrdecl.h"
#include "matalloc.h"

/* No checking of input arguments in this code! */
#define min(i,j) ((i)<(j) ? (i) : (j))
#define max(i,j) ((i)>(j) ? (i) : (j))



/*------------------------------------------------------*/
/*                                                      */
/*                                                      */
/*      Find connections for each pair of pillars       */
/*                                                      */
/*                                                      */
/*                                                      */
/*------------------------------------------------------*/


/* Determine face geometry first, then compute intersections. */
/* All intersections that occur are present in the final face geometry.*/
int *computeFaceTopology(int *a1, 
			 int *a2,
			 int *b1, 
			 int *b2,
			 int intersect[4],
			 int *faces)
{
  int mask[8];
  
  /* Which pillar points should we use? */
  if (a1[1] > b1[1]){ mask[0] = b1[1]; } else { mask[0] = a1[1]; }
  if (a2[1] > b2[1]){ mask[2] = b2[1]; } else { mask[2] = a2[1]; }
  if (a2[0] > b2[0]){ mask[4] = a2[0]; } else { mask[4] = b2[0]; }
  if (a1[0] > b1[0]){ mask[6] = a1[0]; } else { mask[6] = b1[0]; }

  /* Get shape of face: */
  /*   each new intersection will be part of the new face, */
  /*   but not all pillar points. This is encoded in mask. */


  mask[1] = intersect[0]; /* top-top */
  mask[3] = 0;
  mask[5] = intersect[3]; /* bottom-bottom*/
  mask[7] = 0;

  /* bottom-top */
  if (intersect[1]){
    if(a1[0] > b1[1]){ /* intersection[1] left of (any) intersection[0] */
      mask[0] = 0;
      mask[6] = 0;
      mask[7] = intersect[1];
    }
    else{
      mask[2] = 0;
      mask[4] = 0;
      mask[3] = intersect[1];
    }
  }

  /* top-bottom */
  if (intersect[2]){
    if(a1[1] < b1[0]){ /* intersection[2] left of (any) intersection[3] */
      mask[0] = 0;
      mask[6] = 0;
      mask[7] = intersect[2];
    }
    else{
      mask[2] = 0;
      mask[4] = 0;
      mask[3] = intersect[2];
    }
  }

  int k;
  int *f = faces;
  for (k=0; k<8; ++k){
    if(mask[k]){
      *f++ = mask[k];
    }
  }
  
  return f;
}




/* values in za and zb must be increasing. */


/* a) If we assume that the index increase when z increase for
   each pillar (but only separately), we can skip the z

   b) We assume no intersections occur on the first and last lines.



*/
/* #define overlap(a1,a2,b1,b2)      max(a1,b1) < min(a2,b2) */
#define intersection(a1,a2,b1,b2)(((a1>b1)&&(a2<b2))||((a1<b1)&&(a2>b2)))

int faceintersection(int *a1, int *a2, int *b1, int *b2)
{
/*     overlap(a1[0], a1[1], b1[0], b1[1]) || */
/*     overlap(a2[0], a2[1], b2[0], b2[1]) || */

  return
    max(a1[0],b1[0]) < min(a1[1], b1[1]) ||
    max(a2[0],b2[0]) < min(a2[1],b2[1]) ||
    intersection(a1[0], a2[0], b1[0], b2[0]);
}

/* work should be pointer to 2n ints initialised to zero . */
void findconnections(int n, int *pts[4],  int *ptnumber, int *intersections, 
		     int *neighbors, 
		     int *faces, int *fptr, int *fpos, int *work)
{
  /* vectors of point numbers for faces a(b) on pillar 1(2) */
  int *a1 = pts[0];
  int *a2 = pts[1];
  int *b1 = pts[2];
  int *b2 = pts[3];

  /* Intersection record for top line and bottomline of a */
  int *itop    = work;     /* calloc(n, sizeof(*itop)); */
  int *ibottom = work + n; /* calloc(n, sizeof(*ibottom)); */
  int *f       = faces     + fptr[*fpos];
  int *c       = neighbors;

  int k1  = 0;
  int k2  = 0;

  int i,j=0;
  int intersect[4];

  for (i = 0; i<n-1; ++i){
    if (a1[i] == a1[i+1] && a2[i] == a2[i+1]) continue;


    

    
    while(j<n-1 && (b1[j] < a1[i+1] || b2[j] < a2[i+1])){

      if (b1[j] == b1[j+1] && b2[j] == b2[j+1]){
	itop[j+1] = itop[j];
	++j;
	continue;
      }


      /* --------------------------------------------------------- */
      /* face a(i,i+1) and face b(j,j+1) have nonzero intersection */
      /* --------------------------------------------------------- */
      if (faceintersection(a1+i, a2+i, b1+j, b2+j)){


	/* Add neighbors to list of neighbors if not any first or  */
	/* last points are involved in face geometry. */
	if (!((i==0) && (j==0)) && !((i==n-2) && (j==n-2))){
	  fprintf(stderr, "here: %d %d\n", i%2 ? (i-1)/2 : -1,
		  j%2 ? (j-1)/2 : -1);
	  *c++ = i%2 ? (i-1)/2 : -1;
	  *c++ = j%2 ? (j-1)/2 : -1;
	}


	/* Completely matching faces */
	if (a1[i]==b1[j] && a1[i+1]==b1[j+1] &&
	    a2[i]==b2[j] && a2[i+1]==b2[j+1]){

	  /* Add face to list of faces if not any first or last points are involved. */
	  if (!((i==0) && (j==0)) && !((i==n-2) && (j==n-2))){
	    *f++ = a1[i];
	    *f++ = a1[i+1];
	    *f++ = a2[i+1];
	    *f++ = a2[i];
	    fptr[++(*fpos)] = f-faces;
	  }
	}

	/* Non-matching faces */
	else{
	  
	  /* Find new intersection */
	  if (intersection(a1[i+1],a2[i+1],b1[j+1],b2[j+1])) {
	    itop[j+1] = (*ptnumber)++;
	    
	    /* store point numbers of intersecting lines */
	    *intersections++ = a1[i+1];
	    *intersections++ = a2[i+1];
	    *intersections++ = b1[j+1];
	    *intersections++ = b2[j+1];


	  }else{
	    itop[j+1] = 0;
	  }

	  /* Update intersection record */
	  intersect[0] = ibottom[j  ];  /* i   x j   */
	  intersect[1] = ibottom[j+1];  /* i   x j+1 */
	  intersect[2] = itop[j  ];     /* i+1 x j   */
	  intersect[3] = itop[j+1];     /* i+1 x j+1 */


	  /* Add face to list of faces if not any first or last points are involved. */
	  if (!((i==0) && (j==0)) && !((i==n-2) && (j==n-2))){
	    f = computeFaceTopology(a1+i, a2+i, b1+j, b2+j, intersect, f);
	    fptr[++(*fpos)] = f-faces;
	  }
	}
      }

      /* Update candidates for restart of j for in next i-iteration */
      if (b1[j] < a1[i+1]) k1 = j;
      if (b2[j] < a2[i+1]) k2 = j;

      j = j+1;
    }

    /* Swap intersection records: top line of a is next bottom line of a */
    int *tmp;
    tmp = itop; itop = ibottom; ibottom = tmp;

    for(;j>min(k1,k2);--j) itop[j-1]=0;
    /* Now, j = min(k1, k2) */

  }
  
/*   free(itop); */
/*   free(ibottom); */
/*   return (cells-begin)/2; */
}




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
  fprintf(stderr, "%d %d %d %d\n", im,ip, jm,jp);
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
      
      int start = numpts;
      findconnections(n, pts, &numpts, intersections, 
		      neighbors+2*numfaces, faces, fptr, &numfaces,
		      work);
      

      
      /* Compute cell numbers from cell k-index (stored in neighbors) */
      for (k=0; k<numfaces; ++k){
	if (neighbors[2*k] != -1){
	  neighbors[2*k] = i-1 + g.dims[0]*(j + g.dims[1]*neighbors[2*k]);
	}
	if (neighbors[2*k+1] != -1){
	  
	  neighbors[2*k+1] = i + g.dims[0]*(j + g.dims[1]*neighbors[2*k+1]);
	}
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
	iptr+=2;
      }
      
      fprintf(stderr, "\nline intersections:\n");
      iptr = intersections;
      int numintersections = numpts - start;
      for(k=0; k<numintersections; ++k){
	fprintf(stderr, " (%d %d %d %d)\n",  iptr[0], iptr[1], iptr[2], iptr[3]);
	iptr+=4;
      }
#endif

      /* compute intersections */
      
    }

    /* Add boundary face */
    

  }
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

