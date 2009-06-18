#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <limits.h>

#include "sparsetable.h"
#include "facetopology.h"

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
static int *computeFaceTopology(int *a1,
				int *a2,
				int *b1,
				int *b2,
				int intersect[4],
				int *faces)
{
  int mask[8] = {-1};

  /* Which pillar points should we use? */
  if (a1[1] > b1[1]){ mask[0] = b1[1]; } else { mask[0] = a1[1]; }
  if (a2[1] > b2[1]){ mask[2] = b2[1]; } else { mask[2] = a2[1]; }
  if (a2[0] > b2[0]){ mask[4] = a2[0]; } else { mask[4] = b2[0]; }
  if (a1[0] > b1[0]){ mask[6] = a1[0]; } else { mask[6] = b1[0]; }

  /* Get shape of face: */
  /*   each new intersection will be part of the new face, */
  /*   but not all pillar points. This is encoded in mask. */


  mask[1] = intersect[3]; /* top-top */
  mask[3] = -1;
  mask[5] = intersect[0]; /* bottom-bottom*/
  mask[7] = -1;

  /* bottom-top */
  if (intersect[1] != -1){
    if(a1[0] > b1[1]){ /* intersection[1] left of (any) intersection[0] */
      mask[0] = -1;
      mask[6] = -1;
      mask[7] = intersect[1];
    }
    else{
      mask[2] = -1;
      mask[4] = -1;
      mask[3] = intersect[1];
    }
  }

  /* top-bottom */
  if (intersect[2] != -1){
    if(a1[1] < b1[0]){ /* intersection[2] left of (any) intersection[3] */
      mask[0] = -1;
      mask[6] = -1;
      mask[7] = intersect[2];
    }
    else{
      mask[2] = -1;
      mask[4] = -1;
      mask[3] = intersect[2];
    }
  }
  int k;
  int *f = faces;
  for (k=0; k<8; ++k){
    if(mask[k] != -1){
      *f++ = mask[k];
    }
  }
  return f;
}




/* a) If we assume that the index increase when z increase for
      each pillar (but only separately), we can use only the point indices.

   b) We assume no intersections occur on the first and last lines.
      This is convenient in the identification of (unique) intersections.

*/

#define lineintersection(a1,a2,b1,b2)(((a1>b1)&&(a2<b2))||((a1<b1)&&(a2>b2)))
static int faceintersection(int *a1, int *a2, int *b1, int *b2)
{
  return
    max(a1[0],b1[0]) < min(a1[1],b1[1]) ||
    max(a2[0],b2[0]) < min(a2[1],b2[1]) ||
    lineintersection(a1[0], a2[0], b1[0], b2[0]);
}


#define meaningful_face !((a1[i]  ==INT_MIN) && (b1[j]  ==INT_MIN)) && \
                        !((a1[i+1]==INT_MAX) && (b1[j+1]==INT_MAX))

/* work should be pointer to 2n ints initialised to zero . */
void findconnections(int n, int *pts[4],
		     int *ptnumber,
		     int *intersectionlist,
		     int *neighbors,
		     int *work,
		     sparse_table_t *ftab)
{
  /* vectors of point numbers for faces a(b) on pillar 1(2) */
  int *a1 = pts[0];
  int *a2 = pts[1];
  int *b1 = pts[2];
  int *b2 = pts[3];

  /* Intersection record for top line and bottomline of a */
  int *itop    = work;
  int *ibottom = work + n;
  int *f       = (int*)ftab->data     + ftab->ptr[ftab->position];
  int *c       = neighbors + 2*ftab->position;

  int k1  = 0;
  int k2  = 0;

  int i,j=0;
  int intersect[4]= {-1};


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


	/* Completely matching faces */
	if (a1[i]==b1[j] && a1[i+1]==b1[j+1] &&
	    a2[i]==b2[j] && a2[i+1]==b2[j+1]){

	  /* Add face to list of faces if not any first or last points are involved. */
	  if (meaningful_face){

	    int cell_a = i%2 ? (i-1)/2 : -1;
	    int cell_b = j%2 ? (j-1)/2 : -1;

	    if (cell_a != -1 || cell_b != -1){
	      *c++ = cell_a;
	      *c++ = cell_b;

	      /* face */
	      *f++ = a1[i];
	      *f++ = a1[i+1];
	      *f++ = a2[i+1];
	      *f++ = a2[i];
	      ftab->ptr[++ftab->position] = f - (int*)ftab->data;
	    }
	    else{
	      ;
	      /*
		fprintf(stderr, 
		      "Warning. For some reason we get face connecting void spaces\n");
	      */
	    }

	  }
	}

	/* Non-matching faces */
	else{

	  /* Find new intersection */
	  if (lineintersection(a1[i+1],a2[i+1],b1[j+1],b2[j+1])) {
	    itop[j+1] = (*ptnumber)++;

	    /* store point numbers of intersecting lines */
	    *intersectionlist++ = a1[i+1];
	    *intersectionlist++ = a2[i+1];
	    *intersectionlist++ = b1[j+1];
	    *intersectionlist++ = b2[j+1];


	  }else{
	    itop[j+1] = -1;
	  }

	  /* Update intersection record */
	  intersect[0] = ibottom[j  ];  /* i   x j   */
	  intersect[1] = ibottom[j+1];  /* i   x j+1 */
	  intersect[2] = itop[j  ];     /* i+1 x j   */
	  intersect[3] = itop[j+1];     /* i+1 x j+1 */


	  /* Add face to list of faces if no INT_MIN or INT_MAX appear in a or b. */
	  if (meaningful_face){

	    /* 
	       Even indices refer to space between cells, 
	       odd indices refer to cells 
	    */
	    int cell_a = i%2 ? (i-1)/2 : -1;
	    int cell_b = j%2 ? (j-1)/2 : -1;



	    if (cell_a != -1 || cell_b != -1){
	      *c++ = cell_a;
	      *c++ = cell_b;
	      
	      f = computeFaceTopology(a1+i, a2+i, b1+j, b2+j, intersect, f);
	      ftab->ptr[++ftab->position] = f - (int*)ftab->data;
	    }
	    else{
	      ;
	      /*
		fprintf(stderr, 
		      "Warning. For some reason we get face connecting void spaces\n");
	      */
	    }
	  }
	}
      }

      /* Update candidates for restart of j for in next i-iteration */
      if (b1[j] < a1[i+1]) k1 = j;
      if (b2[j] < a2[i+1]) k2 = j;

      j = j+1;
    }


    
    /* Swap intersection records: top line of a[i,i+1] is bottom line of a[i+1,i+2] */
    int *tmp; tmp = itop; itop = ibottom; ibottom = tmp;

    /* Zero out the "new" itop */
    for(j=0;j<n; ++j) itop[j]=-1;

    /* Set j to appropriate start position for next i */
    j = min(k1, k2);
  }
}
