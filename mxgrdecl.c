#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>


#include "grdecl.h"



/* Get COORD, ZCORN, ACTNUM and DIMS from mxArray.       */
/*-------------------------------------------------------*/
void mxInitGrdecl(struct grdecl *g, const mxArray *prhs[])
{
  int i,j,k,n;


  g->coord    = mxGetPr(mxGetField(prhs[0], 0, "COORD"));
  double *tmp = mxGetPr(mxGetField(prhs[0], 0, "cartDims"));
  n = 1;
  for (i=0; i<3; ++i){
    g->dims[i] = tmp[i];
    n      *= tmp[i];
  }
/*   mexPrintf("dimensions: %d %d %d\n", */
/* 	    g->dims[0], */
/* 	    g->dims[1], */
/* 	    g->dims[2]); */



  /* grdecl.actnum = permute(actnum, [3,1,2]);   */
  int *actnum  = mxGetData(mxGetField(prhs[0], 0, "ACTNUM"));
  
  int *a = malloc(n*  sizeof(*g->actnum));
  int *iptr = a;
  for (j=0; j<g->dims[1]; ++j){
    for (i=0; i<g->dims[0]; ++i){
      for (k=0; k<g->dims[2]; ++k){
	*iptr++ = actnum[i+g->dims[0]*(j+g->dims[1]*k)];
      }
    }
  }
  g->actnum = a;

  /* grdecl.zcorn = permute(zcorn, [3,1,2]);   */
  double *zcorn = mxGetPr(mxGetField(prhs[0], 0, "ZCORN"));
  double *z = malloc(n*8*sizeof(*g->zcorn));
  double *dptr = z;
  for (j=0; j<2*g->dims[1]; ++j){
    for (i=0; i<2*g->dims[0]; ++i){
      for (k=0; k<2*g->dims[2]; ++k){
	*dptr++ = zcorn[i+2*g->dims[0]*(j+2*g->dims[1]*k)];
      }
    }
  }
  g->zcorn = z;
}



/* Free stuff that was allocated in initgrdecl.          */
/*-------------------------------------------------------*/
void freeGrdecl(struct grdecl *g)
{
  free((double*)g->zcorn);  g->zcorn  = NULL;
  free((double*)g->actnum); g->actnum = NULL;
}

