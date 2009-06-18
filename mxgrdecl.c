#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>


#include "grdecl.h"



/* Get COORD, ZCORN, ACTNUM and DIMS from mxArray.       */
/*-------------------------------------------------------*/
void mx_init_grdecl(struct grdecl *g, const mxArray *s)
{
  int i,j,k,n;
  mxArray *field;
  int numel;

  field = mxGetField(s, 0, "cartDims");
  numel = mxGetNumberOfElements(field);
  double *tmp = mxGetPr(field);
  if (numel != 3){
    mexErrMsgTxt("cartDims field must be 3 numbers");
  }

  n = 1;
  for (i=0; i<3; ++i){
    g->dims[i] = tmp[i];
    n      *= tmp[i];
  }


  field = mxGetField(s, 0, "ACTNUM");
  numel = mxGetNumberOfElements(field);
  if (mxGetClassID(field) != mxINT32_CLASS ||
      numel != g->dims[0]*g->dims[1]*g->dims[2] ){
    mexErrMsgTxt("ACTNUM field must be nx*ny*nz numbers int32");
  }
  /* grdecl.actnum = permute(actnum, [3,1,2]);   */
  int *actnum  = mxGetData(field);
  
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
  
  
  field = mxGetField(s, 0, "COORD");
  numel = mxGetNumberOfElements(field);
  if (mxGetClassID(field) != mxDOUBLE_CLASS || 
      numel != 6*(g->dims[0]+1)*(g->dims[1]+1)){
    mexErrMsgTxt("COORD field must have 6*(nx+1)*(ny+1) doubles.");
  }
  g->coord = mxGetPr(field);
  

  field = mxGetField(s, 0, "ZCORN");
  numel = mxGetNumberOfElements(field);
  if (mxGetClassID(field) != mxDOUBLE_CLASS || 
      numel != 8*g->dims[0]*g->dims[1]*g->dims[2]){
    mexErrMsgTxt("ZCORN field must have 8*nx*ny*nz doubles.");
  }
  double *zcorn = mxGetPr(field);

  /* grdecl.zcorn = permute(zcorn, [3,1,2]);   */
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
void free_grdecl(struct grdecl *g)
{
  free((double*)g->zcorn);  g->zcorn  = NULL;
  free((double*)g->actnum); g->actnum = NULL;
}

