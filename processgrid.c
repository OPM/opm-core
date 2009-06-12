#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>

#include "preprocess.h"
#include "mxgrdecl.h"

/* Gateway routine for Matlab mex function.              */
/*-------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Set up data passed from Matlab */
  struct grdecl g;
  struct processed_grid out;

  mxInitGrdecl(&g, prhs);
  processGrdecl(&g, 0.0, &out);


  if (plhs == 0){
    ;/* write to matlab struct */
  }

  /* Free whatever was allocated in initGrdecl. */
  freeGrdecl(&g);
  free_processed_grid(&out);
}
