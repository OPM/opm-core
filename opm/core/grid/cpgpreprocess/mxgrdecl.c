//===========================================================================
//
// File: mxgrdecl.c
//
// Created: Fri Jun 19 08:48:21 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
Copyright 2009, 2010 Statoil ASA.

This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>


#include <grdecl.h>



/* Get COORD, ZCORN, ACTNUM and DIMS from mxArray.       */
/*-------------------------------------------------------*/
void mx_init_grdecl(struct grdecl *g, const mxArray *s)
{
  int i,n;
  mxArray *field;
  int numel;

  mxArray *cartdims=NULL, *actnum=NULL, *coord=NULL, *zcorn=NULL;
  
  if (!mxIsStruct(s) 
   || !(cartdims = mxGetField(s, 0, "cartDims"))
   || !(actnum   = mxGetField(s, 0, "ACTNUM"))
   || !(coord    = mxGetField(s, 0, "COORD"))
   || !(zcorn     = mxGetField(s, 0, "ZCORN"))
   )
  {
    char str[]="Input must be a single Matlab struct with fields\n"
               "cartDims, ACTNUM, COORD and ZCORN\n";
    mexErrMsgTxt(str);
  }
  
  
  numel = mxGetNumberOfElements(cartdims);
  if (!mxIsNumeric(cartdims) || numel != 3){
     mexErrMsgTxt("cartDims field must be 3 numbers");
  }
  
  double *tmp = mxGetPr(cartdims);
  n = 1;
  for (i=0; i<3; ++i){
    g->dims[i] = tmp[i];
    n      *= tmp[i];
  }


  numel = mxGetNumberOfElements(actnum);
  if (mxGetClassID(actnum) != mxINT32_CLASS ||
      numel != g->dims[0]*g->dims[1]*g->dims[2] ){
    mexErrMsgTxt("ACTNUM field must be nx*ny*nz numbers int32");
  }
  g->actnum = mxGetData(actnum);


  
  field = mxGetField(s, 0, "COORD");
  numel = mxGetNumberOfElements(coord);
  if (mxGetClassID(coord) != mxDOUBLE_CLASS || 
      numel != 6*(g->dims[0]+1)*(g->dims[1]+1)){
    mexErrMsgTxt("COORD field must have 6*(nx+1)*(ny+1) doubles.");
  }
  g->coord = mxGetPr(coord);
  

  numel = mxGetNumberOfElements(zcorn);
  if (mxGetClassID(zcorn) != mxDOUBLE_CLASS || 
      numel != 8*g->dims[0]*g->dims[1]*g->dims[2]){
    mexErrMsgTxt("ZCORN field must have 8*nx*ny*nz doubles.");
  }
  g->zcorn = mxGetPr(zcorn);
}
