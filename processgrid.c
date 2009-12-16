/*===========================================================================
//
// File: processgrid.c
//
// Created: Fri Jun 19 08:46:53 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================*/

/*
Copyright 2009 SINTEF ICT, Applied Mathematics.
Copyright 2009 Statoil ASA.

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

/* Copyright 2009 SINTEF ICT, Applied Mathematics. */
/* Mex gateway by Jostein R. Natvig, SINTEF ICT.   */


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>

#include "preprocess.h"
#include "mxgrdecl.h"


void fill_grid(mxArray **out, struct processed_grid *grid)
{
  const char *names[] = {"nodes", "faces", "cells", "cellFaces", "faceNodes", "cartDims"};
  mxArray *G = mxCreateStructMatrix(1,1,sizeof names / sizeof names[0],names);

  int i,j;
  double *ptr;


  /* nodes */
  const char *n2[] = {"num", "coords"};
  mxArray *nodes  = mxCreateStructMatrix(1,1,2,n2);
  mxSetField(nodes, 0, "num", mxCreateDoubleScalar(grid->number_of_nodes));

  mxArray *nodecoords = mxCreateDoubleMatrix(grid->number_of_nodes, 3, mxREAL);
  ptr = mxGetPr(nodecoords);

  for (j=0;j<3;++j)
  {
    for(i=0; i<grid->number_of_nodes; ++i)
    {
      ptr[i+grid->number_of_nodes*j] = grid->node_coordinates[3*i+j];
    }
  }
  mxSetField(nodes, 0, "coords", nodecoords);
  mxSetField(G, 0, "nodes", nodes);


  /* faces */
  const char *n3[] = {"num", "neighbors", "numNodes", "nodePos", "tag"};
  mxArray *faces = mxCreateStructMatrix(1,1,5,n3);


  mxSetField(faces, 0, "num", mxCreateDoubleScalar(grid->number_of_faces));

  mxArray *faceneighbors = mxCreateNumericMatrix(grid->number_of_faces, 2,
						 mxINT32_CLASS, mxREAL);
  {
    int *iptr = mxGetData(faceneighbors);
    for(j=0; j<2; ++j)
    {
      for (i=0; i<grid->number_of_faces; ++i)
      {
        int ix = grid->face_neighbors[2*i+j];
        if (ix == -1)
        {
          iptr[i+grid->number_of_faces*j] = 0;
        }
	else
        {
          iptr[i+grid->number_of_faces*j] = ix+1;
        }
      }
    }
  }
  mxSetField(faces, 0, "neighbors", faceneighbors);
  mxArray *numnodes = mxCreateNumericMatrix(grid->number_of_faces,   1,
					    mxINT32_CLASS, mxREAL);
  mxArray *nodepos  = mxCreateNumericMatrix(grid->number_of_faces+1, 1,
					    mxINT32_CLASS, mxREAL);
  {
    int *iptr1 = mxGetData(numnodes);
    int *iptr2 = mxGetData(nodepos);
    iptr2[0] = 1;
    for (i=0; i<grid->number_of_faces; ++i)
    {
      iptr1[i]   = grid->face_ptr[i+1] - grid->face_ptr[i];
      iptr2[i+1] = iptr2[i] + iptr1[i];
    }
  }
  mxSetField(faces, 0, "numNodes", numnodes);
  mxSetField(faces, 0, "nodePos",  nodepos);

  mxArray *tags = mxCreateNumericMatrix(grid->number_of_faces, 1,
					mxINT32_CLASS, mxREAL);
  {
    int *iptr = mxGetData(tags);
    for (i = 0; i < grid->number_of_faces; ++i)
    {
        iptr[i] = grid->face_tag[i] + 1;
    }
  }
  mxSetField(faces, 0, "tag", tags);

  mxSetField(G, 0, "faces", faces);

  const char *n4[] = {"num", "facePos", "indexMap"};
  mxArray *cells = mxCreateStructMatrix(1,1,3,n4);

  mxSetField(cells, 0, "num", mxCreateDoubleScalar(grid->number_of_cells));

  mxArray *map = mxCreateNumericMatrix(grid->number_of_cells, 1,
				       mxINT32_CLASS, mxREAL);
  {
    int *iptr = mxGetData(map);
    for(i=0; i<grid->number_of_cells; ++i)
    {
      iptr[i] = grid->local_cell_index[i]+1;
    }
  }
  mxSetField(cells, 0, "indexMap", map);

  mxArray *facepos  = mxCreateNumericMatrix(grid->number_of_cells+1, 1,
					    mxINT32_CLASS, mxREAL);
  {
    int *iptr = mxGetData(facepos);
    for(i=0; i<grid->number_of_cells+1; ++i)
    {
      iptr[i] = 0;
    }
    for (i=0; i<2*grid->number_of_faces; ++i)
    {
      int c=grid->face_neighbors[i];
      if(c != -1)
      {
        iptr[c+1]++;
      }
    }
    iptr[0] = 1;
    for(i=0; i<grid->number_of_cells; ++i)
    {
      iptr[i+1] += iptr[i];
    }
  }
  mxSetField(cells, 0, "facePos",  facepos);

  mxSetField(G, 0, "cells", cells);


  int *counter       = calloc(grid->number_of_cells, sizeof(*counter));
  int num_half_faces = 0;
  {
    int *iptr = mxGetData(facepos);
    for(i=0; i<grid->number_of_cells; ++i)
    {
      counter[i] = num_half_faces;
      num_half_faces += iptr[i+1]-iptr[i];
    }
  }

  mxArray *cellfaces = mxCreateNumericMatrix(num_half_faces, 1,
					     mxINT32_CLASS, mxREAL);
  {
    int *iptr = mxGetData(cellfaces);
    for (i=0; i<grid->number_of_faces; ++i)
    {
      int c1 = grid->face_neighbors[2*i];
      int c2 = grid->face_neighbors[2*i+1];
      if(c1 != -1) iptr[counter[c1]++] = i+1;
      if(c2 != -1) iptr[counter[c2]++] = i+1;
    }
  }
  mxSetField(G, 0, "cellFaces", cellfaces);


  int n = grid->face_ptr[grid->number_of_faces];

  mxArray *facenodes = mxCreateNumericMatrix(n, 1,
					     mxINT32_CLASS, mxREAL);
  {
    int *iptr = mxGetData(facenodes);
    for (i=0; i<n; ++i)
    {
      iptr[i] = grid->face_nodes[i]+1;
    }
  }
  mxSetField(G, 0, "faceNodes", facenodes);
  free(counter);

  mxArray *cartDims = mxCreateDoubleMatrix(1, 3, mxREAL);
  ptr    = mxGetPr(cartDims);
  ptr[0] = grid->dimensions[0];
  ptr[1] = grid->dimensions[1];
  ptr[2] = grid->dimensions[2];

  mxSetField(G, 0, "cartDims", cartDims);

  out[0] = G;
}


/* Gateway routine for Matlab mex function.              */
/*-------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Set up data passed from Matlab */
  struct grdecl g;
  struct processed_grid out;
  double tolerance = 0.0;

  mx_init_grdecl(&g, prhs[0]);
  if (nrhs == 2)
  {
    tolerance = mxGetScalar (prhs[1]);
  }

  process_grdecl(&g, tolerance, &out);


  if (plhs >0)
  {
    /* write to matlab struct */
    fill_grid(plhs, &out);
  }


  /* Free whatever was allocated in initGrdecl. */
  free_processed_grid(&out);
}
