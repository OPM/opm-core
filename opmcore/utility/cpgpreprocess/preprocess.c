/*===========================================================================
//
// File: preprocess.c
//
// Created: Fri Jun 19 08:42:39 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//==========================================================================*/

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
#include <stdio.h>
#include <float.h>
#include "preprocess.h"
#include "uniquepoints.h"
#include "facetopology.h"


#define min(i,j) ((i)<(j) ? (i) : (j))
#define max(i,j) ((i)>(j) ? (i) : (j))



/*-----------------------------------------------------------------
  Given a vector <field> with k index running faster than i running
  faster than j, and Cartesian dimensions <dims>, find pointers to the
  (i-1, j-1, 0), (i-1, j, 0), (i, j-1, 0) and (i, j, 0) elements of
  field.
 */
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






/*-----------------------------------------------------------------
  Special purpose

  Convert from k-index to Cartesian index i+nx*(j+ny*k) for every other
  element in neighbors.

 */
static void compute_cell_index(const int dims[3], int i, int j, int *neighbors, int len)
{

  int k;
  if (i<0 || i>=dims[0] || j<0 || j >= dims[1]){
    for(k=0; k<len; k+=2){
      neighbors[k] = -1;
    }
  }else{
    for(k=0; k<len; k+=2){
      if (neighbors[k] != -1){
	int tmp = i + dims[0]*(j + dims[1]*neighbors[k]);
	neighbors[k] = tmp;
      }
    }
  }
}


/*-----------------------------------------------------------------
  Ensure there's sufficient memory
 */
static int checkmemeory(int nz, struct processed_grid *out, int **intersections)
{

  /* Ensure there is enough space */
  int r = (2*nz+2)*(2*nz+2);
  int m = out->m;
  int n = out->n;

  if(out->number_of_faces +  r > m){
    m += max(m*0.5, 2*r);
  }
  if (out->face_ptr[out->number_of_faces] + 6*r > n){
    n += max(n*0.5, 12*r);
  }

  if (m != out->m){
    void *p1 = realloc(out->face_neighbors, 2*m * sizeof(*out->face_neighbors));
    void *p2 = realloc(*intersections, 4*m * sizeof(**intersections));
    if (p1 && p2) {
      out->face_neighbors = p1;
      *intersections = p2;
    } else {
      return 0;
    }
  }

  if (m != out->m || n != out->n) {
    void *p1 = realloc(out->face_nodes, n * sizeof *out->face_nodes);
    void *p2 = realloc(out->face_ptr, (m+1) * sizeof *out->face_ptr);
    void *p3 = realloc(out->face_tag, m * sizeof *out->face_tag);

    if (p1 && p2 && p3) {
      out->face_nodes = p1;
      out->face_ptr = p2;
      out->face_tag = p3;
    } else {
      return 0;
    }

    out->m = m;
    out->n = n;
  }

  return 1;
}

/*-----------------------------------------------------------------
  For each vertical face (i.e. i or j constant),
       -find point numbers for the corners and
       -cell neighbors.
       -new points on faults defined by two intgersecting lines.

  direction == 0 : constant-i faces.
  direction == 1 : constant-j faces.
 */
static void process_vertical_faces(int direction,
			   int **intersections,
			    int *plist, int *work,
			    struct processed_grid *out)
{
  int i,j;
  int *cornerpts[4];

  enum face_tag tag[] = { LEFT, BACK };

  int nx = out->dimensions[0];
  int ny = out->dimensions[1];
  int nz = out->dimensions[2];

  for (j=0; j<ny+direction; ++j) {
    for (i=0; i<nx+1-direction; ++i){

      if (!checkmemeory(nz, out, intersections)){
	fprintf(stderr, "Could not allocat enough space in process_vertical_faces\n");
	exit(1);
      }

      /* Vectors of point numbers */
      int d[] = {2*nx, 2*ny, 2+2*nz};
      igetvectors(d, 2*i+direction, 2*j+1-direction, plist, cornerpts);

      if(direction==1){
	/* 1   3       0   1    */
	/*       --->           */
	/* 0   2       2   3    */
	/* rotate clockwise     */
	int *tmp     = cornerpts[1];
	cornerpts[1] = cornerpts[0];
	cornerpts[0] = cornerpts[2];
	cornerpts[2] = cornerpts[3];
	cornerpts[3] = tmp;
      }

      /* int startface = ftab->position; */
      int startface = out->number_of_faces;
      /* int num_intersections = *npoints - npillarpoints; */
      int num_intersections = out->number_of_nodes -
	                      out->number_of_nodes_on_pillars;

      findconnections(2*nz+2, cornerpts,
		      *intersections+4*num_intersections,
		      work, out);

      int *ptr = out->face_neighbors + 2*startface;
      int len = 2*out->number_of_faces - 2*startface;

      compute_cell_index(out->dimensions, i-1+direction, j-direction, ptr,   len);
      compute_cell_index(out->dimensions, i,   j, ptr+1, len);

      /* Tag the new faces */
      int f = startface;
      for (; f < out->number_of_faces; ++f) {
        out->face_tag[f] = tag[direction];
      }
    }
  }
}

static int linearindex(const int dims[3], int i, int j, int k)
{
  return i+dims[0]*(j+dims[1]*k);
}


/*-----------------------------------------------------------------
  For each horizontal face (i.e. k constant),
       -find point numbers for the corners and
       -cell neighbors.

  Also define map from logically Cartesian
  cell index to local cell index 0, ..., #<active cells>.   Exclude
  cells that are have collapsed coordinates. (This includes cells with
  ACTNUM==0)

 */
static void process_horizontal_faces(int **intersections,
			      int *plist,
			      struct processed_grid *out)
{
  int i,j,k;

  int nx = out->dimensions[0];
  int ny = out->dimensions[1];
  int nz = out->dimensions[2];

  int *cell  = out->local_cell_index;
  int cellno = 0;

  /* dimensions of plist */
  int  d[] = {2*nx, 2*ny, 2+2*nz};


  for(j=0; j<ny; ++j) {
    for (i=0; i<nx; ++i) {


      if (!checkmemeory(nz, out, intersections)){
	fprintf(stderr, "Could not allocat enough space in process_horizontal_faces\n");
	exit(1);
      }


      int *f = out->face_nodes     + out->face_ptr[out->number_of_faces];
      int *n = out->face_neighbors + 2*out->number_of_faces;


      /* Vectors of point numbers */
      int *c[4];
      igetvectors(d, 2*i+1, 2*j+1, plist, c);

      int prevcell = -1;


      for (k = 1; k<nz*2+1; ++k){

	/* Skip if space between face k and face k+1 is collapsed. */
	/* Note that inactive cells (with ACTNUM==0) have all been  */
	/* collapsed in finduniquepoints.                           */
	if (c[0][k] == c[0][k+1] && c[1][k] == c[1][k+1] &&
	    c[2][k] == c[2][k+1] && c[3][k] == c[3][k+1]){

	  /* If the pinch is a cell: */
	  if (k%2){
	    int idx = linearindex(out->dimensions, i,j,(k-1)/2);
	    cell[idx] = -1;
	  }
	}
	else{

	  if (k%2){
	    /* Add face */
	    *f++ = c[0][k];
	    *f++ = c[2][k];
	    *f++ = c[3][k];
	    *f++ = c[1][k];

	    out->face_tag[  out->number_of_faces] = TOP;
	    out->face_ptr[++out->number_of_faces] = f - out->face_nodes;

	    int thiscell = linearindex(out->dimensions, i,j,(k-1)/2);
	    *n++ = prevcell;
	    *n++ = prevcell = thiscell;

	    cell[thiscell] = cellno++;

	  }
	  else{
	    if (prevcell != -1){
	      /* Add face */
	      *f++ = c[0][k];
	      *f++ = c[2][k];
	      *f++ = c[3][k];
	      *f++ = c[1][k];

              out->face_tag[  out->number_of_faces] = TOP;
	      out->face_ptr[++out->number_of_faces] = f - out->face_nodes;

	      *n++ = prevcell;
	      *n++ = prevcell = -1;
	    }
	  }
	}
      }
    }
  }
  out->number_of_cells = cellno;
}


/*-----------------------------------------------------------------
  On input,
  L points to 4 ints that indirectly refers to points in c.
  c points to array of coordinates [x0,y0,z0,x1,y1,z1,...,xn,yn,zn].
  pt points to array of 3 doubles.

  On output,
  pt holds coordinates to intersection between lines given by point
  numbers L[0]-L[1] and L[2]-L[3].
*/
static void approximate_intersection_pt(int *L, double *c, double *pt)
{

  double z0 = c[3*L[0]+2];
  double z1 = c[3*L[1]+2];
  double z2 = c[3*L[2]+2];
  double z3 = c[3*L[3]+2];

  double a = (z2-z0)/(z1-z0 - (z3-z2));
  if (isinf(a) || isnan(a)){
    a = 0;
  }


  pt[0] = c[3*L[0]+0]* (1.0-a) + c[3*L[1]+0]* a;
  pt[1] = c[3*L[0]+1]* (1.0-a) + c[3*L[1]+1]* a;
  pt[2] = c[3*L[0]+2]* (1.0-a) + c[3*L[1]+2]* a;

}

/*-----------------------------------------------------------------
  Compute x,y and z coordinates for points on each pillar.
  Then, append x,y and z coordinates for extra points on faults.
 */
static void
compute_intersection_coordinates(int                   *intersections,
				 struct processed_grid *out)
{
  int n  = out->number_of_nodes;
  int np = out->number_of_nodes_on_pillars;


  /* Make sure the space allocated for nodes match the number of node. */
  void *p = realloc (out->node_coordinates, 3*n*sizeof(double));
  if (p) {
    out->node_coordinates = p;
  }
  else {
    fprintf(stderr, "Could not allocate extra space for intersections\n");
  }


  /* Append intersections */
  int    k;
  double *pt    = out->node_coordinates + 3*np;
  int    *itsct = intersections;

  for (k=np; k<n; ++k){
    approximate_intersection_pt(itsct, out->node_coordinates, pt);
    pt    += 3;
    itsct += 4;

  }
}


/*-----------------------------------------------------------------
  Public interface
 */
void process_grdecl(const struct grdecl   *in,
		    double                tolerance,
		    struct processed_grid *out)
{
  int i,j,k;

  int nx = in->dims[0];
  int ny = in->dims[1];
  int nz = in->dims[2];

  int nc = nx*ny*nz;
  struct grdecl g;

  g.dims[0] = nx;
  g.dims[1] = ny;
  g.dims[2] = nz;
  int    *actnum  = malloc (nc *     sizeof *actnum);
  double *zcorn   = malloc (nc * 8 * sizeof *zcorn);

  /* Permute actnum */
  int *iptr = actnum;
  for (j=0; j<ny; ++j){
    for (i=0; i<nx; ++i){
      for (k=0; k<nz; ++k){
	*iptr++ = in->actnum[i+nx*(j+ny*k)];
      }
    }
  }
  g.actnum = actnum;


  /* HACK */
  /* Check that ZCORN is strictly nodecreaing along pillars.  If */
  /* not, check if -ZCORN is strictly nondecreasing.             */
  int sign, error;
  for (sign = 1; sign>-2; sign = sign - 2){
    error = 0;

    /* Ensure that zcorn is strictly nondecreasing in k-direction */
    for (j=0; j<2*ny; ++j){
      for (i=0; i<2*nx; ++i){
	for (k=0; k<2*nz-1; ++k){
	  double z1 = sign*in->zcorn[i+2*nx*(j+2*ny*(k))];
	  double z2 = sign*in->zcorn[i+2*nx*(j+2*ny*(k+1))];
	  
	  if (z2 < z1){
	    fprintf(stderr, "\nZCORN should be strictly nondecreasing along pillars!\n");	 
	    error = 1;
	    goto end;
	  }
	}
      }
    }

  end:
    if (!error){
      break;
    }
  }

  if (error){
    fprintf(stderr, "Attempt to reverse sign in ZCORN failed.\n"
		    "Grid definition may be broken\n");
  }
  
  /* Permute zcorn */
  double *dptr = zcorn;
  for (j=0; j<2*ny; ++j){
    for (i=0; i<2*nx; ++i){
      for (k=0; k<2*nz; ++k){
	*dptr++ = sign*in->zcorn[i+2*nx*(j+2*ny*k)];
      }
    }
  }




  g.zcorn = zcorn;
  g.coord = in->coord;


  /* Code below assumes k index runs fastests, ie. that dimensions of
     table is permuted to (dims[2], dims[0], dims[1]) */

  const int BIGNUM = 64;
  out->m                = BIGNUM/3;
  out->n                = BIGNUM;

  out->face_neighbors   = malloc( BIGNUM      * sizeof *out->face_neighbors);
  out->face_nodes       = malloc( out->n      * sizeof *out->face_nodes);
  out->face_ptr         = malloc((out->m + 1) * sizeof *out->face_ptr);
  out->face_tag         = malloc( out->m      * sizeof *out->face_tag);
  out->face_ptr[0]      = 0;

  out->dimensions[0]    = g.dims[0];
  out->dimensions[1]    = g.dims[1];
  out->dimensions[2]    = g.dims[2];
  out->number_of_faces  = 0;
  out->number_of_nodes  = 0;
  out->number_of_cells  = 0;

  out->node_coordinates = NULL;
  out->local_cell_index = malloc(nx*ny*nz * sizeof *out->local_cell_index);


  /* internal */
  int *intersections = malloc(BIGNUM* sizeof(*intersections));

  /* internal */
  int *work          = malloc(2* (2*nz+2)*   sizeof(*work));
  for(i=0; i<4*(nz+1); ++i)   work[i] = -1;

  /* Allocate space for cornerpoint numbers plus INT_MIN (INT_MAX) padding */
  int    *plist = malloc( 4*nx*ny*(2*nz+2) * sizeof *plist);


  /* Do actual work here:*/

  finduniquepoints(&g, plist, tolerance, out);

  free (zcorn);
  free (actnum);

  process_vertical_faces   (0, &intersections, plist, work, out);
  process_vertical_faces   (1, &intersections, plist, work, out);
  process_horizontal_faces (&intersections, plist, out);

  free (plist);
  free (work);

  compute_intersection_coordinates(intersections, out);

  free (intersections);


  /* Convert to local cell numbers in face_neighbors */
  int *ptr=out->face_neighbors;;
  for (i=0; i<out->number_of_faces*2; ++i, ++ptr){
    if (*ptr != -1){
      *ptr = out->local_cell_index[*ptr];
    }
  }

  /* Invert global-to-local map */
  int *global_cell_index = malloc(out->number_of_cells *
				  sizeof (*global_cell_index));
  for (i=0; i<nx*ny*nz; ++i){
    if(out->local_cell_index[i]!=-1){
      global_cell_index[out->local_cell_index[i]] = i;
    }
  }
  free(out->local_cell_index);
  out->local_cell_index = global_cell_index;

  /* HACK: undo sign reversal in ZCORN prepprocessing */
  for (i=2; i<3*out->number_of_nodes; i = i+3)
    out->node_coordinates[i]=sign*out->node_coordinates[i];

}

/*-------------------------------------------------------*/
void free_processed_grid(struct processed_grid *g)
{
  if( g ){
    free ( g->face_nodes       );
    free ( g->face_ptr         );
    free ( g->face_tag         );
    free ( g->face_neighbors   );
    free ( g->node_coordinates );
    free ( g->local_cell_index );
  }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
