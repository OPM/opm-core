#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "preprocess.h"
#include "sparsetable.h"
#include "uniquepoints.h"
#include "facetopology.h"

/* No checking of input arguments in this code! */
#define min(i,j) ((i)<(j) ? (i) : (j))
#define max(i,j) ((i)>(j) ? (i) : (j))



static void interpolate_pillar(const double *coord, double *pt)
{
  double a = (pt[2]-coord[2])/(coord[5]-coord[2]);
  pt[0]       = coord[0] + a*(coord[3]-coord[0]);
  pt[1]       = coord[1] + a*(coord[4]-coord[1]);
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


/*-------------------------------------------------------*/
void free_processed_grid(struct processed_grid *g)
{
  if( g ){
    free ( g->face_nodes       ); 
    free ( g->face_ptr         );
    free ( g->face_neighbors   );
    free ( g->node_coordinates );
    free ( g->local_cell_index );
  }
}





/*-------------------------------------------------------*/
void compute_cell_index(const int dims[3], int i, int j, int *neighbors, int len)   
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


/* Ensure there's sufficient memory */
/*-------------------------------------------------------*/
int checkmemeory(int nz, sparse_table_t *ftab, int **neighbors, int **intersections)
{

  /* Ensure there is enough space */
  int r = (2*nz+2)*(2*nz+2);
  int m = ftab->m;
  int n = ftab->n;
  if(ftab->position +  r> m){
    m += max(m*0.5, 2*r);
  }
  if (ftab->ptr[ftab->position] + 6*r > n){
    n += max(n*0.5, 12*r);
  }
      
  if (m != ftab->m){
    void *p1 = realloc(*neighbors,     2*m * sizeof(**neighbors));
    void *p2 = realloc(*intersections, 4*m * sizeof(**neighbors));
    if (p1 && p2){
      *neighbors     = p1;
      *intersections = p2;
    }else{
      return 0;
    }
  }

  if (m != ftab->m || n != ftab->n){
    void *p = realloc_sparse_table(ftab, m, n, sizeof(int));
    if (p){
      ftab = p;
    }else{
      return 0;
    }
  }

  return 1;
}

/*-------------------------------------------------------*/
void process_vertical_faces(const int dims[3], int direction, sparse_table_t *ftab,  
			    int **neighbors, int **intersections, int *npoints, 
			    int npillarpoints, int *plist, int *work)
{
  /* direction == 0 : constant i-faces
     direction == 1 : constant j-faces */

  int i,j;
  int *cornerpts[4];

  /* constant i- or j-faces */
  for (j=0; j<dims[1]+direction; ++j) {
    for (i=0; i<dims[0]+1-direction; ++i){

      if (!checkmemeory(dims[2], ftab, neighbors, intersections)){
	fprintf(stderr, "Could not allocat enough space\n");
	exit(1);      
      }

      
      /* Vectors of point numbers */
      int d[] = {2*dims[0], 2*dims[1], 2+2*dims[2]};
      igetvectors(d, 2*i+direction, 2*j+1-direction, plist, cornerpts);
      
      if(direction==1){
	/* swap */
	int *tmp     = cornerpts[2];
	cornerpts[2] = cornerpts[1];
	cornerpts[1] = tmp;
      }

      int startface = ftab->position;
      int num_intersections = *npoints - npillarpoints;
      findconnections(2*dims[2]+2, cornerpts, npoints, 
		      *intersections+4*num_intersections, 
		      *neighbors, work, ftab);
      
	
 

      int *ptr = *neighbors + 2*startface;
      int len  = 2*ftab->position - 2*startface;

      compute_cell_index(dims, i-1+direction, j-direction, ptr,   len);
      compute_cell_index(dims, i,   j, ptr+1, len);
    }
  }
}

int global_cell(const int dims[3], int i, int j, int k)
{
  return i+dims[0]*(j+dims[1]*k);
}
#include <limits.h>
/*-------------------------------------------------------*/
void process_horizontal_faces(const struct grdecl *g, 
			      int                 *cell_index,
			      int                 *ncells,
			      sparse_table_t      *faces, 
			      int                 **neighbors, 
			      int                 **intersections, 
			      int                 *plist)
{
  int i,j,k;
  int *cell  = cell_index;
  int cellno = 0;
  
  /* compute constant-k faces */
  for(j=0; j<g->dims[1]; ++j){
    for (i=0; i<g->dims[0]; ++i) {


      if (!checkmemeory(g->dims[2], faces, neighbors, intersections)){
	fprintf(stderr, "Could not allocat enough space\n");
	exit(1);
      }

      
      int *f       = (int*)faces->data;
      int *newface = (int*)faces->data + faces->ptr[faces->position];
      int *newneighbors = *neighbors   + 2*faces->position;

 
      /* Vectors of point numbers */
      int *c[4];
      int  d[] = {2*g->dims[0], 2*g->dims[1], 2+2*g->dims[2]};
      igetvectors(d, 2*i+1, 2*j+1, plist, c);
      
      int previous_global_cell = -1;
      for (k = 1; k<g->dims[2]*2+1; ++k){

      	
	/* Skip if space between face k and face k+1 is collapsed. */
	/* Note that inactive cells (with ACTNUM==0) have all been  */
	/* collapsed in finduniquepoints.                           */
	if (c[0][k] == c[0][k+1] && c[1][k] == c[1][k+1] &&
	    c[2][k] == c[2][k+1] && c[3][k] == c[3][k+1]){
	  
	  if (k%2) cell[global_cell(g->dims, i,j,(k-1)/2)] = -1;  
	  

	}
	else{

	  if (k%2){
	    /* Add face */
	    *newface++ = c[0][k];
	    *newface++ = c[1][k];
	    *newface++ = c[3][k];
	    *newface++ = c[2][k];
	    
	    faces->ptr[++faces->position] = newface - f;
	    *newneighbors++ = previous_global_cell;
	    *newneighbors++ = previous_global_cell = global_cell(g->dims, i,j,(k-1)/2);
	     
	    cell[global_cell(g->dims, i,j,(k-1)/2)] = cellno++;
	    
	  }
	  else{
	    if (previous_global_cell != -1){
	      /* Add face */
	      *newface++ = c[0][k];
	      *newface++ = c[1][k];
	      *newface++ = c[3][k];
	      *newface++ = c[2][k];
	      
	      faces->ptr[++faces->position] = newface - f;
	      *newneighbors++ = previous_global_cell;
	      *newneighbors++ = previous_global_cell = -1;
	    }
	  }
	}
      }
    }
  }
  *ncells = cellno;
}


static void approximate_intersection_pt(int     *intersection, 
					double  *coordinates, double *pt)
{
  double z0 = coordinates[3*intersection[0]+2];
  double z1 = coordinates[3*intersection[1]+2];
  double z2 = coordinates[3*intersection[2]+2];
  double z3 = coordinates[3*intersection[3]+2];

  double alpha = (z2-z0)/(z1-z0 - (z3-z2));

  pt[0] = coordinates[3*intersection[0]+0]*(1.0-alpha) 
         +coordinates[3*intersection[1]+0]*     alpha;
  pt[1] = coordinates[3*intersection[0]+1]*(1.0-alpha) 
         +coordinates[3*intersection[1]+1]*     alpha;
  pt[2] = coordinates[3*intersection[0]+2]*(1.0-alpha) 
         +coordinates[3*intersection[1]+2]*     alpha;

}

static void compute_node_coordinates(const struct grdecl *g, 
				     double              *coordinates, 
				     int                 *intersections, 
				     sparse_table_t      *pillarz, 
				     int                 npillarpoints,
				     int                 npoints)
{
  int i,k;
  int nx = g->dims[0];
  int ny = g->dims[1];

  double       *pt          = coordinates;
  const double *c           = g->coord;

  /* Loop over pillars */
  int pillar = 0;
  for (i=0; i< (nx+1)*(ny+1); ++i){


    /* Loop over unique zcorn values - may be none */
    for (k=pillarz->ptr[pillar]; k<pillarz->ptr[pillar+1]; ++k) {
	
      /* Assign z-coordinate */
      pt[2] = ((double*)pillarz->data)[k];

      /* Compute x- and y- coordinate */
      interpolate_pillar(c, pt);
	
      pt += 3;
    }

    ++pillar;
    c += 6;
  }

  
  /* Append intersections */
  int *itsct = intersections;
  for (k=npillarpoints; k<npoints; ++k){

    approximate_intersection_pt(itsct, coordinates, pt);

    pt    += 3;
    itsct += 4;
  }
}


/* Gateway routine.              */
/*-------------------------------*/
void processGrdecl(const struct grdecl *g, double tol, struct processed_grid *out)
{
  int i;

  /* Code below assumes k index runs fastests, ie. that dimensions of
     table is permuted to (dims[2], dims[0], dims[1]) */

  int nx = g->dims[0]; 
  int ny = g->dims[1];
  int nz = g->dims[2];


  /* ztab->data may need extra space temporarily due to simple boundary treatement  */
  int            npillarpoints = 8*(nx+1)*(ny+1)*nz; 
  int            npillars      = (nx+1)*(ny+1);
  sparse_table_t *pillarz      = malloc_sparse_table(npillars, 
						     npillarpoints, 
						     sizeof(double));

  /* Allocate space for cornerpoint numbers plus INT_MIN (INT_MAX) padding */
  int    *plist = malloc( 4*nx*ny*(2*nz+2) * sizeof(int));

  /* Fill plist of cornerpoint numbers and ztab of unique zcorn values. */
  finduniquepoints(g, plist, pillarz);

  npillarpoints = pillarz->ptr[npillars];


  void *p = realloc_sparse_table (pillarz, npillars, npillarpoints, sizeof(double));
  if (p) {
    pillarz = p;
  }else{
    fprintf(stderr, "Could not reallocate space\n");
    exit(1);
  }



  /* Process face geometry and cell-face topology on constant-i pillar pairs */


  const int BIGNUM = 64;
  /* Unstructured storage of face nodes */
  sparse_table_t *faces = malloc_sparse_table(BIGNUM/3,
					     BIGNUM,
					     sizeof(int));

  int *neighbors     = malloc(BIGNUM* sizeof(*neighbors));
  int *intersections = malloc(BIGNUM* sizeof(*intersections));

  
  int *work          = malloc(2* (2*nz+2)*   sizeof(*work));
  for(i=0; i<2* (2*nz+2); ++i){
    work[i] = -1;
  }
  int npoints        = npillarpoints;
  
  faces->position = 0;
  faces->ptr[0]   = 0;

  /* faces with constant i */
  process_vertical_faces(g->dims, 0, faces,
			 &neighbors, &intersections, &npoints,
			 npillarpoints, plist, work);
  /* faces with constant j */
  process_vertical_faces(g->dims, 1, faces,  
			 &neighbors, &intersections, &npoints, 
			 npillarpoints, plist, work);
  
  free(work);
  

  int *cell_index = calloc(nx*ny*nz,sizeof(*cell_index));
  int ncells;
  
  process_horizontal_faces(g, cell_index, &ncells, faces, &neighbors, 
			   &intersections, plist);


  /* Convert to local cell numbers in face_neighbors */
  int *ptr=neighbors;;
  for (i=0; i<faces->position*2; ++i, ++ptr){
    if (*ptr != -1){
      *ptr = cell_index[*ptr];
    }
  }
 
  /* Invert global-to-local map */
  ptr = cell_index;
  for (i=0; i<nx*ny*nz; ++i){
    if(cell_index[i]!=-1){
      *ptr++ = i;
    }
  }

  /* Shrink memory allocated for cell_index */
  p = realloc(cell_index, (ptr-cell_index)*sizeof(*cell_index));
  if (p){
    cell_index = p;
  }
  else{
    fprintf(stderr, "Could not reallocate space for index map\n");
    exit(1);
  }


  /* compute node coordinates on pillars and new intersections */  
  double       *coordinates = malloc(3*npoints * sizeof(*coordinates));

  compute_node_coordinates(g, coordinates, intersections, pillarz, 
			   npillarpoints, npoints);  
  
  
   free_sparse_table(pillarz);

  free (intersections);

  free (plist);


  out->number_of_faces  = faces->position;
  out->face_nodes       = faces->data;
  out->face_ptr         = faces->ptr;
  out->face_neighbors   = neighbors;
  out->number_of_nodes  = npoints;
  out->node_coordinates = coordinates;
  out->number_of_cells  = ncells;
  out->local_cell_index = cell_index;
}
