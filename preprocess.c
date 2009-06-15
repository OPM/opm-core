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



static void interpolate_pillar(double *coord, double *x, double *y, double *z)
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
	neighbors[k] = i + dims[0]*(j + dims[1]*neighbors[k]);
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
  /* constant i-faces */
  for (j=0; j<dims[1]+direction; ++j) {
    for (i=0; i<dims[0]+1-direction; ++i){
      /* fprintf(stderr, "%p %p\n  ", neighbors, intersections); */
      if (!checkmemeory(dims[2], ftab, neighbors, intersections)){
	fprintf(stderr, "Could not allocat enough space\n");
	exit(1);      
      }
      /* fprintf(stderr, "%p %p\n\n", neighbors, intersections); */

      int startface = ftab->position;
      int num_intersections = *npoints - npillarpoints;

      /* Vectors of point numbers */
      int d[] = {2*dims[0], 2*dims[1], 2+2*dims[2]};
      igetvectors(d, 2*i+direction, 2*j+1-direction, plist, cornerpts);

      if(direction==1){
	/* swap */
	int *tmp     = cornerpts[2]; 
	cornerpts[2] = cornerpts[1];
	cornerpts[1] = tmp;
      }
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



void process_horizontal_faces(const struct grdecl *g, int *cell_index, 
			      sparse_table_t *faces, int **neighbors, 
			      int **intersections, int *plist)
{
  int nx = g->dims[0]; 
  int ny = g->dims[1];

  int i,j,k;
  int *cell = cell_index;
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
      int *newneighbors = *neighbors + 2*faces->position;


      
      /* Vectors of point numbers */
      int *c[4];
      int d[] = {2*g->dims[0], 2*g->dims[1], 2+2*g->dims[2]};
      igetvectors(d, 2*i+1, 2*j+1, plist, c);
      
      int kprev = -1;
      for (k = 1; k<g->dims[2]*2+1; ++k){

	/* Skip if space between face k and face k+1 is collapsed. */
	/* Note that inactive cells (with ACNUM==0) have all been  */
	/* collased in finduniquepoints.                           */
	if (c[0][k] == c[0][k+1] &&
	    c[1][k] == c[1][k+1] &&
	    c[2][k] == c[2][k+1] &&
	    c[3][k] == c[3][k+1]){
	  if (k%2)*cell++ = -1;
	  continue;
	}else{

	  /* Add face */
	  *newface++ = c[0][k];
	  *newface++ = c[1][k];
	  *newface++ = c[2][k];
	  *newface++ = c[3][k];
	  
	  faces->ptr[++faces->position] = newface - f;
	  *newneighbors++ = kprev != -1       ? i+nx*(j+ny*(kprev-1)/2) : -1;
	  *newneighbors++ = k != g->dims[2]*2 ? i+nx*(j+ny*(k    -1)/2) : -1;

	  kprev = k;  

	  if (k%2)*cell++ = cellno++;
	}
      }
    }
  }

}

/* Gateway routine.              */
/*-------------------------------*/
void processGrdecl(const struct grdecl *g, double tol, struct processed_grid *out)
{
  int i,k;

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


  fprintf(stderr, "process face geometry\n");
  /* Process face geometry and cell-face topology on constant-i pillar pairs */


  const int BIGNUM = 64;
  /* Unstructured storage of face nodes */
  sparse_table_t *faces = malloc_sparse_table(BIGNUM/3,
					     BIGNUM,
					     sizeof(int));

  int *neighbors     = malloc(BIGNUM* sizeof(*neighbors));
  int *intersections = malloc(BIGNUM* sizeof(*intersections));

  
  int *work          = calloc(2* (2*nz+2),   sizeof(*work));
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
  
  fprintf(stderr, "processing k-constant faces and cell numbering \n");
  fprintf(stderr, "Cell numbering is warped due to k running fastest\n");
  int *cell_index = malloc(nx*ny*nz*sizeof(*cell_index));

  process_horizontal_faces(g, cell_index, faces, &neighbors, &intersections, plist);


  /*                                                       */
  /*                                                       */
  /*                                                       */
  /*                 D E B U G    CODE                     */
  /*                                                       */
  /*                                                       */
  /*                                                       */

#if 1
      fprintf(stderr, "\nfaces\nnumfaces %d\n", faces->position);
      for (i=0; i<faces->position; ++i){
	for (k=faces->ptr[i]; k<faces->ptr[i+1]; ++k){
	  fprintf(stderr, "%d ", ((int*)faces->data)[k]);
	}
	fprintf(stderr, "\n");
      }


      fprintf(stderr, "\nneighbors\n");
      int *iptr = neighbors;
      for(k=0; k<faces->position; ++k){
	fprintf(stderr, " (%d %d)\n",  iptr[0], iptr[1]);
	++iptr;
	++iptr;
      }

      fprintf(stderr, "\nline intersections:\n");
      iptr = intersections;
      int numintersections = npoints - npillarpoints;
      for(k=0; k<numintersections; ++k){
	fprintf(stderr, " (%d %d %d %d)\n",  iptr[0], iptr[1], iptr[2], iptr[3]);
	iptr+=4;
      }


      /* fprintf(stderr, "number of cells : %d\n", nx*ny*nz); */
      /* for (k=0; k<nx*ny*nz; ++k){ */
      /* fprintf(stderr, "cell index %d is %d\n", k, cell_index[k]); */
      /* } */
#endif
  /* compute node coordinates */
  free_sparse_table(pillarz);

  /* compute intersections */
  free (intersections);

  /* compute local cell numbering */
  free (plist);


  out->number_of_faces  = faces->position;
  out->face_nodes       = faces->data;
  out->face_ptr         = faces->ptr;
  out->face_neighbors   = neighbors;
  out->number_of_nodes  = 0;
  out->node_coordinates = NULL;
  out->number_of_cells  = 0;
  out->local_cell_index = NULL;

}
