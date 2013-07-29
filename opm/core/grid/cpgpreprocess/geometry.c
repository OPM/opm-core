/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */
#include "config.h"
#include <math.h>
#include <stdio.h>
#include "geometry.h"
#include <assert.h>

/* ------------------------------------------------------------------ */
static void
cross(const double u[3], const double v[3], double w[3])
/* ------------------------------------------------------------------ */
{
   w[0] = u[1]*v[2]-u[2]*v[1];
   w[1] = u[2]*v[0]-u[0]*v[2];
   w[2] = u[0]*v[1]-u[1]*v[0];
}

/* ------------------------------------------------------------------ */
static double
norm(const double w[3])
/* ------------------------------------------------------------------ */
{
   return sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
}


/* ------------------------------------------------------------------ */
static void
compute_face_geometry_3d(double *coords, int nfaces,
                         int *nodepos, int *facenodes, double *fnormals,
                         double *fcentroids, double *fareas)
/* ------------------------------------------------------------------ */
{

   /* Assume 3D for now */
   const int ndims = 3;
   int f;
   double x[3];
   double u[3];
   double v[3];
   double w[3];

   int i,k;
   int node;

   double cface[3]  = {0};
   double n[3]  = {0};
   double twothirds = 0.666666666666666666666666666667;
   double a;
   int    num_face_nodes;
   double area;
   /*#pragma omp parallel for  */

   /*#pragma omp parallel for shared(fnormals,fcentroids,fareas)*/
#pragma omp parallel for default(none) 			   \
    private(f,x,u,v,w,i,k,node,cface,n,a,num_face_nodes,area)		\
    shared(fnormals,fcentroids,fareas  \
	   ,coords, nfaces, nodepos, facenodes) \
    firstprivate(ndims, twothirds)
   for (f=0; f<nfaces; ++f)
   {
      for(i=0; i<ndims; ++i) x[i] = 0.0;
      for(i=0; i<ndims; ++i) n[i] = 0.0;
      for(i=0; i<ndims; ++i) cface[i] = 0.0;

      /* average node */
      for(k=nodepos[f]; k<nodepos[f+1]; ++k)
      {
         node = facenodes[k];
         for (i=0; i<ndims; ++i) x[i] += coords[3*node+i];
      }
      num_face_nodes = nodepos[f+1] - nodepos[f];
      for(i=0; i<ndims; ++i) x[i] /= num_face_nodes;



      /* compute first vector u (to the last node in the face) */
      node = facenodes[nodepos[f+1]-1];
      for(i=0; i<ndims; ++i) u[i] = coords[3*node+i] - x[i];

      area=0.0;
      /* Compute triangular contrib. to face normal and face centroid*/
      for(k=nodepos[f]; k<nodepos[f+1]; ++k)
      {


         node = facenodes[k];
         for (i=0; i<ndims; ++i) v[i] = coords[3*node+i] - x[i];

         cross(u,v,w);
         a = 0.5*norm(w);
         area += a;
	 /*         if(!(a>0))
         {
            fprintf(stderr, "Internal error in compute_face_geometry.");
         }
	 */
         /* face normal */
         for (i=0; i<ndims; ++i) n[i] += w[i];

         /* face centroid */
         for (i=0; i<ndims; ++i)
            cface[i] += a*(x[i]+twothirds*0.5*(u[i]+v[i]));

         /* Store v in u for next iteration */
         for (i=0; i<ndims; ++i) u[i] = v[i];
      }

      /* Store face normal and face centroid */
      for (i=0; i<ndims; ++i)
      {
         /* normal is scaled with face area */
         fnormals  [3*f+i] = 0.5*n[i];
         fcentroids[3*f+i] = cface[i]/area;
      }
      fareas[f] = area;
   }
}

/* ------------------------------------------------------------------ */
static void
compute_edge_geometry_2d(
      /* in  */ double *node_coords,
      /* in  */ int     num_edges,
      /* in  */ int    *edge_node_pos,
      /* in  */ int    *edge_nodes,
      /* out */ double *edge_normals,
      /* out */ double *edge_midpoints,
      /* out */ double *edge_lengths)
{
   const int num_dims = 2;

   /* offsets to each of the nodes in a compacted edge */
   const int a_ofs = 0;
   const int b_ofs = 1;

   /* offsets to each dimension is a compacted point */
   const int x_ofs = 0;
   const int y_ofs = 1;

   int edge;                     /* edge index       */
   int a_nod, b_nod;             /* node indices     */
   double a_x, a_y, b_x, b_y;    /* node coordinates */
   double v_x, v_y;              /* vector elements  */

   /* decompose each edge into a tuple (a,b) between two points and
    * compute properties for that face. hopefully the host has enough
    * cache pages to keep both input and output at the same time, and
    * registers for all the local variables */
   for (edge = 0; edge < num_edges; ++edge)
   {
      /* an edge in 2D can only have starting and ending point
       * check that there are exactly two nodes till the next edge */
      assert (edge_node_pos[edge + 1] - edge_node_pos[edge] == num_dims);

      /* get the first and last point on the edge */
      a_nod = edge_nodes[edge_node_pos[edge] + a_ofs];
      b_nod = edge_nodes[edge_node_pos[edge] + b_ofs];

      /* extract individual coordinates for the points */
      a_x = node_coords[a_nod * num_dims + x_ofs];
      a_y = node_coords[a_nod * num_dims + y_ofs];
      b_x = node_coords[b_nod * num_dims + x_ofs];
      b_y = node_coords[b_nod * num_dims + y_ofs];

      /* compute edge center -- average of node coordinates */
      edge_midpoints[edge * num_dims + x_ofs] = (a_x + b_x) * 0.5;
      edge_midpoints[edge * num_dims + y_ofs] = (a_y + b_y) * 0.5;

      /* vector from first to last point */
      v_x = b_x - a_x;
      v_y = b_y - a_y;

      /* two-dimensional (unary) cross product analog that makes the
       * "triple" (dot-cross) product zero, i.e. it's a normal; the
       * direction of this vector is such that it will be pointing
       * inwards when enumerating nodes clock-wise */
      edge_normals[edge * num_dims + x_ofs] = +v_y;
      edge_normals[edge * num_dims + y_ofs] = -v_x;

      /* Euclidian norm in two dimensions is magnitude of edge */
      edge_lengths[edge] = sqrt(v_x*v_x + v_y*v_y);
   }
}

/* ------------------------------------------------------------------ */
void
compute_face_geometry(int ndims, double *coords, int nfaces,
                      int *nodepos, int *facenodes, double *fnormals,
                      double *fcentroids, double *fareas)
/* ------------------------------------------------------------------ */
{
   if (ndims == 3)
   {
      compute_face_geometry_3d(coords, nfaces, nodepos, facenodes,
                               fnormals, fcentroids, fareas);
   }
   else if (ndims == 2)
   {
      /* two-dimensional interfaces are called 'edges' */
      compute_edge_geometry_2d(coords, nfaces, nodepos, facenodes,
                               fnormals, fcentroids, fareas);
   }
   else
   {
      assert(0);
   }
}


/* ------------------------------------------------------------------ */
static void
compute_cell_geometry_3d(double *coords,
                         int *nodepos, int *facenodes, int *neighbors,
                         double *fnormals,
                         double *fcentroids,
                         int ncells, int *facepos, int *cellfaces,
                         double *ccentroids, double *cvolumes)
/* ------------------------------------------------------------------ */
{
   const int ndims = 3;
   int i,k, f,c;
   int face,node;
   double x[3];
   double u[3];
   double v[3];
   double w[3];
   double xcell[3];
   double ccell[3];
   double cface[3]  = {0};
   int num_faces;
   double volume;
   double tet_volume, subnormal_sign;
   double twothirds = 0.666666666666666666666666666667;
#pragma omp parallel for default(none)     \
    private(i,k,f,c,face,node,x,u,v,w,xcell				\
	    ,ccell ,cface,num_faces,volume, tet_volume, subnormal_sign) \
   shared(coords,nodepos,facenodes,neighbors,				\
	  fnormals,fcentroids,facepos,cellfaces,ccentroids,cvolumes)	\
    firstprivate(ncells,ndims,twothirds)
   for (c=0; c<ncells; ++c)
   {


      for(i=0; i<ndims; ++i) xcell[i] = 0.0;
      for(i=0; i<ndims; ++i) ccell[i] = 0.0;


      /*
       * Approximate cell center as average of face centroids
       */
      for(f=facepos[c]; f<facepos[c+1]; ++f)
      {
         face = cellfaces[f];
         for (i=0; i<ndims; ++i) xcell[i] += fcentroids[3*face+i];
      }
      num_faces = facepos[c+1] - facepos[c];

      for(i=0; i<ndims; ++i) xcell[i] /= num_faces;



      /*
       * For all faces, add tetrahedron's volume and centroid to
       * 'cvolume' and 'ccentroid'.
       */
      volume=0.0;
      for(f=facepos[c]; f<facepos[c+1]; ++f)
      {
         int num_face_nodes;

         for(i=0; i<ndims; ++i) x[i] = 0.0;
         for(i=0; i<ndims; ++i) cface[i] = 0.0;

         face = cellfaces[f];

         /* average face node x */
         for(k=nodepos[face]; k<nodepos[face+1]; ++k)
         {
            node = facenodes[k];
            for (i=0; i<ndims; ++i) x[i] += coords[3*node+i];
         }
         num_face_nodes = nodepos[face+1] - nodepos[face];
         for(i=0; i<ndims; ++i) x[i] /= num_face_nodes;



         /* compute first vector u (to the last node in the face) */
         node = facenodes[nodepos[face+1]-1];
         for(i=0; i<ndims; ++i) u[i] = coords[3*node+i] - x[i];


         /* Compute triangular contributions to face normal and face centroid */
         for(k=nodepos[face]; k<nodepos[face+1]; ++k)
         {

            node = facenodes[k];
            for (i=0; i<ndims; ++i) v[i] = coords[3*node+i] - x[i];

            cross(u,v,w);



            tet_volume = 0.0;
            for(i=0; i<ndims; ++i){
		tet_volume += w[i]*(x[i]-xcell[i]);
	    }
            tet_volume *= 0.5 / 3;

	    subnormal_sign=0.0;
	    for(i=0; i<ndims; ++i){
		subnormal_sign += w[i]*fnormals[3*face+i];
	    }

	    if(subnormal_sign < 0.0){
		tet_volume =- tet_volume;
	    }
	    if(!(neighbors[2*face+0]==c)){
		tet_volume = -tet_volume;
	    }
	    volume += tet_volume;
            /* face centroid of triangle  */
            for (i=0; i<ndims; ++i) cface[i] = (x[i]+(twothirds)*0.5*(u[i]+v[i]));

            /* Cell centroid */
            for (i=0; i<ndims; ++i) ccell[i] += tet_volume * 3/4.0*(cface[i] - xcell[i]);


            /* Store v in u for next iteration */
            for (i=0; i<ndims; ++i) u[i] = v[i];
         }
      }
      for (i=0; i<ndims; ++i) ccentroids[3*c+i] = xcell[i] + ccell[i]/volume;
      cvolumes[c] = volume;
   }
}

/* ------------------------------------------------------------------ */
static void
compute_cell_geometry_2d(
      /* in  */ double *node_coords,
      /* in  */ int     *edge_node_pos,
      /* in  */ int    *edge_nodes,
      /* in  */ double *edge_midpoints,
      /* in  */ int     num_cells,
      /* in  */ int    *cell_edge_pos,
      /* in  */ int    *cell_edges,
      /* out */ double *cell_centers,
      /* out */ double *cell_areas)
{
   const int num_dims = 2;

   /* offsets to each of the nodes in a compacted edge */
   const int a_ofs = 0;
   const int b_ofs = 1;

   /* offsets to each dimension is a compacted point */
   const int x_ofs = 0;
   const int y_ofs = 1;

   int cell;            /* cell index */
   int num_nodes;       /* number of vertices in current cell */
   int edge_ndx;        /* relative edge index within cell */
   int edge;            /* absolute cell index */
   double center_x;     /* x-coordinate for cell barycenter */
   double center_y;     /* y-coordinate for cell barycenter */
   double area;         /* (accumulated) cell area */
   int a_nod, b_nod;    /* node indices for edge start and end points */
   double a_x, a_y,
          b_x, b_y;     /* vectors from center to edge points */

   for (cell = 0; cell < num_cells; ++cell)
   {
      /* since the cell is a closed polygon, each point serves as the starting
       * point of one edge and the ending point of another; thus there is as
       * many vertices as there are edges */
      num_nodes = cell_edge_pos[cell + 1] - cell_edge_pos[cell];

      /* to enumerate all vertices of a cell, we would have to expand the
       * edges and then remove duplicates. however, the centroid of each
       * edge contains half of the two vertices that are incident on it. if
       * we instead sum all the face centroids, we get the sum of all the
       * vertices */
      center_x = 0.;
      center_y = 0.;
      for (edge_ndx = cell_edge_pos[cell];
           edge_ndx < cell_edge_pos[cell + 1]; ++edge_ndx)
      {
         edge = cell_edges[edge_ndx];
         center_x += edge_midpoints[edge * num_dims + x_ofs];
         center_y += edge_midpoints[edge * num_dims + y_ofs];
      }
      center_x /= (double) num_nodes;
      center_y /= (double) num_nodes;
      cell_centers[cell * num_dims + x_ofs] = center_x;
      cell_centers[cell * num_dims + y_ofs] = center_y;

      /* triangulate the polygon by introducing the cell center and then new
       * internal edges from this center to the vertices. the total area of
       * the cell is the sum of area of these sub-triangles */
      area = 0.;
      for (edge_ndx = cell_edge_pos[cell];
           edge_ndx < cell_edge_pos[cell + 1]; ++edge_ndx)
      {
         /* indirect lookup of edge index (from array that contains all the
          * edge indices for a certain cell) */
         edge = cell_edges[edge_ndx];

         /* get the first and last point on the edge */
         a_nod = edge_nodes[edge_node_pos[edge] + a_ofs];
         b_nod = edge_nodes[edge_node_pos[edge] + b_ofs];

         /* vector from center to each of the nodes */
         a_x = node_coords[a_nod * num_dims + x_ofs] - center_x;
         a_y = node_coords[a_nod * num_dims + y_ofs] - center_y;
         b_x = node_coords[b_nod * num_dims + x_ofs] - center_x;
         b_y = node_coords[b_nod * num_dims + y_ofs] - center_y;

         /* two-dimensional (binary) cross product analog that has length
          * equal to the parallelogram spanned by the two vectors (but which
          * is a scalar). the sign tells us the orientation between the nodes
          * a and b, but we are not interested in that, just the area */
         area += fabs(a_x * b_y - a_y * b_x);
      }
      /* we summed parallelograms which are twice the size of the triangles
       * that make up the cell; divide out the half for all terms here */
      area *= 0.5;
      cell_areas[cell] = area;
   }
}

/* ------------------------------------------------------------------ */
void
compute_cell_geometry(int ndims, double *coords,
                      int *nodepos, int *facenodes, int *neighbors,
                      double *fnormals,
                      double *fcentroids,
                      int ncells, int *facepos, int *cellfaces,
                      double *ccentroids, double *cvolumes)
/* ------------------------------------------------------------------ */
{
   if (ndims == 3)
   {
      compute_cell_geometry_3d(coords, nodepos, facenodes,
                               neighbors, fnormals, fcentroids, ncells,
                               facepos, cellfaces, ccentroids, cvolumes);
   }
   else if (ndims == 2)
   {
      compute_cell_geometry_2d(coords, nodepos, facenodes, fcentroids,
                               ncells, facepos, cellfaces, ccentroids,
                               cvolumes);
   }
   else
   {
      assert(0);
   }
}
