/* Copyright 2013 Uni Research AS
 * This file is licensed under GPL3, see http://www.opm-project.org/
*/
#include <config.h>

/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE CompGeo2DTest
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

/* --- our own headers --- */
#include <algorithm>
#include <vector>
#include <opm/core/grid.h>
#include <opm/core/grid/cornerpoint_grid.h>  /* compute_geometry */

using namespace std;

/* Test properties on a grid that looks like this:
 *
 *       /------
 *      / |    |
 *     <  |    |
 *      \ |    |
 *       \------
 *
 */
struct BackspaceGrid {
   UnstructuredGrid* g;

   BackspaceGrid()
   {
      g = allocate_grid(
               2,     /* ndims */
               2,     /* number of cells */
               6,     /* number of edges */
               6*2,   /* six edges with two nodes each */
               3+4,   /* one triangle and one quadrilateral */
               5);    /* five nodes total */
      /* nodes */
      double nodes[] = {
         -1.,  0.,      /* node 0 */
          0.,  1.,      /* node 1 */
          1.,  1.,      /* node 2 */
          1., -1.,      /* node 3 */
          0., -1.,      /* node 4 */
      };
      copy(nodes, nodes+sizeof(nodes)/sizeof(nodes[0]), g->node_coordinates);
      /* edges */
      int edges[] = {
         0, 1,          /* edge 0 */
         1, 2,          /* edge 1 */
         2, 3,          /* edge 2 */
         3, 4,          /* edge 3 */
         4, 0,          /* edge 4 */
         4, 1,          /* edge 5 */
      };
      copy(edges, edges+sizeof(edges)/sizeof(edges[0]), g->face_nodes);
      /* starting index in map for each edge */
      int edge_pos[] = {
         0, 2, 4, 6, 8, 10, 12,
      };
      copy(edge_pos, edge_pos+sizeof(edge_pos)/sizeof(edge_pos[0]), g->face_nodepos);
      /* topology, clock-wise ordering */
      int neighbours[] = {
         -1, 0,   /* edge 0, between boundary and cell 0 */
         -1, 1,   /* edge 1, between boundary and cell 1 */
         -1, 1,   /* edge 2, between boundary and cell 1 */
         -1, 1,   /* edge 3, between boundary and cell 1 */
         -1, 0,   /* edge 4, between boundary and cell 0 */
          0, 1,   /* edge 5, between cell 0 and cell 1 */
      };
      copy(neighbours, neighbours+sizeof(neighbours)/sizeof(neighbours[0]), g->face_cells);
      /* cells */
      int cells[] = {
         0, 5, 4,    /* cell 0, clockwise */
         1, 2, 3, 5, /* cell 1, clockwise */
      };
      copy(cells, cells+sizeof(cells)/sizeof(cells[0]), g->cell_faces);
      /* starting index in map for each cell */
      int cell_pos[] = {
         0, 3, 7,
      };
      copy(cell_pos, cell_pos+sizeof(cell_pos)/sizeof(cell_pos[0]), g->cell_facepos);
      /* everything interesting actually happens here, for all tests */
      compute_geometry(g);
   }

   ~BackspaceGrid()
   {
      destroy_grid(g);
   }
};

BOOST_FIXTURE_TEST_SUITE(CompGeo2D, BackspaceGrid)

BOOST_AUTO_TEST_CASE(edgeMidpoints)
{
   double midpoints[] = {
      -0.5,  0.5,    /* edge 0 */
       0.5,  1.,     /* edge 1 */
       1.,   0.,     /* edge 2 */
       0.5, -1.,     /* edge 3 */
      -0.5, -0.5,    /* edge 4 */
       0.,   0.,     /* edge 5 */
   };
   BOOST_REQUIRE (sizeof(midpoints)/sizeof(midpoints[0]) ==
                  g->number_of_faces * g->dimensions);
   for (int edge = 0; edge < g->number_of_faces; ++edge)
   {
      for (int dim = 0; dim < g->dimensions; ++dim)
      {
         BOOST_REQUIRE_CLOSE (g->face_centroids[edge*g->dimensions+dim],
               midpoints[edge*g->dimensions+dim], 0.001);
      }
   }
}

BOOST_AUTO_TEST_CASE(edgeNormals)
{
   /* inward normal when doing clockwise enumeration */
   double normals[] = {
       1., -1.,      /* edge 0 */
       0., -1.,      /* edge 1 */
      -2.,  0.,      /* edge 2 */
       0.,  1.,      /* edge 3 */
       1.,  1.,      /* edge 4 */
       2.,  0.,      /* edge 5 */
   };
   BOOST_REQUIRE (sizeof(normals)/sizeof(normals[0]) ==
                  g->number_of_faces * g->dimensions);
   for (int edge = 0; edge < g->number_of_faces; ++edge)
   {
      for (int dim = 0; dim < g->dimensions; ++dim)
      {
         BOOST_REQUIRE_CLOSE (g->face_normals[edge*g->dimensions+dim],
               normals[edge*g->dimensions+dim], 0.001);
      }
   }
}

BOOST_AUTO_TEST_CASE(edgeLengths)
{
   double lengths[] = {
      1.4142,     /* edge 0 */
      1.,         /* edge 1 */
      2.,         /* edge 2 */
      1.,         /* edge 3 */
      1.4142,     /* edge 4 */
      2.,         /* edge 5 */
   };
   BOOST_REQUIRE (sizeof(lengths)/sizeof(lengths[0]) == g->number_of_faces);
   for (int edge = 0; edge < g->number_of_faces; ++edge)
   {
      BOOST_REQUIRE_CLOSE (g->face_areas[edge], lengths[edge], 0.001);
   }
}

BOOST_AUTO_TEST_CASE(cellAreas)
{
   double areas[] = {
      1., 2.,
   };
   BOOST_REQUIRE (sizeof(areas)/sizeof(areas[0]) == g->number_of_cells);
   for (int cell = 0; cell < g->number_of_cells; ++cell)
   {
      BOOST_REQUIRE_CLOSE (g->cell_volumes[cell], areas[cell], 0.001);
   }
}

BOOST_AUTO_TEST_CASE(cellCenters)
{
   double cellCenters[] = {
      -1./3., 0.0,     /* cell 0 */
       1./2., 0.0,     /* cell 1 */
   };
   BOOST_REQUIRE (sizeof(cellCenters)/sizeof(cellCenters[0]) ==
                  g->number_of_cells * g->dimensions);
   for (int cell = 0; cell < g->number_of_cells; ++cell)
   {
      for (int dim = 0; dim < g->dimensions; ++dim)
      {
         BOOST_REQUIRE_CLOSE (g->cell_centroids[cell*g->dimensions+dim],
               cellCenters[cell*g->dimensions+dim], 0.001);
      }
   }
}

BOOST_AUTO_TEST_SUITE_END()
