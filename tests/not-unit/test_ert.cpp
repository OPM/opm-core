/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "config.h"
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid/cart_grid.h>
#include <opm/core/grid.h>
#include <cstdio>
#include <boost/scoped_ptr.hpp>
#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>


#include <ert/ecl/ecl_grid.h>

using namespace std;

#if 0
static
void cell_nodes(const UnstructuredGrid * c_grid , int cell , std::vector<int>& nodes) {
  int face_offset = c_grid->cell_facepos[cell];
  int num_faces   = c_grid->cell_facepos[cell + 1] - face_offset;
  
  nodes.clear();
  //printf("cell: %d \n",cell);
  for (int iface = 0; iface < num_faces; iface++) {
    int face = c_grid->cell_faces[ face_offset + iface];
    //printf("face[%d] = %d \n",iface , face );
    
    {
      int node_offset = c_grid->face_nodepos[ face ];
      int num_nodes   = c_grid->face_nodepos[ face + 1] - node_offset;
      for (int inode = 0; inode < num_nodes; inode++) {
        int node = c_grid->face_nodes[ inode + node_offset ];
        //printf("   node[%d] = %d \n",inode , node);
        nodes.push_back( node );
      }
    }
  }

  {
    /*for (int i =0; i < nodes.size(); i++)
      std::cout << nodes[i] << " ";
    std::cout << "\n";
    */
    sort( nodes.begin() , nodes.end());
    /*for (int i =0; i < nodes.size(); i++)
      std::cout << nodes[i] << " ";
    std::cout << "\n";
    */
    unique( nodes.begin() , nodes.end() );
    /*for (int i =0; i < nodes.size(); i++)
      std::cout << nodes[i] << " ";
    std::cout << "\n";
    */
    nodes.resize( 8 );
    /*for (int i =0; i < nodes.size(); i++)
      std::cout << nodes[i] << " ";
    std::cout << "\n";
    */
  }
}
#endif

/*
static
void eclExport(Opm::GridManager& grid) {
  const UnstructuredGrid * c_grid = grid.c_grid();
  
  printf("dimensions         : %d \n",c_grid->dimensions);
  printf("number of cells    : %d \n",c_grid->number_of_cells);
  printf("number of nodes    : %d \n",c_grid->number_of_nodes);
  printf("number of faces    : %d \n",c_grid->number_of_faces);
  printf("length(face_nodes) : %d \n",c_grid->face_nodepos[ c_grid->number_of_faces ]);
  printf("cartdims         : %d %d %d \n",
         c_grid->cartdims[0] , 
         c_grid->cartdims[1] , 
         c_grid->cartdims[2]); 

  printf("global_cell      : %d %d %d %d %d\n",
         c_grid->global_cell[0] , 
         c_grid->global_cell[1] , 
         c_grid->global_cell[2] , 
         c_grid->global_cell[3] , 
         c_grid->global_cell[4]);

  {
    std::vector<int> nodes;
    cell_nodes( c_grid , 10 , nodes );
    cell_nodes( c_grid , 15 , nodes );
    cell_nodes( c_grid , 20 , nodes );
    cell_nodes( c_grid , 25 , nodes );
  }
 
  {
    ecl_grid_type * ecl_grid;
    int num_coords  = c_grid->number_of_cells;
    int coords_size = 4;
    int nx          = c_grid->cartdims[0];
    int ny          = c_grid->cartdims[1];
    int nz          = c_grid->cartdims[2]; 
    
    int   ** coords;
    float ** corners;
    // float  * mapaxes = NULL;
    std::vector<int> nodes;
    
    corners = (float **) malloc( num_coords * sizeof * corners );
    coords  = (int **) malloc( num_coords * sizeof * coords );

    {
      int c;
      for (c=0; c < num_coords; c++) {
        corners[c] = (float *) malloc( 24 * sizeof * corners[c] );
        coords[c]  = (int *)   malloc( coords_size * sizeof * coords[c] );
      }
      
      for (c=0; c < num_coords; c++) {
        cell_nodes( c_grid , c , nodes );
        for (int p=0; p < 8; p++) {
          int n = nodes[p];
          for (int d=0; d < c_grid->dimensions; d++) 
            corners[c][3*p + d] = c_grid->node_coordinates[ c_grid->dimensions * n + d ];
        }

        {
          int i,j,k;
          {
            int g = c_grid->global_cell[ c ];
            k =  g / nx*ny; g -= k * nx*ny;
            j =  g / nx;    g -= j * nx; 
            i  = g;
          }

          coords[c][0] = i + 1;
          coords[c][1] = j + 1;
          coords[c][2] = k + 1;
          coords[c][3] = c_grid->global_cell[ c ] + 1;
        }
      }
    }

    ecl_grid = ecl_grid_alloc( "/private/joaho/ERT/NR/libenkf/src/Gurbat/EXAMPLE_01_BASE.EGRID" );
    printf("Grid loaded ... \n");
    ecl_grid_free( ecl_grid );
    
    printf("Grid discarded ... \n");

    ecl_grid = ecl_grid_alloc_GRID_data( num_coords , nx , ny , nz , coords_size , coords , corners , NULL );
    ecl_grid_fwrite_GRID( ecl_grid , "/tmp/test.GRID" );

    {
      FILE * stream = fopen( "/tmp/test.grdecl" , "w");
      ecl_grid_fprintf_grdecl( ecl_grid , stream );
      fclose( stream );
    }

    ecl_grid_free( ecl_grid );
    
    {
      for (int c=0; c < num_coords; c++) {
        free(corners[c]);
        free(coords[c]);
      }
    }
    free( corners );
    free( coords );
  }
}
*/


/*
  
  #ifdef HAVE_ERT
ecl_grid_type * create_ecl_grid( const struct UnstructuredGrid * g) {
  int num_coords  = g->number_of_cells;
  int nx          = g->cartdims[0];
  int ny          = g->cartdims[1];
  int nz          = g->cartdims[2];
  int coords_size = 4;
  int   ** coords;
  float ** corners;
  float  * mapaxes = NULL;
  
  corners = malloc( num_coords * sizeof * corners );
  coords  = malloc( num_coords * sizeof * coords );

  {
    for (int c=0; c < num_coords; c++) {
      corners[c] = malloc( 24 * sizeof * corners[0] );
      coords[c]  = malloc( coords_size * sizeof * coords[0] );
    }
  }

  {
    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
          int global_index = i + j*nx + k*nx*ny;

          coords[global_index][0] = i;
          coords[global_index][1] = j;
          coords[global_index][2] = k;
          coords[global_index][3] = 1;
          
        }
      }
    }
  }
  {
    for (int c=0; c < num_coords; c++) {
      free(corners[c]);
      free(coords[c]);
    }
  }
  free( corners );
  free( coords );
}
#endif


  */

// struct grdecl         : opm/core/grid/cpgpreprocess/preprocess.h
// struct processes_grid : opm/core/grid/cpgpreprocess/preprocess.h



int main(int /*argc*/ , char **argv)
{
  std::string filename( argv[1] );
  boost::scoped_ptr<Opm::GridManager> grid;
  boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
  Opm::EclipseGridParser eclParser(filename , false);

  //eclParser.saveEGRID_INIT("/tmp" , "OPM" );

  grid.reset(new Opm::GridManager(eclParser));
  
  props.reset(new Opm::IncompPropertiesFromDeck(eclParser , *grid->c_grid()));
}
