//===========================================================================
//
// File: GridCplusplus.hpp
//
// Created: Fri Sep 24 11:01:04 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This file is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the code.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SINTEF_GRIDCPLUSPLUS_HEADER
#define SINTEF_GRIDCPLUSPLUS_HEADER


#include "grid.h"
#include <stdexcept>

class GridCplusplus
{
public:
    /// @brief
    /// Initialize the grid.
    /// @tparam Grid This must conform to the SimpleGrid concept.
    /// @param grid The grid object.
    template <class Grid>
    void init(const Grid& grid)
    {
        buildTopology(grid);
        buildGeometry(grid);
    }

    /// Access the underlying C grid.
    grid_t* c_grid()
    {
        return &g_;
    }

private:
    grid_t g_;
    // Topology storage.
    std::vector<int> face_nodes_;
    std::vector<int> face_nodepos_;
    std::vector<int> face_cells_;
    std::vector<int> cell_faces_;
    std::vector<int> cell_facepos_;
    // Geometry storage.
    std::vector<double> node_coordinates_;
    std::vector<double> face_centroids_;
    std::vector<double> face_areas_;
    std::vector<double> face_normals_;
    std::vector<double> cell_centroids_;
    std::vector<double> cell_volumes_;

    /// Build (copy of) topological structure from grid.
    template <class Grid>
    void buildTopology(const Grid& grid)
    {
        // Face topology.
        int num_cells = grid.numCells();
        int num_nodes = grid.numVertices();
        int num_faces = grid.numFaces();
        face_nodepos_.resize(num_faces + 1);
        int facenodecount = 0;
        for (int f = 0; f < num_faces; ++f) {
            face_nodepos_[f] = facenodecount;
            facenodecount += grid.numFaceVertices(f);
        }
        face_nodepos_.back() = facenodecount;
        face_nodes_.resize(facenodecount);
        for (int f = 0; f < num_faces; ++f) {
            for (int local = 0; local < grid.numFaceVertices(f); ++local) {
                face_nodes_[face_nodepos_[f] + local] = grid.faceVertex(f, local);
            }
        }
        face_cells_.resize(2*num_faces);
        for (int f = 0; f < num_faces; ++f) {
            face_cells_[2*f] = grid.faceCell(f, 0);
            face_cells_[2*f + 1] = grid.faceCell(f, 1);
        }

        // Cell topology.
        int cellfacecount = 0;
        cell_facepos_.resize(num_cells + 1);
        for (int c = 0; c < num_cells; ++c) {
            cell_facepos_[c] = cellfacecount;
            cellfacecount += grid.numCellFaces(c);
        }
        cell_facepos_.back() = cellfacecount;
        cell_faces_.resize(cellfacecount);
        for (int c = 0; c < num_cells; ++c) {
            for (int local = 0; local < grid.numCellFaces(c); ++local) {
                cell_faces_[cell_facepos_[c] + local] = grid.cellFace(c, local);
            }
        }

        // Set C grid members.
        g_.dimensions = Grid::dimension;
        g_.number_of_cells = grid.numCells();
        g_.number_of_faces = grid.numFaces();
        g_.number_of_nodes = grid.numVertices();
        g_.face_nodes = &face_nodes_[0];
        g_.face_nodepos = &face_nodepos_[0];
        g_.face_cells = &face_cells_[0];
        g_.cell_faces = &cell_faces_[0];
        g_.cell_facepos = &cell_facepos_[0];
    }


    /// Build (copy of) geometric properties of grid.
    /// Assumes that buildTopology() has been called.
    template <class Grid>
    void buildGeometry(const Grid& grid)
    {
        // Node geometry.
        int num_cells = grid.numCells();
        int num_nodes = grid.numVertices();
        int num_faces = grid.numFaces();
        int dim = Grid::dimension;
        node_coordinates_.resize(dim*num_nodes);
        for (int n = 0; n < num_nodes; ++n) {
            for (int dd = 0; dd < dim; ++dd) {
                node_coordinates_[dim*n + dd] = grid.vertexPosition(n)[dd];
            }
        }

        // Face geometry.
        face_centroids_.resize(dim*num_faces);
        face_areas_.resize(num_faces);
        face_normals_.resize(dim*num_faces);
        for (int f = 0; f < num_faces; ++f) {
            face_areas_[f] = grid.faceArea(f);
            for (int dd = 0; dd < dim; ++dd) {
                face_centroids_[dim*f + dd] = grid.faceCentroid(f)[dd];
                face_normals_[dim*f + dd] = grid.faceNormal(f)[dd];
            }
        }

        // Cell geometry.
        cell_centroids_.resize(dim*num_cells);
        cell_volumes_.resize(num_cells);
        for (int c = 0; c < num_cells; ++c) {
            cell_volumes_[c] = grid.cellVolume(c);
            for (int dd = 0; dd < dim; ++dd) {
                cell_centroids_[dim*c + dd] = grid.cellCentroid(c)[dd];
            }
        }

        // Set C grid members.
        g_.node_coordinates = &node_coordinates_[0];
        g_.face_centroids = &face_centroids_[0];
        g_.face_areas = &face_areas_[0];
        g_.face_normals = &face_normals_[0];
        g_.cell_centroids = &cell_centroids_[0];
        g_.cell_volumes = &cell_volumes_[0];
    }
};


#endif // SINTEF_GRIDCPLUSPLUS_HEADER
