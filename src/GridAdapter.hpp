/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_GRIDADAPTER_HEADER_INCLUDED
#define OPM_GRIDADAPTER_HEADER_INCLUDED


#include "grid.h"
#include <stdexcept>

class GridAdapter
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

    grid_t* c_grid()
    {
       return &g_;
    }
    /// Access the underlying C grid.
    const grid_t* c_grid() const
    {
        return &g_;
    }

    // ------ Forwarding the same interface that init() expects ------
    //
    // This is only done in order to verify that init() works correctly.
    //

    enum { dimension = 3 }; // This is actually a hack used for testing (dim is a runtime parameter).

    struct Vector
    {
        explicit Vector(const double* source)
        {
            for (int i = 0; i < dimension; ++i) {
                data[i] = source[i];
            }
        }
        double& operator[] (const int ix)
        {
            return data[ix];
        }
        double operator[] (const int ix) const
        {
            return data[ix];
        }
        double data[dimension];
    };

    // Topology
    int numCells() const
    {
        return g_.number_of_cells;
    }
    int numFaces() const
    {
        return g_.number_of_faces;
    }
    int numVertices() const
    {
        return g_.number_of_nodes;
    }

    int numCellFaces(int cell) const
    {
        return cell_facepos_[cell + 1] - cell_facepos_[cell];
    }
    int cellFace(int cell, int local_index) const
    {
        return cell_faces_[cell_facepos_[cell] + local_index];
    }
    int faceCell(int face, int local_index) const
    {
        return face_cells_[2*face + local_index];
    }
    int numFaceVertices(int face) const
    {
        return face_nodepos_[face + 1] - face_nodepos_[face];
    }
    int faceVertex(int face, int local_index) const
    {
        return face_nodes_[face_nodepos_[face] + local_index];
    }

    // Geometry
    Vector vertexPosition(int vertex) const
    {
        return Vector(&node_coordinates_[g_.dimensions*vertex]);
    }
    double faceArea(int face) const
    {
        return face_areas_[face];
    }
    Vector faceCentroid(int face) const
    {
        return Vector(&face_centroids_[g_.dimensions*face]);
    }
    Vector faceNormal(int face) const
    {
        Vector fn(&face_normals_[g_.dimensions*face]);
        // We must renormalize since the stored normals are
        // 'unit normal * face area'.
        double invfa = 1.0 / faceArea(face);
        for (int i = 0; i < dimension; ++i) {
            fn[i] *= invfa;
        }
        return fn;
    }
    double cellVolume(int cell) const
    {
        return cell_volumes_[cell];
    }
    Vector cellCentroid(int cell) const
    {
        return Vector(&cell_centroids_[g_.dimensions*cell]);
    }

    bool operator==(const GridAdapter& other)
    {
        return face_nodes_ == other.face_nodes_
            && face_nodepos_ == other.face_nodepos_
            && face_cells_ == other.face_cells_
            && cell_faces_ == other.cell_faces_
            && cell_facepos_ == other.cell_facepos_
            && node_coordinates_ == other.node_coordinates_
            && face_centroids_ == other.face_centroids_
            && face_areas_ == other.face_areas_
            && face_normals_ == other.face_normals_
            && cell_centroids_ == other.cell_centroids_
            && cell_volumes_ == other.cell_volumes_;
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
                face_normals_[dim*f + dd] = grid.faceNormal(f)[dd]*face_areas_[f];
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


#endif // OPM_GRIDADAPTER_HEADER_INCLUDED
