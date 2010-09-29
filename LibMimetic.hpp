//===========================================================================
//
// File: LibMimetic.hpp
//
// Created: Thu Sep 23 20:00:49 2010
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

  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the code.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SINTEF_LIBMIMETIC_HEADER
#define SINTEF_LIBMIMETIC_HEADER

#include "ifsh.h"
#include "mimetic.h"
#include "GridCplusplus.hpp"
#include <stdexcept>



/// @brief
/// Encapsulates the ifsh (= incompressible flow solver hybrid) solver modules.
class Ifsh
{
public:
    /// @brief
    /// Default constructor, does nothing.
    Ifsh()
        : data_(0)
    {
    }

    /// @brief
    /// Destructor.
    ~Ifsh()
    {
        ifsh_destroy(data_);
    }

    /// @brief
    /// Initialize the solver's structures for a given grid (at some point also well pattern).
    /// @tparam Grid This must conform to the SimpleGrid concept.
    /// @param grid The grid object.
    template <class Grid>
    void init(const Grid& grid, const double* perm)
    {
        // Build C grid structure.
        grid_.init(grid);

        // Build (empty for now) C well structure.
        well_t* w = 0;

        // Initialize ifsh data.
        data_ = ifsh_construct(grid_.c_grid(), w);
        if (!data_) {
            throw std::runtime_error("Failed to initialize ifsh solver.");
        }

        // Compute inner products.
        int num_cells = grid.numCells();
        ncf_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            ncf_[cell] = grid.numCellFaces(cell);
        }
        // Zero gravity (for now)
        int ngconn  = grid_.c_grid()->cell_facepos[num_cells];
        int ngconn2 = data_->sum_ngconn2;
        Binv_.resize(ngconn2);
        gpress_.clear();
        gpress_.resize(ngconn, 0.0); // Zero gravity for now!

        grid_t* g = grid_.c_grid();
        mim_ip_simple_all(g->number_of_cells, g->dimensions,
                          data_->max_ngconn,
                          &ncf_[0], g->cell_facepos, g->cell_faces,
                          g->face_cells, g->face_centroids,
                          g->face_normals, g->face_areas,
                          g->cell_centroids, g->cell_volumes,
                          const_cast<double*>(perm), &Binv_[0]);
    }

    /// @brief
    /// Assemble the sparse system.
    void assemble(const std::vector<double>& sources,
                  const std::vector<double>& total_mobilities,
                  const std::vector<double>& omegas)
    {
        // Noflow conditions for now.
        int num_faces = grid_.c_grid()->number_of_faces;
        std::vector<flowbc_type> bc_types(num_faces, FLUX);
        std::vector<double> bc_vals(num_faces, 0);
        flowbc_t bc = { &bc_types[0], &bc_vals[0] };

        // Source terms from user.
        double* src = const_cast<double*>(&sources[0]); // Ugly? Yes. Safe? I think so.

        // Inner products are precomputed.
        double* Binv = &Binv_[0];

        // Gravity contribs are precomputed.
        double* gpress = &gpress_[0];

        // All well related things are zero.
        well_control_t* wctrl = 0;
        double* WI = 0;
        double* wdp = 0;

        double* totmob = const_cast<double*>(&total_mobilities[0]);
        double* omega = const_cast<double*>(&omega[0]);

        // Zero the linalg structures.
        csrmatrix_zero(data_->A);
        for (std::size_t i = 0; i < data_->A->m; i++) {
            data_->b[i] = 0.0;
        }

        ifsh_assemble(&bc, src, Binv, gpress, wctrl, WI, wdp, totmob, omega, data_);
    }

    /// Encapsulate a sparse linear system in CSR format.
    struct LinearSystem
    {
        int n;
        int nnz;
        int* ia;
        int* ja;
        double* sa;
        double* b;
        double* x;
    };

    /// Access the linear system assembled.
    void linearSystem(LinearSystem& s)
    {
        s.n = data_->A->n;
        s.nnz = data_->A->nnz;
        s.ia = data_->A->ia;
        s.ja = data_->A->ja;
        s.sa = data_->A->sa;
        s.b = data_->b;
        s.x = data_->x;
    }

    /// Compute cell pressures and face fluxes.
    void computePressuresAndFluxes(std::vector<double>& cell_pressures,
                                   std::vector<double>& face_fluxes)
    {
        int num_cells = grid_.c_grid()->number_of_cells;
        int num_faces = grid_.c_grid()->number_of_faces;
        cell_pressures.clear();
        cell_pressures.resize(num_cells, 0.0);
        face_fluxes.clear();
        face_fluxes.resize(num_faces, 0.0);
        ifsh_press_flux(grid_.c_grid(), data_, &cell_pressures[0], &face_fluxes[0], 0, 0);
    }

    /// Compute cell fluxes from face fluxes.
    void faceFluxToCellFlux(const std::vector<double>& face_fluxes,
                            std::vector<double>& cell_fluxes)
    {
        const grid_t& g = *(grid_.c_grid());
        int num_cells = g.number_of_cells;
        cell_fluxes.resize(g.cell_facepos[num_cells]);
        for (int cell = 0; cell < num_cells; ++cell) {
            for (int hface = g.cell_facepos[cell]; hface < g.cell_facepos[cell + 1]; ++hface) {
                int face = g.cell_faces[hface];
                bool pos = (g.face_cells[2*face] == cell);
                cell_fluxes[hface] = pos ? face_fluxes[face] : -face_fluxes[face];
            }
        }
    }

    /// Access the number of connections (faces) per cell.
    const std::vector<int>& numCellFaces()
    {
        return ncf_;
    }

private:
    // Disabling copy and assigment for now.
    Ifsh(const Ifsh&);
    Ifsh& operator=(const Ifsh&);

    // Solver data.
    ifsh_data* data_;
    // Grid.
    GridCplusplus grid_;
    // Number of faces per cell.
    std::vector<int> ncf_;
    // B^{-1} storage.
    std::vector<double> Binv_;
    // Gravity contributions.
    std::vector<double> gpress_;
    // Total mobilities.
    std::vector<double> totmob_;
    // Gravity coefficients (\omega = sum_{i = 1}^{num phases}f_i \rho_i[TODO: check this]).
    std::vector<double> omega_;
};




#endif // SINTEF_LIBMIMETIC_HEADER
