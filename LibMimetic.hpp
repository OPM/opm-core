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
#include "GridCplusplus.hpp"
#include <stdexcept>


/// @brief
/// Encapsulates the ifsh (= incompressible flow solver hybrid) solver modules.
class Ifsh
{
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
        if (data_) {
            ifsh_destroy(data_);
        }
    }

    /// @brief
    /// Initialize the solver's structures for a given grid (at some point also well pattern).
    /// @tparam Grid This must conform to the SimpleGrid concept.
    /// @param grid The grid object.
    template <class Grid>
    void init(const Grid& grid)
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
    }

    /// @brief
    /// Assemble the sparse system.
    void assemble(const std::vector<double>& sources)
    {
        flowbc_t* bc = 0; // TODO
        double* src = 0; // &sources[0]; // TODO
        double* Binv = 0; // TODO
        double* gpress = 0; // TODO

        // All well related things are zero.
        well_control_t* wctrl = 0;
        double* WI = 0;
        double* wdp = 0;

        double* totmob = 0; // TODO
        double* omega = 0; // TODO
        ifsh_assemble(bc, src, Binv, gpress, wctrl, WI, wdp, totmob, omega, data_);
    }

private:
    // Disabling copy and assigment for now.
    Ifsh(const Ifsh&);
    Ifsh& operator=(const Ifsh&);

    // Solver data.
    ifsh_data* data_;
    // Grid.
    GridCplusplus grid_;
};

#if 0
void
ifsh_assemble(flowbc_t         *bc,
              double           *src,
              double           *Binv,
              double           *gpress,
              well_control_t   *wctrl,
              double           *WI,
              double           *wdp,
              double           *totmob, /* \sum_i \lambda_i */
              double           *omega,  /* \sum_i \rho_i f_i */
              struct ifsh_data *h);

void
ifsh_press_flux(grid_t *G, struct ifsh_data *h, double *src,
                double *cpress, double *fflux);
#endif



#endif // SINTEF_LIBMIMETIC_HEADER
