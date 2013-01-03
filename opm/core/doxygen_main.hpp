/*
  Copyright 2012, 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_DOXYGEN_MAIN_HEADER_INCLUDED
#define OPM_DOXYGEN_MAIN_HEADER_INCLUDED


/** \mainpage Documentation for the opm-core library.

The following are the main library features:

<h3>Grid handling</h3>

The library defines a simple grid interface through the struct
UnstructuredGrid. This can be used both from C and C++ code, and has a
structure that is similar to MRST grids.

Cells can have arbitrary shapes, and arbitrary numbers of neighbours.
This flexibility allows us to handle complex grids such as faulted
corner-point grids or PEBI grids. The structure is suited for
computation with (for example) finite volume methods or mimetic finite
differences. It is less ideal for finite element computations, as it
does not provide any reference element mappings.

There are multiple construction functions for UnstructuredGrid, for
example create_grid_cart2d(), create_grid_hexa3d(), read_grid() and
create_grid_cornerpoint(). The function destroy_grid() frees the
resources used by a grid.

For C++ users, the Opm::GridManager class can be used to encapsulate
creation and destruction of an UnstructuredGrid. The method
Opm::GridManager::c_grid() provides access to the underlying
UnstructuredGrid struct.


<h3>Well handling</h3>

<h3>Pressure solvers</h3>

<h3>Transport solvers</h3>

<h3>Various utilities</h3>

*/

#endif // OPM_DOXYGEN_MAIN_HEADER_INCLUDED
