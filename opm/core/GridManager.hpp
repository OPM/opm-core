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

#ifndef OPM_GRIDMANAGER_HEADER_INCLUDED
#define OPM_GRIDMANAGER_HEADER_INCLUDED


struct UnstructuredGrid;


namespace Opm
{

    class EclipseGridParser;

    /// This class manages an Opm::UnstructuredGrid in the sense that it
    /// encapsulates creation and destruction of the grid.
    /// The following grid types can be constructed:
    ///   - 3d corner-point grids (from deck input)
    ///   - 3d tensor grids (from deck input)
    ///   - 2d cartesian grids
    ///   - 3d cartesian grids
    /// The resulting UnstructuredGrid is available through the c_grid() method.
    class GridManager
    {
    public:
        /// Construct a 3d corner-point grid or tensor grid from a deck.
        GridManager(const Opm::EclipseGridParser& deck);

        /// Construct a 2d cartesian grid with cells of unit size.
        GridManager(int nx, int ny);

        /// Construct a 3d cartesian grid with cells of unit size.
        GridManager(int nx, int ny, int nz);

        /// Construct a 3d cartesian grid with cells of size [dx, dy, dz].
        GridManager(int nx, int ny, int nz,
                    double dx, double dy, double dz);

        /// Destructor.
        ~GridManager();

        /// Access the managed UnstructuredGrid.
        /// The method is named similarly to c_str() in std::string,
        /// to make it clear that we are returning a C-compatible struct.
        const UnstructuredGrid* c_grid() const;

    private:
        // Disable copying and assignment.
        GridManager(const GridManager& other);
        GridManager& operator=(const GridManager& other);

        // Construct corner-point grid from deck.
        void initFromDeckCornerpoint(const Opm::EclipseGridParser& deck);
        // Construct tensor grid from deck.
        void initFromDeckTensorgrid(const Opm::EclipseGridParser& deck);

        // The managed UnstructuredGrid.
        UnstructuredGrid* ug_;
    };

} // namespace Opm

#endif // OPM_GRIDMANAGER_HEADER_INCLUDED
