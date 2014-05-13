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

#ifndef OPM_BLACKOILPVTPROPERTIES_HEADER_INCLUDED
#define OPM_BLACKOILPVTPROPERTIES_HEADER_INCLUDED

#include <opm/core/props/pvt/PvtInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Utility/PvtoTable.hpp>

#include <string>
#include <memory>

namespace Opm
{

    /// Class collecting the pvt properties for all active phases.
    /// For all the methods, the following apply:
    /// - p and z are expected to be of size n and n*num_phases, respectively.
    /// - pvtTableIdx specifies the PVT table to be used for each data
    ///               point and is thus expected to be an array of size n
    /// - Output arrays shall be of size n*num_phases, and must be valid
    ///   before calling the method.
    /// NOTE: The difference between this interface and the one defined
    /// by PvtInterface is that this collects all phases' properties,
    /// and therefore the output arrays are of size n*num_phases as opposed
    /// to size n in PvtInterface.
    class BlackoilPvtProperties : public BlackoilPhases
    {
    public:
        /// Default constructor.
        BlackoilPvtProperties();

        /// Initialize from deck.
        ///
        /// \param deck     An input deck from the opm-parser module.
        void init(Opm::DeckConstPtr deck,
                  int samples);

        /// \return   Object describing the active phases.
        PhaseUsage phaseUsage() const;

        /// Number of active phases.
        int numPhases() const;

        /// For each canonical phase, indicates if it is
        /// active or not (boolean usage of int).
        /// \return  Array of size MaxNumPhases
        const int* phaseUsed() const;

        /// Positions of canonical phases in arrays of phase
        /// properties (saturations etc.).
        /// \return  Array of size MaxNumPhases
        const int* phasePosition() const;

        /// Densities of stock components at surface conditions.
        /// \return  Array of size numPhases().
        const double* surfaceDensities(int regionIdx = 0) const;

        /// Viscosity as a function of p and z.
        void mu(const int n,
                const int *pvtTableIdx,
                const double* p,
                const double* z,
                double* output_mu) const;

        /// Formation volume factor as a function of p and z.
        void B(const int n,
               const int *pvtTableIdx,
               const double* p,
               const double* z,
               double* output_B) const;

        /// Formation volume factor and p-derivative as functions of p and z.
        void dBdp(const int n,
                  const int *pvtTableIdx,
                  const double* p,
                  const double* z,
                  double* output_B,
                  double* output_dBdp) const;

        /// Solution factor as a function of p and z.
        void R(const int n,
               const int *pvtTableIdx,
               const double* p,
               const double* z,
               double* output_R) const;

        /// Solution factor and p-derivative as functions of p and z.
        void dRdp(const int n,
                  const int *pvtTableIdx,
                  const double* p,
                  const double* z,
                  double* output_R,
                  double* output_dRdp) const;

    private:
        // Disabling copying (just to avoid surprises, since we use shared_ptr).
        BlackoilPvtProperties(const BlackoilPvtProperties&);
        BlackoilPvtProperties& operator=(const BlackoilPvtProperties&);

        PhaseUsage phase_usage_;

        // The PVT properties. We need to store one object per PVT
        // region per active fluid phase.
        std::vector<std::shared_ptr<PvtInterface> > props_;
        std::vector<std::array<double, MaxNumPhases> > densities_;

        mutable std::vector<double> data1_;
        mutable std::vector<double> data2_;
    };

}



#endif // OPM_BLACKOILPVTPROPERTIES_HEADER_INCLUDED
