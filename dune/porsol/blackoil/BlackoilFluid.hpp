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

#ifndef OPM_BLACKOILFLUID_HEADER_INCLUDED
#define OPM_BLACKOILFLUID_HEADER_INCLUDED


#include <dune/porsol/blackoil/fluid/FluidMatrixInteractionBlackoil.hpp>
#include <dune/porsol/blackoil/fluid/FluidSystemBlackoil.hpp>
#include <dune/porsol/blackoil/fluid/FluidStateBlackoil.hpp>
#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/fvector.hh>
#include <vector>


namespace Opm
{


    /// Forward declaration for typedef i BlackoilFluid.
    class BlackoilFluidData;


    /// Class responsible for computing all fluid properties from
    /// face pressures and composition.
    class BlackoilFluid : public BlackoilDefs
    {
    public:
        typedef FluidStateBlackoil FluidState;
        typedef BlackoilFluidData FluidData;

        void init(const Dune::EclipseGridParser& parser)
        {
            fmi_params_.init(parser);
            FluidSystemBlackoil<>::init(parser);
        }
        FluidState computeState(PhaseVec phase_pressure, CompVec z) const
        {
            FluidState state;
            state.temperature_ = 300;
            state.phase_pressure_ = phase_pressure;
            state.surface_volume_ = z;
            FluidSystemBlackoil<>::computeEquilibrium(state); // Sets everything but relperm and mobility.
            FluidMatrixInteractionBlackoil<double>::kr(state.relperm_, fmi_params_, state.saturation_, state.temperature_);
            for (int phase = 0; phase < numPhases; ++phase) {
                state.mobility_[phase] = state.relperm_[phase]/state.viscosity_[phase];
            }
            return state;
        }
    private:
        FluidMatrixInteractionBlackoilParams<double> fmi_params_;
    };




    /// Container for all fluid data needed by solvers.
    struct BlackoilFluidData : public BlackoilDefs
    {
        // Per-cell data.
        std::vector<double> totcompr;     // Total compressibility.
        std::vector<double> totphasevol;  // Total volume filled by fluid phases.
        std::vector<double> cellA;        // A = RB^{-1}. Fortran ordering, flat storage.
        std::vector<PhaseVec> saturation; // Saturation.
        std::vector<PhaseVec> frac_flow;  // Fractional flow.
        std::vector<PhaseVec> rel_perm;   // Relative permeability.
        std::vector<PhaseVec> viscosity;  // Viscosity.
        // Per-face data.
        std::vector<double> faceA;        // A = RB^{-1}. Fortran ordering, flat storage.
        std::vector<double> phasemobf;    // Phase mobilities. Flat storage (numPhases per face).
    private:
        std::vector<PhaseVec> phasemobc;  // Just a helper. Mobilities per cell.

    public:
        template <class Grid>
        void compute(const Grid& grid,
                     const BlackoilFluid& fluid,
                     const std::vector<PhaseVec>& phase_pressure,
                     const std::vector<PhaseVec>& phase_pressure_face,
                     const std::vector<CompVec>& z,
                     const CompVec& bdy_z)
        {
            int num_cells = z.size();
            ASSERT(num_cells == grid.numCells());
            int num_faces = phase_pressure_face.size();
            ASSERT(num_faces == grid.numFaces());
            const int np = numPhases;
            const int nc = numComponents;
            BOOST_STATIC_ASSERT(np == nc);
            totcompr.resize(num_cells);
            totphasevol.resize(num_cells);
            saturation.resize(num_cells);
            frac_flow.resize(num_cells);
            rel_perm.resize(num_cells);
            viscosity.resize(num_cells);
            cellA.resize(num_cells*nc*np);
            faceA.resize(num_faces*nc*np);
            phasemobf.resize(np*num_faces);
            phasemobc.resize(num_cells);
            PhaseVec mob;
            BOOST_STATIC_ASSERT(np == 3);
            for (int cell = 0; cell < num_cells; ++cell) {
                FluidStateBlackoil state = fluid.computeState(phase_pressure[cell], z[cell]);
                totcompr[cell] = state.total_compressibility_;
                totphasevol[cell] = state.total_phase_volume_;
                saturation[cell] = state.saturation_;
                rel_perm[cell] = state.relperm_;
                viscosity[cell] = state.viscosity_;
                phasemobc[cell] = state.mobility_;
                std::copy(state.phase_to_comp_, state.phase_to_comp_ + nc*np, &cellA[cell*nc*np]);
                // Fractional flow must be calculated.
                double total_mobility = 0.0;
                for (int phase = 0; phase < numPhases; ++phase) {
                    total_mobility += state.mobility_[phase];
                }
                frac_flow[cell] = state.mobility_;
                frac_flow[cell] /= total_mobility;
            }
            // Set phasemobf to average of cells' phase mobs, if pressures are equal, else use upwinding.
            // Set faceA by using average of cells' z and face pressures.
            for (int face = 0; face < num_faces; ++face) {
                int c[2] = { grid.faceCell(face, 0), grid.faceCell(face, 1) };
                PhaseVec phase_p[2];
                CompVec z_face(0.0);
                int num = 0;
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        phase_p[j] = phase_pressure[c[j]];
                        z_face += z[c[j]];
                        ++num;
                    } else {
                        // Boundaries get essentially -inf pressure for upwinding purpose. \TODO handle BCs.
                        phase_p[j] = PhaseVec(-1e100);
                        // \TODO The two lines below are wrong for outflow faces.
                        z_face += bdy_z;
                        ++num;
                    }
                }
                z_face /= double(num);
                for (int phase = 0; phase < np; ++phase) {
                    if (phase_p[0][phase] == phase_p[1][phase]) {
                        // Average mobilities.
                        double aver = 0.5*(phasemobc[c[0]][phase] + phasemobc[c[1]][phase]);
                        phasemobf[np*face + phase] = aver;
                    } else {
                        // Upwind mobilities.
                        int upwind = (phase_p[0][phase] > phase_p[1][phase]) ? 0 : 1;
                        phasemobf[np*face + phase] = phasemobc[c[upwind]][phase];
                    }
                }
                FluidStateBlackoil face_state = fluid.computeState(phase_pressure_face[face], z_face);
                std::copy(face_state.phase_to_comp_, face_state.phase_to_comp_ + nc*np, &faceA[face*nc*np]);
            }

        }
    };




} // namespace Opm

#endif // OPM_BLACKOILFLUID_HEADER_INCLUDED
