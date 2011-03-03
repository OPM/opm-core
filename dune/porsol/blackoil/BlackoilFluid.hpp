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
            const std::vector<double>& dens = parser.getDENSITY().densities_[0];
            surface_densities_[Oil] = dens[0];
            surface_densities_[Water] = dens[1];
            surface_densities_[Gas] = dens[2];
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
        const CompVec& surfaceDensities() const
        {
            return surface_densities_;
        }
    private:
        FluidMatrixInteractionBlackoilParams<double> fmi_params_;
        CompVec surface_densities_;
    };




    /// Container for all fluid data needed by solvers.
    struct BlackoilFluidData : public BlackoilDefs
    {
        // Per-cell data.
        std::vector<double> totcompr;     // Total compressibility.
        std::vector<double> totphasevol;  // Total volume filled by fluid phases.
        std::vector<double> voldiscr;     // Volume discrepancy = (totphasevol - porevol)/dt
        std::vector<double> relvoldiscr;  // Relative volume discrepancy = (totphasevol - porevol)/porevol
        std::vector<double> cellA;        // A = RB^{-1}. Fortran ordering, flat storage.
        std::vector<PhaseVec> saturation; // Saturation.
        std::vector<PhaseVec> frac_flow;  // Fractional flow.
        std::vector<PhaseVec> rel_perm;   // Relative permeability.
        std::vector<PhaseVec> viscosity;  // Viscosity.
        // Per-face data.
        std::vector<double> faceA;        // A = RB^{-1}. Fortran ordering, flat storage.
        std::vector<double> phasemobf;    // Phase mobilities. Flat storage (numPhases per face).
        std::vector<PhaseVec> phasemobc;  // Phase mobilities per cell.

    public:
        template <class Grid, class Rock>
        void compute(const Grid& grid,
                     const Rock& rock,
                     const BlackoilFluid& fluid,
                     const std::vector<PhaseVec>& cell_pressure,
                     const std::vector<PhaseVec>& face_pressure,
                     const std::vector<CompVec>& cell_z,
                     const CompVec& bdy_z,
                     const double dt)
        {
            int num_cells = cell_z.size();
            ASSERT(num_cells == grid.numCells());
            int num_faces = face_pressure.size();
            ASSERT(num_faces == grid.numFaces());
            const int np = numPhases;
            const int nc = numComponents;
            BOOST_STATIC_ASSERT(np == nc);
            totcompr.resize(num_cells);
            totphasevol.resize(num_cells);
            voldiscr.resize(num_cells);
            relvoldiscr.resize(num_cells);
            saturation.resize(num_cells);
            frac_flow.resize(num_cells);
            rel_perm.resize(num_cells);
            viscosity.resize(num_cells);
            cellA.resize(num_cells*nc*np);
            faceA.resize(num_faces*nc*np);
            phasemobf.resize(np*num_faces);
            phasemobc.resize(num_cells);
            BOOST_STATIC_ASSERT(np == 3);
#pragma omp parallel for
            for (int cell = 0; cell < num_cells; ++cell) {
                FluidStateBlackoil state = fluid.computeState(cell_pressure[cell], cell_z[cell]);
                totcompr[cell] = state.total_compressibility_;
                totphasevol[cell] = state.total_phase_volume_;
                double pv = rock.porosity(cell)*grid.cellVolume(cell);
                voldiscr[cell] = (totphasevol[cell] - pv)/dt;
                relvoldiscr[cell] = std::fabs(totphasevol[cell] - pv)/pv;
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
#pragma omp parallel for
            for (int face = 0; face < num_faces; ++face) {
                int c[2] = { grid.faceCell(face, 0), grid.faceCell(face, 1) };
                PhaseVec phase_p[2];
                PhaseVec phase_mob[2];
                CompVec face_z(0.0);
                bool bdy = false;
                bool inflow_bdy = false;
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        phase_p[j] = cell_pressure[c[j]];
                        phase_mob[j] = phasemobc[c[j]];
                        face_z += cell_z[c[j]];
                    } else {
                        bdy = true;
                        phase_p[j] = face_pressure[face];
                        /// \TODO with capillary pressures etc., what is an inflow bdy.
                        /// Using Liquid phase pressure here.
                        inflow_bdy = face_pressure[face][Liquid]
                            > cell_pressure[c[(j+1)%2]][Liquid];
                        if (inflow_bdy) {
                            FluidStateBlackoil bdy_state = fluid.computeState(face_pressure[face], bdy_z);
                            phase_mob[j] = bdy_state.mobility_;
                            face_z += bdy_z;
                        } else {
                            phase_p[j] = -1e100; // To ensure correct upwinding.
                            // No need to set phase_mob[j].
                        }
                    }
                }
                if (!bdy || inflow_bdy) {
                    face_z *= 0.5;
                }
                for (int phase = 0; phase < np; ++phase) {
                    if (phase_p[0][phase] == phase_p[1][phase]) {
                        // Average mobilities.
                        double aver = 0.5*(phase_mob[0][phase] + phase_mob[1][phase]);
                        phasemobf[np*face + phase] = aver;
                    } else {
                        // Upwind mobilities.
                        int upwind = (phase_p[0][phase] > phase_p[1][phase]) ? 0 : 1;
                        phasemobf[np*face + phase] = phase_mob[upwind][phase];
                    }
                }
                FluidStateBlackoil face_state = fluid.computeState(face_pressure[face], face_z);
                std::copy(face_state.phase_to_comp_, face_state.phase_to_comp_ + nc*np, &faceA[face*nc*np]);
            }

        }
    };




} // namespace Opm

#endif // OPM_BLACKOILFLUID_HEADER_INCLUDED
