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
            FluidMatrixInteractionBlackoil<double>::dkr(state.drelperm_, fmi_params_, state.saturation_, state.temperature_);
            for (int phase = 0; phase < numPhases; ++phase) {
                state.mobility_[phase] = state.relperm_[phase]/state.viscosity_[phase];
                for (int p2 = 0; p2 < numPhases; ++p2) {
                    state.dmobility_[phase][p2] = state.drelperm_[phase][p2]/state.viscosity_[phase];
                }
            }
            return state;
        }
        const CompVec& surfaceDensities() const
        {
            return surface_densities_;
        }

        /// \param[in] A state matrix in fortran ordering
        PhaseVec phaseDensities(const double* A) const
        {
            PhaseVec phase_dens(0.0);
            for (int phase = 0; phase < numPhases; ++phase) {
                for (int comp = 0; comp < numComponents; ++comp) {
                    phase_dens[phase] += A[numComponents*phase + comp]*surface_densities_[comp];
                }
            }
            return phase_dens;
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
        std::vector<double> totphasevol_density;  // Total volume filled by fluid phases
                                                  // per pore volume.
        std::vector<double> voldiscr;     // Volume discrepancy = (totphasevol_dens - 1)*pv/dt
        std::vector<double> relvoldiscr;  // Relative volume discrepancy = |totphasevol_dens - 1|
        std::vector<double> cellA;        // A = RB^{-1}. Fortran ordering, flat storage.
        std::vector<PhaseVec> saturation; // Saturation.
        std::vector<PhaseVec> frac_flow;  // Fractional flow.
        std::vector<PhaseVec> rel_perm;   // Relative permeability.
        std::vector<PhaseVec> viscosity;  // Viscosity.

        // Extra data.
        std::vector<double> expjacterm;

        // Per-face data.
        std::vector<double> faceA;           // A = RB^{-1}. Fortran ordering, flat storage.
        std::vector<double> phasemobf;       // Phase mobilities. Flat storage (numPhases per face).
        std::vector<PhaseVec> phasemobc;     // Phase mobilities per cell.
        std::vector<double> phasemobf_deriv; // Phase mobility derivatives. Flat storage (numPhases^2 per face).
        typedef Dune::FieldVector<Dune::FieldVector<Scalar, numPhases>, numPhases> PhaseMat;
        std::vector<PhaseMat> phasemobc_deriv; // Phase mobilities derivatives per cell.
        std::vector<double> gravcapf;      // Gravity (\rho g \delta z-ish) contribution per face, flat storage.

    public:
        template <class Grid, class Rock>
        void compute(const Grid& grid,
                     const Rock& rock,
                     const BlackoilFluid& fluid,
                     const typename Grid::Vector gravity,
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
            bool nonzero_gravity = gravity.two_norm() > 0.0;
            BOOST_STATIC_ASSERT(np == nc);
            totcompr.resize(num_cells);
            totphasevol_density.resize(num_cells);
            voldiscr.resize(num_cells);
            relvoldiscr.resize(num_cells);
            saturation.resize(num_cells);
            frac_flow.resize(num_cells);
            rel_perm.resize(num_cells);
            viscosity.resize(num_cells);
            expjacterm.resize(num_cells);
            cellA.resize(num_cells*nc*np);
            faceA.resize(num_faces*nc*np);
            phasemobf.resize(np*num_faces);
            phasemobc.resize(num_cells);
            phasemobf_deriv.resize(np*np*num_faces);
            phasemobc_deriv.resize(np*np*num_cells);
            gravcapf.resize(np*num_faces);
            BOOST_STATIC_ASSERT(np == 3);
#pragma omp parallel for
            for (int cell = 0; cell < num_cells; ++cell) {
                FluidStateBlackoil state = fluid.computeState(cell_pressure[cell], cell_z[cell]);
                totcompr[cell] = state.total_compressibility_;
                totphasevol_density[cell] = state.total_phase_volume_density_;
                double pv = rock.porosity(cell)*grid.cellVolume(cell);
                voldiscr[cell] = (totphasevol_density[cell] - 1.0)*pv/dt;
                relvoldiscr[cell] = std::fabs(totphasevol_density[cell] - 1.0);
                saturation[cell] = state.saturation_;
                rel_perm[cell] = state.relperm_;
                viscosity[cell] = state.viscosity_;
                phasemobc[cell] = state.mobility_;
                phasemobc_deriv[cell] = state.dmobility_;
                std::copy(state.phase_to_comp_, state.phase_to_comp_ + nc*np, &cellA[cell*nc*np]);
                // Fractional flow must be calculated.
                double total_mobility = 0.0;
                for (int phase = 0; phase < numPhases; ++phase) {
                    total_mobility += state.mobility_[phase];
                }
                frac_flow[cell] = state.mobility_;
                frac_flow[cell] /= total_mobility;

                // Experimental
                expjacterm[cell] = state.experimental_term_;
            }

            // Obtain properties from both sides of the face.
#pragma omp parallel for
            for (int face = 0; face < num_faces; ++face) {
                typedef typename Grid::Vector Vec;
                Vec fc = grid.faceCentroid(face);
                int c[2] = { grid.faceCell(face, 0), grid.faceCell(face, 1) };

                // Get pressures and compute gravity contributions,
                // to decide upwind directions.
                PhaseVec phase_p[2];
                PhaseVec gravcontrib[2];
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        // Pressures
                        phase_p[j] = cell_pressure[c[j]];
                        // Gravity contribution.
                        if (nonzero_gravity) {
                            Vec cc = grid.cellCentroid(c[j]);
                            cc -= fc;
                            gravcontrib[j] = fluid.phaseDensities(&cellA[np*nc*c[j]]);;
                            gravcontrib[j] *= (cc*gravity);
                        } else {
                            gravcontrib[j] = 0.0;
                        }
                    } else {
                        // Pressures
                        phase_p[j] = face_pressure[face];
                        // Gravity contribution.
                        gravcontrib[j] = 0.0;
                    }
                }

                // Gravity contribution.
                //    gravcapf = rho_1*g*(z_1 - z_12) - rho_2*g*(z_2 - z_12)
                // where _1 and _2 refers to two neigbour cells, z is the
                // z coordinate of the centroid, and z_12 is the face centroid.
                for (int phase = 0; phase < np; ++phase) {
                    gravcapf[np*face + phase] = gravcontrib[0][phase] - gravcontrib[1][phase];
                }

                // Now we can easily find the upwind direction for every phase,
                // we can also tell which boundary faces are inflow bdys.


                PhaseVec phase_mob[2];
                PhaseMat phasemob_deriv[2];
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        // Pressures, mobilities.
                        phase_p[j] = cell_pressure[c[j]];
                        phase_mob[j] = phasemobc[c[j]];
                        phasemob_deriv[j] = phasemobc_deriv[c[j]];
                        // Gravity contribution.
                        if (nonzero_gravity) {
                            Vec cc = grid.cellCentroid(c[j]);
                            cc -= fc;
                            gravcontrib[j] = fluid.phaseDensities(&cellA[np*nc*c[j]]);;
                            gravcontrib[j] *= (cc*gravity);
                        } else {
                            gravcontrib[j] = 0.0;
                        }
                    } else {
                        // Pressures, mobilities.
                        phase_p[j] = face_pressure[face];
                        FluidStateBlackoil bdy_state = fluid.computeState(face_pressure[face], bdy_z);
                        phase_mob[j] = bdy_state.mobility_;
                        phasemob_deriv[j] = bdy_state.dmobility_;
                        // Gravity contribution.
                        gravcontrib[j] = 0.0;
                    }
                }

                // Compute face_z, which is averaged from the cells, unless on outflow or noflow bdy.
                CompVec face_z(0.0);
                double face_z_factor = 0.5;
                double pot_l[2] = { phase_p[0][Liquid], phase_p[1][Liquid] + gravcapf[np*face + Liquid] };
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        face_z += cell_z[c[j]];
                    } else if (pot_l[j] > pot_l[(j + 1)%2]) {
                        // Inflow boundary.
                        face_z += bdy_z;
                    } else {
                        // For outflow or noflow boundaries, only cell z is used.
                        face_z_factor = 1.0;
                        // Also, make sure the boundary data are not used for mobilities.
                        phase_p[j] = -1e100;
                    }
                }
                face_z *= face_z_factor;

                // Computing upwind mobilities and derivatives
                for (int phase = 0; phase < np; ++phase) {
                    double pot[2] = { phase_p[0][phase], phase_p[1][phase] + gravcapf[np*face + phase] };
                    if (pot[0] == pot[1]) {
                        // Average.
                        double aver = 0.5*(phase_mob[0][phase] + phase_mob[1][phase]);
                        phasemobf[np*face + phase] = aver;
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[np*np*face + np*phase + p2] = phasemob_deriv[0][phase][p2]
                                + phasemob_deriv[1][phase][p2];
                        }
                    } else {
                        // Upwind.
                        int upwind = pot[0] > pot[1] ? 0 : 1;
                        phasemobf[np*face + phase] = phase_mob[upwind][phase];
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[np*np*face + np*phase + p2] = phasemob_deriv[upwind][phase][p2];
                        }
                    }
                }

                // Find faceA.
                FluidStateBlackoil face_state = fluid.computeState(face_pressure[face], face_z);
                std::copy(face_state.phase_to_comp_, face_state.phase_to_comp_ + nc*np, &faceA[face*nc*np]);
            }
        }
    };


} // namespace Opm

#endif // OPM_BLACKOILFLUID_HEADER_INCLUDED
