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
            // FluidSystemBlackoil<>::init(parser);
            pvt_.init(parser);
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
            // FluidSystemBlackoil<>::computeEquilibrium(state); // Sets everything but relperm and mobility.
            computeEquilibrium(state);
            FluidMatrixInteractionBlackoil<double>::kr(state.relperm_, fmi_params_, state.saturation_, state.temperature_);
            FluidMatrixInteractionBlackoil<double>::dkr(state.drelperm_, fmi_params_, state.saturation_, state.temperature_);
            for (int phase = 0; phase < numPhases; ++phase) {
                state.mobility_[phase] = state.relperm_[phase]/state.viscosity_[phase];
                for (int p2 = 0; p2 < numPhases; ++p2) {
                    // Ignoring pressure variation in viscosity for this one.
                    state.dmobility_[phase][p2] = state.drelperm_[phase][p2]/state.viscosity_[phase];
                }
            }
            return state;
        }

        void computeManyStates(const std::vector<PhaseVec>& phase_pressure,
                               const std::vector<CompVec>& z,
                               ManyFluidStatesBlackoil& state) const
        {
            int num = phase_pressure.size();
            state.temperature_.clear();
            state.temperature_.resize(num, 300.0);
            state.phase_pressure_ = phase_pressure;
            state.surface_volume_ = z;
            // FluidSystemBlackoil<>::computeManyEquilibria(state); // Sets everything but relperm and mobility.
            computeManyEquilibria(state);
            state.relperm_.resize(num);
            state.drelperm_.resize(num);
            state.mobility_.resize(num);
            state.dmobility_.resize(num);
#pragma omp parallel for
            for (int i = 0; i < num; ++i) {
                FluidMatrixInteractionBlackoil<double>::kr(state.relperm_[i], fmi_params_,
                                                           state.saturation_[i], state.temperature_[i]);
                FluidMatrixInteractionBlackoil<double>::dkr(state.drelperm_[i], fmi_params_,
                                                            state.saturation_[i], state.temperature_[i]);
                for (int phase = 0; phase < numPhases; ++phase) {
                    state.mobility_[i][phase] = state.relperm_[i][phase]/state.viscosity_[i][phase];
                    for (int p2 = 0; p2 < numPhases; ++p2) {
                        // Ignoring pressure variation in viscosity for this one.
                        state.dmobility_[i][phase][p2] = state.drelperm_[i][phase][p2]/state.viscosity_[i][phase];
                    }
                }
            }
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


        // ---- New interface ----

        /// Input: p, z
        /// Output: B, R
        template <class States>
        void computeBAndR(States& states) const
        {
            const std::vector<PhaseVec>& p = states.phase_pressure;
            const std::vector<CompVec>& z = states.surface_volume_density;
            std::vector<PhaseVec>& B = states.formation_volume_factor;
            std::vector<PhaseVec>& R = states.solution_factor;
            pvt_.B(p, z, B);
            pvt_.R(p, z, R);
        }

        /// Input: p, z
        /// Output: B, R, mu
        template <class States>
        void computePvtNoDerivs(States& states) const
        {
            computeBAndR(states);
            const std::vector<PhaseVec>& p = states.phase_pressure;
            const std::vector<CompVec>& z = states.surface_volume_density;
            std::vector<PhaseVec>& mu = states.viscosity;
            pvt_.getViscosity(p, z, mu);
        }

        /// Input: p, z
        /// Output: B, dB/dp, R, dR/dp, mu
        template <class States>
        void computePvt(States& states) const
        {
            const std::vector<PhaseVec>& p = states.phase_pressure;
            const std::vector<CompVec>& z = states.surface_volume_density;
            std::vector<PhaseVec>& B = states.formation_volume_factor;
            std::vector<PhaseVec>& dB = states.formation_volume_factor_deriv;
            std::vector<PhaseVec>& R = states.solution_factor;
            std::vector<PhaseVec>& dR = states.solution_factor_deriv;
            std::vector<PhaseVec>& mu = states.viscosity;
            pvt_.dBdp(p, z, B, dB);
            pvt_.dRdp(p, z, R, dR);
            pvt_.getViscosity(p, z, mu);
        }

        /// Input: B, R
        /// Output: A
        template <class States>
        void computeStateMatrix(States& states) const
        {
            int num = states.formation_volume_factor.size();
            states.state_matrix.resize(num);
#pragma omp parallel for
            for (int i = 0; i < num; ++i) {
                const PhaseVec& B = states.formation_volume_factor[i];
                const PhaseVec& R = states.solution_factor[i];
                PhaseToCompMatrix& At = states.state_matrix[i];
                // Set the A matrix (A = RB^{-1})
                // Using A transposed (At) since we really want Fortran ordering:
                // ultimately that is what the opmpressure C library expects.
                At = 0.0;
                At[Aqua][Water] = 1.0/B[Aqua];
                At[Vapour][Gas] = 1.0/B[Vapour];
                At[Liquid][Gas] = R[Liquid]/B[Liquid];
                At[Vapour][Oil] = R[Vapour]/B[Vapour];
                At[Liquid][Oil] = 1.0/B[Liquid];
            }
        }


        /// Input: z, B, dB/dp, R, dR/dp
        /// Output: A, u, sum(u), s, c, cT, ex
        template <class States>
        void computePvtDepending(States& states) const
        {
            int num = states.formation_volume_factor.size();
            states.state_matrix.resize(num);
            states.phase_volume_density.resize(num);
            states.total_phase_volume_density.resize(num);
            states.saturation.resize(num);
            states.phase_compressibility.resize(num);
            states.total_compressibility.resize(num);
            states.experimental_term.resize(num);
#pragma omp parallel for
            for (int i = 0; i < num; ++i) {
                const CompVec& z = states.surface_volume_density[i];
                const PhaseVec& B = states.formation_volume_factor[i];
                const PhaseVec& dB = states.formation_volume_factor_deriv[i];
                const PhaseVec& R = states.solution_factor[i];
                const PhaseVec& dR = states.solution_factor_deriv[i];
                PhaseToCompMatrix& At = states.state_matrix[i];
                PhaseVec& u = states.phase_volume_density[i];
                double& tot_phase_vol_dens = states.total_phase_volume_density[i];
                PhaseVec& s = states.saturation[i];
                PhaseVec& cp = states.phase_compressibility[i];
                double& tot_comp = states.total_compressibility[i];
                double& exp_term = states.experimental_term[i];
                computeSingleEquilibrium(B, dB, R, dR, z,
                                         At, u, tot_phase_vol_dens,
                                         s, cp, tot_comp, exp_term);
            }
        }

        /// Input: s, mu
        /// Output: kr, lambda
        template <class States>
        void computeMobilitiesNoDerivs(States& states) const
        {
            int num = states.saturation.size();
            states.relperm.resize(num);
            states.mobility.resize(num);
#pragma omp parallel for
            for (int i = 0; i < num; ++i) {
                const CompVec& s = states.saturation[i];
                const PhaseVec& mu = states.viscosity[i];
                PhaseVec& kr = states.relperm[i];
                PhaseVec& lambda = states.mobility[i];
                FluidMatrixInteractionBlackoil<double>::kr(kr, fmi_params_, s, 300.0);
                for (int phase = 0; phase < numPhases; ++phase) {
                    lambda[phase] = kr[phase]/mu[phase];
                }

            }
        }

        /// Input: s, mu
        /// Output: kr, dkr/ds, lambda, dlambda/ds
        template <class States>
        void computeMobilities(States& states) const
        {
            int num = states.saturation.size();
            states.relperm.resize(num);
            states.relperm_deriv.resize(num);
            states.mobility.resize(num);
            states.mobility_deriv.resize(num);
#pragma omp parallel for
            for (int i = 0; i < num; ++i) {
                const CompVec& s = states.saturation[i];
                const PhaseVec& mu = states.viscosity[i];
                PhaseVec& kr = states.relperm[i];
                PhaseJacobian& dkr = states.relperm_deriv[i];
                PhaseVec& lambda = states.mobility[i];
                PhaseJacobian& dlambda = states.mobility_deriv[i];
                FluidMatrixInteractionBlackoil<double>::kr(kr, fmi_params_, s, 300.0);
                FluidMatrixInteractionBlackoil<double>::dkr(dkr, fmi_params_, s, 300.0);
                for (int phase = 0; phase < numPhases; ++phase) {
                    lambda[phase] = kr[phase]/mu[phase];
                    for (int p2 = 0; p2 < numPhases; ++p2) {
                        // Ignoring pressure variation in viscosity for this one.
                        dlambda[phase][p2] = dkr[phase][p2]/mu[phase];
                    }
                }

            }
        }


    private:
        BlackoilPVT pvt_;
        FluidMatrixInteractionBlackoilParams<double> fmi_params_;
        CompVec surface_densities_;


    /*!
     * \brief Assuming the surface volumes and the pressures of all
     *        phases are known, compute everything except relperm and
     *        mobility.
     */
    template <class FluidState>
    void computeEquilibrium(FluidState& fluid_state) const
    {
        // Get B and R factors.
        const PhaseVec& p = fluid_state.phase_pressure_;
        const CompVec& z = fluid_state.surface_volume_;
        PhaseVec& B = fluid_state.formation_volume_factor_;
        B[Aqua]   = pvt_.B(p[Aqua],   z, Aqua);
        B[Vapour] = pvt_.B(p[Vapour], z, Vapour);
        B[Liquid] = pvt_.B(p[Liquid], z, Liquid);
        PhaseVec& R = fluid_state.solution_factor_; 
        R[Aqua]   = 0.0;
        R[Vapour] = pvt_.R(p[Vapour], z, Vapour);
        R[Liquid] = pvt_.R(p[Liquid], z, Liquid);
        PhaseVec dB;
        dB[Aqua]   = pvt_.dBdp(p[Aqua],   z, Aqua);
        dB[Vapour] = pvt_.dBdp(p[Vapour], z, Vapour);
        dB[Liquid] = pvt_.dBdp(p[Liquid], z, Liquid);
        PhaseVec dR;
        dR[Aqua]   = 0.0;
        dR[Vapour] = pvt_.dRdp(p[Vapour], z, Vapour);
        dR[Liquid] = pvt_.dRdp(p[Liquid], z, Liquid);

        // Convenience vars.
        PhaseToCompMatrix& At = fluid_state.phase_to_comp_;
        PhaseVec& u = fluid_state.phase_volume_density_;
        double& tot_phase_vol_dens = fluid_state.total_phase_volume_density_;
        PhaseVec& s = fluid_state.saturation_;
        PhaseVec& cp = fluid_state.phase_compressibility_;
        double& tot_comp = fluid_state.total_compressibility_;
        double& exp_term = fluid_state.experimental_term_;

        computeSingleEquilibrium(B, dB, R, dR, z,
                                 At, u, tot_phase_vol_dens,
                                 s, cp, tot_comp, exp_term);

        // Compute viscosities.
        PhaseVec& mu = fluid_state.viscosity_;
        mu[Aqua]   = pvt_.getViscosity(p[Aqua],   z, Aqua);
        mu[Vapour] = pvt_.getViscosity(p[Vapour], z, Vapour);
        mu[Liquid] = pvt_.getViscosity(p[Liquid], z, Liquid);
    }

    /*!
     * \brief Assuming the surface volumes and the pressures of all
     *        phases are known, compute everything except relperm and
     *        mobility.
     */
    template <class ManyFluidStates>
    void computeManyEquilibria(ManyFluidStates& fluid_state) const
    {
        // Get B and R factors, viscosities. Vectorized at lower level.
        const std::vector<PhaseVec>& pv = fluid_state.phase_pressure_;
        const std::vector<CompVec>& zv = fluid_state.surface_volume_;
        std::vector<PhaseVec>& Bv = fluid_state.formation_volume_factor_;
        std::vector<PhaseVec> dBv;
        pvt_.dBdp(pv, zv, Bv, dBv);
        std::vector<PhaseVec>& Rv = fluid_state.solution_factor_;
        std::vector<PhaseVec> dRv;
        pvt_.dRdp(pv, zv, Rv, dRv);
        std::vector<PhaseVec>& muv = fluid_state.viscosity_;
        pvt_.getViscosity(pv, zv, muv);

        // The rest is vectorized in this function.
        int num = pv.size();
        fluid_state.phase_to_comp_.resize(num);
        fluid_state.phase_volume_density_.resize(num);
        fluid_state.total_phase_volume_density_.resize(num);
        fluid_state.phase_to_comp_.resize(num);
        fluid_state.saturation_.resize(num);
        fluid_state.phase_compressibility_.resize(num);
        fluid_state.total_compressibility_.resize(num);
        fluid_state.experimental_term_.resize(num);
#pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            // Convenience vars.
            const PhaseVec& B = Bv[i];
            const PhaseVec& dB = dBv[i];
            const PhaseVec& R = Rv[i];
            const PhaseVec& dR = dRv[i];
            const CompVec& z = zv[i];

            PhaseToCompMatrix& At = fluid_state.phase_to_comp_[i];
            PhaseVec& u = fluid_state.phase_volume_density_[i];
            double& tot_phase_vol_dens = fluid_state.total_phase_volume_density_[i];
            PhaseVec& s = fluid_state.saturation_[i];
            PhaseVec& cp = fluid_state.phase_compressibility_[i];
            double& tot_comp = fluid_state.total_compressibility_[i];
            double& exp_term = fluid_state.experimental_term_[i];

            computeSingleEquilibrium(B, dB, R, dR, z,
                                     At, u, tot_phase_vol_dens,
                                     s, cp, tot_comp, exp_term);
        }
    }

    static void computeSingleEquilibrium(const PhaseVec& B,
                                         const PhaseVec& dB,
                                         const PhaseVec& R,
                                         const PhaseVec& dR,
                                         const CompVec& z,
                                         PhaseToCompMatrix& At,
                                         PhaseVec& u,
                                         double& tot_phase_vol_dens,
                                         PhaseVec& s,
                                         PhaseVec& cp,
                                         double& tot_comp,
                                         double& exp_term)
    {
        // Set the A matrix (A = RB^{-1})
        // Using At since we really want Fortran ordering
        // (since ultimately that is what the opmpressure
        //  C library expects).
        // PhaseToCompMatrix& At = fluid_state.phase_to_comp_[i];
        At = 0.0;
        At[Aqua][Water] = 1.0/B[Aqua];
        At[Vapour][Gas] = 1.0/B[Vapour];
        At[Liquid][Gas] = R[Liquid]/B[Liquid];
        At[Vapour][Oil] = R[Vapour]/B[Vapour];
        At[Liquid][Oil] = 1.0/B[Liquid];

        // Update phase volumes. This is the same as multiplying with A^{-1}
        // PhaseVec& u = fluid_state.phase_volume_density_[i];
        double detR = 1.0 - R[Vapour]*R[Liquid];
        u[Aqua] = B[Aqua]*z[Water];
        u[Vapour] = B[Vapour]*(z[Gas] - R[Liquid]*z[Oil])/detR;
        u[Liquid] = B[Liquid]*(z[Oil] - R[Vapour]*z[Gas])/detR;
        tot_phase_vol_dens = u[Aqua] + u[Vapour] + u[Liquid];

        // PhaseVec& s = fluid_state.saturation_[i];
        for (int phase = 0; phase < numPhases; ++phase) {
            s[phase] = u[phase]/tot_phase_vol_dens;
        }

        // Phase compressibilities.
        // PhaseVec& cp = fluid_state.phase_compressibility_[i];
        // Set the derivative of the A matrix (A = RB^{-1})
        PhaseToCompMatrix dAt(0.0);
        dAt[Aqua][Water] = -dB[Aqua]/(B[Aqua]*B[Aqua]);
        dAt[Vapour][Gas] = -dB[Vapour]/(B[Vapour]*B[Vapour]);
        dAt[Liquid][Oil] = -dB[Liquid]/(B[Liquid]*B[Liquid]); // Different order than above.
        dAt[Liquid][Gas] = dAt[Liquid][Oil]*R[Liquid] + dR[Liquid]/B[Liquid];
        dAt[Vapour][Oil] = dAt[Vapour][Gas]*R[Vapour] + dR[Vapour]/B[Vapour];

        PhaseToCompMatrix Ait;
        Dune::FMatrixHelp::invertMatrix(At, Ait);

        PhaseToCompMatrix Ct;
        Dune::FMatrixHelp::multMatrix(dAt, Ait, Ct);

        cp[Aqua] = Ct[Aqua][Water];
        cp[Liquid] = Ct[Liquid][Oil] + Ct[Liquid][Gas];
        cp[Vapour] = Ct[Vapour][Gas] + Ct[Vapour][Oil];
        tot_comp = cp*s;

        // Experimental term.
        PhaseVec tmp1, tmp2, tmp3;
        Ait.mtv(z, tmp1);
        dAt.mtv(tmp1, tmp2);
        Ait.mtv(tmp2, tmp3);
        exp_term = tmp3[Aqua] + tmp3[Liquid] + tmp3[Gas];
    }


    };





    struct AllFluidData : public BlackoilDefs
    {
        // Per-cell data
        AllFluidStates cell_data;
        std::vector<double> voldiscr;
        std::vector<double> relvoldiscr;

        // Per-face data.
        std::vector<PhaseToCompMatrix> faceA;        // A = RB^{-1}. Fortran ordering.
        std::vector<PhaseVec> phasemobf;             // Phase mobilities. Flat storage (numPhases per face).
        std::vector<PhaseJacobian> phasemobf_deriv;  // Phase mobility derivatives. Flat storage (numPhases^2 per face).
        std::vector<PhaseVec> gravcapf;              // Gravity (\rho g \delta z-ish) contribution per face, flat storage.

    public:
        template <class Grid, class Rock>
        void computeNew(const Grid& grid,
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
            BOOST_STATIC_ASSERT(np == 3);

            // p, z -> B, dB, R, dR, mu, A, dA, u, sum(u), s, c, cT, ex, kr, dkr, lambda, dlambda.
            cell_data.phase_pressure = cell_pressure;
            cell_data.surface_volume_density = cell_z;
            fluid.computePvt(cell_data);
            fluid.computePvtDepending(cell_data);
            fluid.computeMobilities(cell_data);

            // Compute volume discrepancies.
            voldiscr.resize(num_cells);
            relvoldiscr.resize(num_cells);
#pragma omp parallel for
            for (int cell = 0; cell < num_cells; ++cell) {
                double pv = rock.porosity(cell)*grid.cellVolume(cell);
                voldiscr[cell] = (cell_data.total_phase_volume_density[cell] - 1.0)*pv/dt;
                relvoldiscr[cell] = std::fabs(cell_data.total_phase_volume_density[cell] - 1.0);
            }

            // Compute face properties.
            faceA.resize(num_faces);
            phasemobf.resize(num_faces);
            phasemobf_deriv.resize(num_faces);
            gravcapf.resize(num_faces);
            // std::vector<CompVec> face_z_vec(num_faces);
#pragma omp parallel for
            for (int face = 0; face < num_faces; ++face) {
            // Obtain properties from both sides of the face.
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
                            Vec cdiff = fc;
                            cdiff -= grid.cellCentroid(c[j]);
                            gravcontrib[j] = fluid.phaseDensities(&cell_data.state_matrix[c[j]][0][0]);
                            gravcontrib[j] *= (cdiff*gravity);
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

                // Gravity contribution:
                //    gravcapf = rho_1*g*(z_12 - z_1) - rho_2*g*(z_12 - z_2)
                // where _1 and _2 refers to two neigbour cells, z is the
                // z coordinate of the centroid, and z_12 is the face centroid.
                // Also compute the potentials.
                PhaseVec pot[2];
                for (int phase = 0; phase < np; ++phase) {
                    gravcapf[face][phase] = gravcontrib[0][phase] - gravcontrib[1][phase];
                    pot[0][phase] = phase_p[0][phase] + gravcapf[face][phase];
                    pot[1][phase] = phase_p[1][phase];
                }

                // Now we can easily find the upwind direction for every phase,
                // we can also tell which boundary faces are inflow bdys.

                // Compute face_z, which is averaged from the cells, unless on outflow or noflow bdy.
                // Get mobilities and derivatives.
                CompVec face_z(0.0);
                double face_z_factor = 0.5;
                PhaseVec phase_mob[2];
                PhaseJacobian phasemob_deriv[2];
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        face_z += cell_z[c[j]];
                        phase_mob[j] = cell_data.mobility[c[j]];
                        phasemob_deriv[j] = cell_data.mobility_deriv[c[j]];
                    } else if (pot[j][Liquid] > pot[(j+1)%2][Liquid]) {
                        // Inflow boundary.
                        face_z += bdy_z;
                        FluidStateBlackoil bdy_state = fluid.computeState(face_pressure[face], bdy_z);
                        phase_mob[j] = bdy_state.mobility_;
                        phasemob_deriv[j] = bdy_state.dmobility_;
                    } else {
                        // For outflow or noflow boundaries, only cell z is used.
                        face_z_factor = 1.0;
                        // Also, make sure the boundary data are not used for mobilities.
                        pot[j] = -1e100;
                    }
                }
                face_z *= face_z_factor;
                // face_z_vec[face] = face_z;

                // Computing upwind mobilities and derivatives
                for (int phase = 0; phase < np; ++phase) {
                    if (pot[0][phase] == pot[1][phase]) {
                        // Average.
                        double aver = 0.5*(phase_mob[0][phase] + phase_mob[1][phase]);
                        phasemobf[face][phase] = aver;
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[face][phase][p2] = phasemob_deriv[0][phase][p2]
                                + phasemob_deriv[1][phase][p2];
                        }
                    } else {
                        // Upwind.
                        int upwind = pot[0][phase] > pot[1][phase] ? 0 : 1;
                        phasemobf[face][phase] = phase_mob[upwind][phase];
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[face][phase][p2] = phasemob_deriv[upwind][phase][p2];
                        }
                    }
                }
                // Find faceA.
                FluidStateBlackoil face_state = fluid.computeState(face_pressure[face], face_z);
                std::copy(&face_state.phase_to_comp_[0][0], &face_state.phase_to_comp_[0][0] + nc*np, &faceA[face][0][0]);
            }

            // Find faceA. Slower, since it does too much.
//             ManyFluidStatesBlackoil face_state;
//             fluid.computeManyStates(face_pressure, face_z_vec, face_state);
//             const double* fA = &face_state.phase_to_comp_[0][0][0];
//             std::copy(fA, fA + num_faces*nc*np, &faceA[0]);
        }


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
        std::vector<PhaseJacobian> phasemobc_deriv; // Phase mobilities derivatives per cell.
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
                std::copy(&state.phase_to_comp_[0][0], &state.phase_to_comp_[0][0] + nc*np, &cellA[cell*nc*np]);
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
                            Vec cdiff = fc;
                            cdiff -= grid.cellCentroid(c[j]);
                            gravcontrib[j] = fluid.phaseDensities(&cellA[np*nc*c[j]]);;
                            gravcontrib[j] *= (cdiff*gravity);
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

                // Gravity contribution:
                //    gravcapf = rho_1*g*(z_12 - z_1) - rho_2*g*(z_12 - z_2)
                // where _1 and _2 refers to two neigbour cells, z is the
                // z coordinate of the centroid, and z_12 is the face centroid.
                // Also compute the potentials.
                PhaseVec pot[2];
                for (int phase = 0; phase < np; ++phase) {
                    gravcapf[np*face + phase] = gravcontrib[0][phase] - gravcontrib[1][phase];
                    pot[0][phase] = phase_p[0][phase] + gravcapf[np*face + phase];
                    pot[1][phase] = phase_p[1][phase];
                }

                // Now we can easily find the upwind direction for every phase,
                // we can also tell which boundary faces are inflow bdys.

                // Compute face_z, which is averaged from the cells, unless on outflow or noflow bdy.
                // Get mobilities and derivatives.
                CompVec face_z(0.0);
                double face_z_factor = 0.5;
                PhaseVec phase_mob[2];
                PhaseJacobian phasemob_deriv[2];
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        face_z += cell_z[c[j]];
                        phase_mob[j] = phasemobc[c[j]];
                        phasemob_deriv[j] = phasemobc_deriv[c[j]];
                    } else if (pot[j][Liquid] > pot[(j+1)%2][Liquid]) {
                        // Inflow boundary.
                        face_z += bdy_z;
                        FluidStateBlackoil bdy_state = fluid.computeState(face_pressure[face], bdy_z);
                        phase_mob[j] = bdy_state.mobility_;
                        phasemob_deriv[j] = bdy_state.dmobility_;
                    } else {
                        // For outflow or noflow boundaries, only cell z is used.
                        face_z_factor = 1.0;
                        // Also, make sure the boundary data are not used for mobilities.
                        pot[j] = -1e100;
                    }
                }
                face_z *= face_z_factor;

                // Computing upwind mobilities and derivatives
                for (int phase = 0; phase < np; ++phase) {
                    if (pot[0][phase] == pot[1][phase]) {
                        // Average.
                        double aver = 0.5*(phase_mob[0][phase] + phase_mob[1][phase]);
                        phasemobf[np*face + phase] = aver;
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[np*np*face + np*phase + p2] = phasemob_deriv[0][phase][p2]
                                + phasemob_deriv[1][phase][p2];
                        }
                    } else {
                        // Upwind.
                        int upwind = pot[0][phase] > pot[1][phase] ? 0 : 1;
                        phasemobf[np*face + phase] = phase_mob[upwind][phase];
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[np*np*face + np*phase + p2] = phasemob_deriv[upwind][phase][p2];
                        }
                    }
                }

                // Find faceA.
                FluidStateBlackoil face_state = fluid.computeState(face_pressure[face], face_z);
                std::copy(&face_state.phase_to_comp_[0][0], &face_state.phase_to_comp_[0][0] + nc*np, &faceA[face*nc*np]);
            }
        }




    public:
        template <class Grid, class Rock>
        void computeNew(const Grid& grid,
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

            // Compute and copy properties.
            ManyFluidStatesBlackoil state;
            fluid.computeManyStates(cell_pressure, cell_z, state);
            totcompr = state.total_compressibility_;
            totphasevol_density = state.total_phase_volume_density_;
            saturation = state.saturation_;
            rel_perm = state.relperm_;
            viscosity = state.viscosity_;
            phasemobc = state.mobility_;
            phasemobc_deriv = state.dmobility_;
            std::copy(&state.phase_to_comp_[0][0][0],
                      &state.phase_to_comp_[0][0][0] + num_cells*nc*np,
                      &cellA[0]);
            expjacterm = state.experimental_term_;

#pragma omp parallel for
            for (int cell = 0; cell < num_cells; ++cell) {
                double pv = rock.porosity(cell)*grid.cellVolume(cell);
                voldiscr[cell] = (totphasevol_density[cell] - 1.0)*pv/dt;
                relvoldiscr[cell] = std::fabs(totphasevol_density[cell] - 1.0);

                // Fractional flow must be calculated.
                double total_mobility = 0.0;
                for (int phase = 0; phase < numPhases; ++phase) {
                    total_mobility += state.mobility_[cell][phase];
                }
                frac_flow[cell] = state.mobility_[cell];
                frac_flow[cell] /= total_mobility;
            }

            // std::vector<CompVec> face_z_vec(num_faces);
#pragma omp parallel for
            for (int face = 0; face < num_faces; ++face) {
            // Obtain properties from both sides of the face.
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
                            Vec cdiff = fc;
                            cdiff -= grid.cellCentroid(c[j]);
                            gravcontrib[j] = fluid.phaseDensities(&cellA[np*nc*c[j]]);;
                            gravcontrib[j] *= (cdiff*gravity);
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

                // Gravity contribution:
                //    gravcapf = rho_1*g*(z_12 - z_1) - rho_2*g*(z_12 - z_2)
                // where _1 and _2 refers to two neigbour cells, z is the
                // z coordinate of the centroid, and z_12 is the face centroid.
                // Also compute the potentials.
                PhaseVec pot[2];
                for (int phase = 0; phase < np; ++phase) {
                    gravcapf[np*face + phase] = gravcontrib[0][phase] - gravcontrib[1][phase];
                    pot[0][phase] = phase_p[0][phase] + gravcapf[np*face + phase];
                    pot[1][phase] = phase_p[1][phase];
                }

                // Now we can easily find the upwind direction for every phase,
                // we can also tell which boundary faces are inflow bdys.

                // Compute face_z, which is averaged from the cells, unless on outflow or noflow bdy.
                // Get mobilities and derivatives.
                CompVec face_z(0.0);
                double face_z_factor = 0.5;
                PhaseVec phase_mob[2];
                PhaseJacobian phasemob_deriv[2];
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        face_z += cell_z[c[j]];
                        phase_mob[j] = phasemobc[c[j]];
                        phasemob_deriv[j] = phasemobc_deriv[c[j]];
                    } else if (pot[j][Liquid] > pot[(j+1)%2][Liquid]) {
                        // Inflow boundary.
                        face_z += bdy_z;
                        FluidStateBlackoil bdy_state = fluid.computeState(face_pressure[face], bdy_z);
                        phase_mob[j] = bdy_state.mobility_;
                        phasemob_deriv[j] = bdy_state.dmobility_;
                    } else {
                        // For outflow or noflow boundaries, only cell z is used.
                        face_z_factor = 1.0;
                        // Also, make sure the boundary data are not used for mobilities.
                        pot[j] = -1e100;
                    }
                }
                face_z *= face_z_factor;
                // face_z_vec[face] = face_z;

                // Computing upwind mobilities and derivatives
                for (int phase = 0; phase < np; ++phase) {
                    if (pot[0][phase] == pot[1][phase]) {
                        // Average.
                        double aver = 0.5*(phase_mob[0][phase] + phase_mob[1][phase]);
                        phasemobf[np*face + phase] = aver;
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[np*np*face + np*phase + p2] = phasemob_deriv[0][phase][p2]
                                + phasemob_deriv[1][phase][p2];
                        }
                    } else {
                        // Upwind.
                        int upwind = pot[0][phase] > pot[1][phase] ? 0 : 1;
                        phasemobf[np*face + phase] = phase_mob[upwind][phase];
                        for (int p2 = 0; p2 < numPhases; ++p2) {
                            phasemobf_deriv[np*np*face + np*phase + p2] = phasemob_deriv[upwind][phase][p2];
                        }
                    }
                }
                // Find faceA.
                FluidStateBlackoil face_state = fluid.computeState(face_pressure[face], face_z);
                std::copy(&face_state.phase_to_comp_[0][0], &face_state.phase_to_comp_[0][0] + nc*np, &faceA[face*nc*np]);
            }

            // Find faceA. Slower, since it does too much.
//             ManyFluidStatesBlackoil face_state;
//             fluid.computeManyStates(face_pressure, face_z_vec, face_state);
//             const double* fA = &face_state.phase_to_comp_[0][0][0];
//             std::copy(fA, fA + num_faces*nc*np, &faceA[0]);
        }


    };


} // namespace Opm

#endif // OPM_BLACKOILFLUID_HEADER_INCLUDED
