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

#ifndef OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED
#define OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED


#include ""


namespace Opm
{
/*!
 * \brief Calcultes the phase state from the primary variables in the
 *        blackoil model.
 */
class FluidStateBlackoil
{
public:
    enum { numPhases = 3 };
    enum { numComponents = 3};

    typedef int FluidSystem; // TODO
    typedef int PrimaryVars; // TODO
    typedef double Scalar;
    typedef FluidMatrixInteractionBlackoil<Scalar> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief Update the phase state from the primary variables.
     */
    void update(const PrimaryVariables &primaryVars,
                const MaterialLawParams &pcParams,
                Scalar temperature)
    {
        // calculate the temperature
        temperature_ = temperature;

        // calculate the saturations
        saturation_[numPhases - 1] = 1.0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            saturation_[phaseIdx] = primaryVars[S0Idx + phaseIdx];
            saturation_[numPhases - 1] -= saturation_[phaseIdx];
        }

        // let the material law calculate the capillary pressures
        MaterialLaw::pC(phasePressure_, 
                        pcParams, 
                        saturation_, 
                        temperature_);

        // convert to phase pressures
        Scalar sumPc = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            sumPc += phasePressure_[phaseIdx];
            phasePressure_[phaseIdx] = primaryVars[p0Idx] - sumPc;
        }
        
        updateComposition_(primaryVars);
    }

private:
    void updateComposition_(const PrimaryVariables &primaryVars)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            meanMolarMass_[phaseIdx] = 0;
        }

        // extract the mole fractions in both phases
        Scalar sumFuga = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            fugacity_[compIdx] = primaryVars[fug0Idx + compIdx];
            Valgrind::CheckDefined(fugacity_[compIdx]);
            sumFuga += fugacity_[compIdx];

            // convert the fugacities into mole fractions
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                moleFrac_[phaseIdx][compIdx] =
                    fugacity_[compIdx] /
                    FluidSystem::activityCoeff(phaseIdx,
                                               compIdx,
                                               temperature_,
                                               phasePressure(phaseIdx),
                                               *this);
                Valgrind::CheckDefined(moleFrac_[phaseIdx][compIdx]);

                // update the mean molar mass of each phase
                meanMolarMass_[phaseIdx]
                    += FluidSystem::molarMass(compIdx)*moleFrac_[phaseIdx][compIdx];
            }
        }

        // make sure that the primary variables do not represent an
        // unphysical state!
        if (sumFuga < 1e-30) {
            DUNE_THROW(NumericalProblem,
                       "Sum of component fugacities is too small: " << sumFuga);
        }

        // calculate the total concentration of each phase
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // make sure that the mean molar mass within a phase is
            // not "way off"
            if (std::abs(meanMolarMass_[phaseIdx]) < 1e-30)
                DUNE_THROW(NumericalProblem,
                           "Mean molar mass of phase "
                           << phaseIdx
                           << " is too small ("
                           << meanMolarMass_[phaseIdx]
                           << ")\n" << primaryVars);

            // calculate the total concentration of the phase
            phaseConcentration_[phaseIdx] =
                FluidSystem::phaseDensity(phaseIdx,
                                          temperature(),
                                          phasePressure(phaseIdx),
                                          *this)
                /
                meanMolarMass_[phaseIdx];
        }

        Valgrind::CheckDefined(fugacity_);
        Valgrind::CheckDefined(moleFrac_);
        Valgrind::CheckDefined(phaseConcentration_);
        Valgrind::CheckDefined(meanMolarMass_);
        Valgrind::CheckDefined(temperature_);
        Valgrind::CheckDefined(phasePressure_);
        Valgrind::CheckDefined(saturation_);
    }

public:
    /*!
     * \brief Returns the saturation of a phase.
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * \brief Returns a vector of all phase saturations.
     */
    const Scalar *saturations() const
    { return saturation_; }

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    { return moleFrac_[phaseIdx][compIdx]; }

    /*!
     * \brief Returns the total concentration of a phase [mol / m^3].
     *
     * This is equivalent to the sum of all component concentrations.
     */
    Scalar phaseConcentration(int phaseIdx) const
    { return phaseConcentration_[phaseIdx]; };

    /*!
     * \brief Returns the concentration of a component in a phase [mol / m^3].
     */
    Scalar concentration(int phaseIdx, int compIdx) const
    { return moleFrac_[phaseIdx][compIdx]*phaseConcentration_[phaseIdx]; };

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    {
        return
            moleFrac_[phaseIdx][compIdx]*
            FluidSystem::molarMass(compIdx)
            /
            meanMolarMass_[phaseIdx];
    }

    /*!
     * \brief Returns the density of a phase [kg / m^3].
     */
    Scalar density(int phaseIdx) const
    { return phaseConcentration_[phaseIdx]*meanMolarMass_[phaseIdx]; }

    /*!
     * \brief Returns mean molar mass of a phase [kg / mol].
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     */
    Scalar meanMolarMass(int phaseIdx) const
    { return meanMolarMass_[phaseIdx]; };

    /*!
     * \brief Returns the fugacity of a component [Pa].
     */
    Scalar fugacity(int compIdx) const
    { return fugacity_[compIdx]; }

    /*!
     * \brief Returns the pressure of a fluid phase [Pa].
     */
    Scalar phasePressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure [Pa]
     */
    Scalar capillaryPressure(int phaseIdx) const
    { return phasePressure_[0] - phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids [K].
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

public:
    Scalar temperature_;
    Scalar surface_volume_[numComponents];
    Scalar phase_pressure_[numPhases];
    Scalar phase_volume_[numPhases];
    Scalar saturation_[numPhases];
};

} // end namespace Opm


#endif // OPM_FLUIDSTATEBLACKOIL_HEADER_INCLUDED
