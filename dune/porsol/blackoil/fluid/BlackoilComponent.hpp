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

#ifndef OPM_BLACKOILCOMPONENT_HEADER_INCLUDED
#define OPM_BLACKOILCOMPONENT_HEADER_INCLUDED


#include <dune/common/stdstreams.hh>

namespace Dumux
{

/*!
 * \ingroup Components

 * \brief
 *     A component class for the black oil model, intended to be used
 *     for all three components.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class BlackoilComponent
{
public:
    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    {
        return "BlackoilComponent";
    }

    /*!
     * \brief The molar mass in [kg] of the component.
     */
    static Scalar molarMass()
    { DUNE_THROW(Dune::NotImplemented, "Component::molarMass()"); }

    /*!
     * \brief Returns the critical temperature in [K] of the component.
     */
    static Scalar criticalTemperature()
    { DUNE_THROW(Dune::NotImplemented, "Component::criticalTemperature()"); }

    /*!
     * \brief Returns the critical pressure in [Pa] of the component.
     */
    static Scalar criticalPressure()
    { DUNE_THROW(Dune::NotImplemented, "Component::criticalPressure()"); }

    /*!
     * \brief Returns the temperature in [K] at the component's triple point.
     */
    static Scalar tripleTemperature()
    { DUNE_THROW(Dune::NotImplemented, "Component::tripleTemperature()"); }

    /*!
     * \brief Returns the pressure in [Pa] at the component's triple point.
     */
    static Scalar triplePressure()
    { DUNE_THROW(Dune::NotImplemented, "Component::triplePressure()"); }

    /*!
     * \brief The vapor pressure in [Pa] of the component at a given
     *        temperature in [K].
     *
     * \param T temperature of the component in [K]
     */
    static Scalar vaporPressure(Scalar T)
    { DUNE_THROW(Dune::NotImplemented, "Component::vaporPressure()"); }

    /*!
     * \brief The density in [kg/m^3] of the component at a given pressure in [Pa] and temperature in [K].
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::density()"); }

    /*!
     * \brief The density [kg/m^3] of the liquid component at a given pressure in [Pa] and temperature in [K].
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::density()"); }

    /*!
     * \brief Specific enthalpy [J/kg] of the pure component in gas.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasEnthalpy()"); }

    /*!
     * \brief Specific enthalpy [J/kg] of the pure component in liquid.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidEnthalpy()"); }

    /*!
     * \brief Specific internal energy [J/kg] of the pure component in gas.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasInternalEnergy()"); }

    /*!
     * \brief Specific internal energy [J/kg] of pure the pure component in liquid.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidInternalEnergy()"); }

    /*!
     * \brief The dynamic viscosity [Pa*s] of the pure component at a given pressure in [Pa] and temperature in [K].
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasViscosity()"); }

    /*!
     * \brief The dynamic liquid viscosity [Pa*s] of the pure component.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidViscosity()"); }

};

} // end namepace

#endif // OPM_BLACKOILCOMPONENT_HEADER_INCLUDED
