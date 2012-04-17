//===========================================================================
//
// File: Units.hpp
//
// Created: Thu Jul  2 09:19:08 2009
//
// Author(s): Halvor M Nilsen <hnil@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2011, 2012 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_UNITS_HEADER
#define OPENRS_UNITS_HEADER

namespace Opm
{
    namespace prefix
    /// Conversion prefix for units.
    {
        const double micro = 1.0e-6;
        const double milli = 1.0e-3;
        const double centi = 1.0e-2;
        const double deci  = 1.0e-1;
        const double kilo  = 1.0e3;
        const double mega  = 1.0e6;
        const double giga  = 1.0e9;
    } // namespace prefix

    namespace unit
    /// Definition of various units.
    /// All the units are defined in terms of international standard
    /// units (SI).  Example of use: We define a variable \c k which
    /// gives a permeability. We want to set \c k to \f$1\,mD\f$.
    /// \code
    /// using namespace Opm::unit
    /// double k = 0.001*darcy;
    /// \endcode
    /// We can also use one of the prefixes defined in Opm::prefix
    /// \code
    /// using namespace Opm::unit
    /// using namespace Opm::prefix
    /// double k = 1.0*milli*darcy;
    /// \endcode
    {
        ///\name Common powers
        /// @{
        inline double square(double v) { return v * v;     }
        inline double cubic (double v) { return v * v * v; }
        /// @}

        // --------------------------------------------------------------
        // Basic (fundamental) units and conversions
        // --------------------------------------------------------------

        /// \name Length
        /// @{
        const double meter =  1;
        const double inch  =  2.54 * prefix::centi*meter;
        const double feet  = 12    * inch;
        /// @}

        /// \name Time
        /// @{
        const double second =   1;
        const double minute =  60 * second;
        const double hour   =  60 * minute;
        const double day    =  24 * hour;
        const double year   = 365 * day;
        /// @}

        /// \name Volume
        /// @{
        const double gallon = 231 * cubic(inch);
        const double stb    =  42 * gallon;
        /// @}

        /// \name Mass
        /// @{
        const double kilogram = 1;
        // http://en.wikipedia.org/wiki/Pound_(mass)#Avoirdupois_pound
        const double pound    = 0.45359237 * kilogram;
        /// @}

        // --------------------------------------------------------------
        // Standardised constants
        // --------------------------------------------------------------

        /// \name Standardised constant
        /// @{
        const double gravity = 9.80665 * meter/square(second);
        /// @}

        // --------------------------------------------------------------
        // Derived units and conversions
        // --------------------------------------------------------------

        /// \name Force
        /// @{
        const double Newton = kilogram*meter / square(second); // == 1
        const double lbf    = pound * gravity; // Pound-force
        /// @}

        /// \name Pressure
        /// @{
        const double Pascal = Newton / square(meter); // == 1
        const double barsa  = 100000 * Pascal;
        const double atm    = 101325 * Pascal;
        const double psia   = lbf / square(inch);
        /// @}

        /// \name Viscosity
        /// @{
        const double Pas   = Pascal * second; // == 1
        const double Poise = prefix::deci*Pas;
        /// @}

        /// \name Permeability
        /// @{
        ///
        /// A porous medium with a permeability of 1 darcy permits a
        /// flow (flux) of \f$1\,cm^3/s\f$ of a fluid with viscosity
        /// \f$1\,cP\f$ (\f$1\,mPa\cdot s\f$) under a pressure
        /// gradient of \f$1\,atm/cm\f$ acting across an area of
        /// \f$1\,cm^2\f$.
        ///
        namespace perm_details {
            const double p_grad   = atm / (prefix::centi*meter);
            const double area     = square(prefix::centi*meter);
            const double flux     = cubic (prefix::centi*meter) / second;
            const double velocity = flux / area;
            const double visc     = prefix::centi*Poise;
            const double darcy    = (velocity * visc) / p_grad;
            //                    == 1e-7 [m^2] / 101325
            //                    == 9.869232667160130e-13 [m^2]
        }
        const double darcy = perm_details::darcy;
        /// @}

        // Unit conversion support.
        //
        // Note: Under the penalty of treason will you be
        //
        //    using namespace Opm::unit::convert;
        //
        // I mean it!
        //
        namespace convert {
            // Convert from external units of measurements to equivalent
            // internal units of measurements.  Note: The internal units
            // of measurements are *ALWAYS*, and exclusively, SI.
            //
            // Example: Convert a double kx, containing a permeability
            // value in units of milli-darcy (mD) to the equivalent
            // value in SI units (m^2).
            //
            //    using namespace Opm::unit;
            //    using namespace Opm::prefix;
            //    convert::from(kx, milli*darcy);
            //
            inline double from(const double q, const double unit)
            {
                return q * unit;
            }

            // Convert from internal units of measurements to equivalent
            // external units of measurements.  Note: The internal units
            // of measurements are *ALWAYS*, and exclusively, SI.
            //
            // Example: Convert a std::vector<double> p, containing
            // pressure values in the SI unit Pascal (i.e., unit::Pascal)
            // to the equivalent values in Psi (unit::psia).
            //
            //    using namespace Opm::unit;
            //    std::transform(p.begin(), p.end(), p.begin(),
            //                   boost::bind(convert::to, _1, psia));
            //
            inline double to(const double q, const double unit)
            {
                return q / unit;
            }
        } // namespace convert
    } // namespace unit

    namespace units {
        //     const double MILLIDARCY = 1.0;//9.86923e-16;
        //     const double VISCOSITY_UNIT = 1.0;//1e-3;
        //     const double DAYS2SECONDS = 1.0;//86400;
        const double MILLIDARCY = 9.86923e-16;
        const double VISCOSITY_UNIT = 1e-3;
        const double DAYS2SECONDS = 86400;
        const double FEET = 0.30479999798832;
        const double WELL_INDEX_UNIT = VISCOSITY_UNIT/(DAYS2SECONDS*1e5);
    } // namespace units
} // namespace Opm
#endif // OPENRS_UNITS_HEADER
