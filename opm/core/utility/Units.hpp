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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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
    {
        // Common powers
        inline double square(double v) { return v * v;     }
        inline double cubic (double v) { return v * v * v; }

        // --------------------------------------------------------------
        // Basic (fundamental) units and conversions
        // --------------------------------------------------------------

        // Length:
        const double meter =  1;
        const double inch  =  2.54 * prefix::centi*meter;
        const double feet  = 12    * inch;

        // Time:
        const double second =   1;
        const double minute =  60 * second;
        const double hour   =  60 * minute;
        const double day    =  24 * hour;
        const double year   = 365 * day;

        // Volume
        const double stb    = 0.158987294928 * cubic(meter);


        // Mass:
        const double kilogram = 1;

        // http://en.wikipedia.org/wiki/Pound_(mass)#Avoirdupois_pound
        const double pound    = 0.45359237 * kilogram;

        // --------------------------------------------------------------
        // Standardised constants
        // --------------------------------------------------------------

        const double gravity = 9.80665 * meter/square(second);

        // --------------------------------------------------------------
        // Derived units and conversions
        // --------------------------------------------------------------

        // Force:
        const double Newton = kilogram*meter / square(second); // == 1
        const double lbf    = pound * gravity; // Pound-force

        // Pressure:
        const double Pascal = Newton / square(meter); // == 1
        const double barsa  = 100000 * Pascal;
        const double atm    = 101325 * Pascal;
        const double psia   = lbf / square(inch);

        // Viscosity:
        const double Pas   = Pascal * second; // == 1
        const double Poise = prefix::deci*Pas;

        // Permeability:
        //
        // A porous medium with a permeability of 1 darcy permits a
        // flow (flux) of 1 cm³/s of a fluid with viscosity 1 cP (1
        // mPa·s) under a pressure gradient of 1 atm/cm acting across
        // an area of 1 cm².
        //
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
