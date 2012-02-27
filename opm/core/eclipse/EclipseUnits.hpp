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

#ifndef OPM_ECLIPSEUNITS_HEADER_INCLUDED
#define OPM_ECLIPSEUNITS_HEADER_INCLUDED

namespace Opm
{
    struct EclipseUnits
    {
        double length;
        double time;
        double density;
	double polymer_density;
        double pressure;
        double compressibility;
        double viscosity;
        double permeability;
        double liqvol_s;
        double liqvol_r;
        double gasvol_s;
        double gasvol_r;
        double transmissibility;

        void setToOne()
        {
            length = 1.0;
            time = 1.0;
            density = 1.0;
	    polymer_density = 1.0;
            pressure = 1.0;
            compressibility = 1.0;
            viscosity = 1.0;
            permeability = 1.0;
            liqvol_s = 1.0;
            liqvol_r = 1.0;
            gasvol_s = 1.0;
            gasvol_r = 1.0;
            transmissibility = 1.0;
        }
    };


} // namespace Opm


#endif // OPM_ECLIPSEUNITS_HEADER_INCLUDED
