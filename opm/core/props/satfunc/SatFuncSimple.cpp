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

#include "config.h"
#include <opm/core/props/satfunc/SatFuncSimple.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/satfunc/SaturationPropsFromDeck.hpp>
#include <opm/core/grid.h>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/utility/buildUniformMonotoneTable.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{

 // SatFuncSimple.cpp can now be removed and the code below moved to SatFuncBase.cpp ...
 
    template<>
    void SatFuncBase<NonuniformTableLinear<double> >::initializeTableType(NonuniformTableLinear<double> & table,
                                                                           const std::vector<double>& arg,
                                                                           const std::vector<double>& value,
                                                                           const int samples)
    {
      table = NonuniformTableLinear<double>(arg, value);
    }

    template<>
    void SatFuncBase<UniformTableLinear<double> >::initializeTableType(UniformTableLinear<double> & table,
                                                                       const std::vector<double>& arg,
                                                                       const std::vector<double>& value,
                                                                       const int samples)
    {
      buildUniformMonotoneTable(arg, value,  samples, table);
    }

    double EPSTransforms::Transform::scaleSat(double s, double s_r, double s_cr, double s_max) const
    {
        if (doNotScale) {
            return s;
        } else if (!do_3pt) { // 2-pt
            if (s <= scr) {
                return s_cr;
            } else {
                return (s >= smax) ? s_max : s_cr + (s-scr)*slope1;
            }
        } else if (s <= sr) {
            return (s <= scr) ? s_cr : s_cr+(s-scr)*slope1;
        } else {
            return (s >= smax) ? s_max : s_r+(s-sr)*slope2;
        }
    }

    double EPSTransforms::Transform::scaleSatInv(double ss, double s_r, double s_cr, double s_max) const
    {
        if (doNotScale) {
            return ss;
        } else if (!do_3pt) { // 2-pt
            if (ss <= s_cr) {
                return scr;
            } else {
                return (ss >= s_max) ? smax : scr + (ss-s_cr)/slope1;
            }
        } else if (ss <= s_r) {
            return (ss <= s_cr) ? scr : scr+(ss-s_cr)/slope1;
        } else {
            return (ss >= s_max) ? smax : sr+(ss-s_r)/slope2;
        }
    }

    double EPSTransforms::Transform::scaleSatDeriv(double s, double s_r, double s_cr, double s_max) const
    {
        if (doNotScale) {
            return 1.0;
        } else if (!do_3pt) { // 2-pt
            if (s <= scr) {
                return 0.0;
            } else {
                return (s >= smax) ? 0.0 : slope1;
            }
        } else if (s <= sr) {
            return (s <= scr) ? 0.0 : slope1;
        } else {
            return (s >= smax) ? 0.0 : slope2;
        }
    }


    double EPSTransforms::Transform::scaleSatPc(double s, double s_min, double s_max) const
    {
        if (doNotScale) {
            return s;
        } else if (s<=smin) {
            return s_min;
        } else if (s <= smax) {
            return s_min + (s-smin)*(s_max-s_min)/(smax-smin);
        } else {
            return s_max;
        }
    }
    double EPSTransforms::Transform::scaleSatDerivPc(double s, double s_min, double s_max) const
    {
        if (doNotScale) {
            return 1.0;
        } else if (s<smin) {
            return 0.0;
        } else if (s <= smax) {
            return (s_max-s_min)/(smax-smin);
        } else {
            return 0.0;
        }
    }

    double EPSTransforms::Transform::scaleKr(double s, double kr, double krsr_tab) const
    {
        if (doKrCrit) {
            if (s <= scr) {
                return 0.0;
            } else if (s <= sr) {
                return kr*krSlopeCrit;
            } else if (s <= smax) {
                if (doSatInterp)
                    return krsr + (s-sr)*krSlopeMax; // Note: Scaling independent of kr-value ...
                else
                    return krsr + (kr-krsr_tab)*krSlopeMax;
            } else {
                return krmax;
            }
        } else if (doKrMax) {
            if (s <= scr) {
                return 0.0;
            } else if (s <= smax) {
                return kr*krSlopeMax;
            } else {
                return krmax;
            }
        } else {
            return kr;
        }
    }


    double EPSTransforms::Transform::scaleKrDeriv(double s, double krDeriv)  const
    {
        if (doKrCrit) {
            if (s <= scr) {
                return 0.0;
            } else if (s <= sr) {
                return krDeriv*krSlopeCrit;
            } else if (s <= smax) {
                if (doSatInterp)
                    return krSlopeMax; // Note: Scaling independent of kr-value ...
                else
                    return krDeriv*krSlopeMax;
            } else {
                return 0.0;
            }
        } else if (doKrMax) {
            if (s <= scr) {
                return 0.0;
            } else if (s <= smax) {
                return krDeriv*krSlopeMax;
            } else {
                return 0.0;
            }
        } else {
            if (s <= scr) {
                return 0.0;
            } else if (s <= smax) {
                return krDeriv;
            } else {
                return 0.0;
            }
        }
   }
   
   
   void EPSTransforms::Transform::printMe(std::ostream & out)
   {
   
       out << "doNotScale: " << doNotScale << std::endl;
       out << "do_3pt: " << do_3pt << std::endl;
       out << "smin: " << smin << std::endl;
       out << "scr: " << scr << std::endl;
       out << "sr: " << sr << std::endl;
       out << "smax: " << smax << std::endl;
       out << "slope1: " << slope1 << std::endl;
       out << "slope2: " << slope2 << std::endl;
       out << "doKrMax: " << doKrMax << std::endl;
       out << "doKrCrit: " << doKrCrit << std::endl;
       out << "doSatInterp: " << doSatInterp << std::endl;
       out << "krsr: " << krsr << std::endl;
       out << "krmax: " << krmax << std::endl;
       out << "krSlopeMax: " << krSlopeMax << std::endl;
       out << "krSlopeCrit: " << krSlopeCrit << std::endl;
    }
    
   void SatHyst::printMe(std::ostream & out) 
   {
        out << "sg_hyst: " << sg_hyst << std::endl;
        out << "sg_shift: " << sg_shift << std::endl;
        out << "sow_hyst: " << sow_hyst << std::endl;
        out << "sow_shift: " << sow_shift << std::endl;
    };
        
   
} // namespace Opm
