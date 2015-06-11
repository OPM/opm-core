/*
  Copyright 2015 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM. If not, see <http://www.gnu.org/licenses/>.
*/

#include "EclipseReader.hpp"
#include <iostream>
#include <opm/core/simulator/WellState.hpp>

#include <ert/ecl/ecl_file.h>

namespace Opm
{
    void restoreOPM_XWELKeyword(const std::string restart_filename, int reportstep, WellState& wellstate)
    {
        const char * keyword = "OPM_XWEL";
        const char* filename = restart_filename.c_str();

        ecl_file_type* file_type = ecl_file_open(filename, 0);
        bool block_selected = ecl_file_select_rstblock_report_step(file_type , reportstep);

        if (block_selected) {
            ecl_kw_type* xwel = ecl_file_iget_named_kw(file_type , keyword, 0);
            const double* xwel_data = ecl_kw_get_double_ptr(xwel);
            size_t offset = 0;

            for (size_t i = 0; i < wellstate.bhp().size(); ++i) {
                wellstate.bhp()[i] = xwel_data[offset++];
            }

            for (size_t i = 0; i < wellstate.perfPress().size(); ++i) {
                wellstate.perfPress()[i] = xwel_data[offset++];
            }

            for (size_t i = 0; i < wellstate.perfRates().size(); ++i) {
                wellstate.perfRates()[i] = xwel_data[offset++];
            }

            for (size_t i = 0; i < wellstate.temperature().size(); ++i) {
                wellstate.temperature()[i] = xwel_data[offset++];
            }

            for (size_t i = 0; i < wellstate.wellRates().size(); ++i) {
                wellstate.wellRates()[i] = xwel_data[offset++];
            }
        }
    }
} // namespace Opm
