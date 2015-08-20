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
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/io/eclipse/EclipseIOUtil.hpp>

#include <algorithm>

#include <ert/ecl/ecl_file.h>

namespace Opm
{
    void restoreOPM_XWELKeyword(const std::string& restart_filename, int reportstep, WellState& wellstate)
    {
        const char * keyword = "OPM_XWEL";
        const char* filename = restart_filename.c_str();
        ecl_file_type* file_type = ecl_file_open(filename, 0);

        if (file_type != NULL) {
            bool block_selected = ecl_file_select_rstblock_report_step(file_type , reportstep);

            if (block_selected) {
                ecl_kw_type* xwel = ecl_file_iget_named_kw(file_type , keyword, 0);
                const double* xwel_data = ecl_kw_get_double_ptr(xwel);

                std::copy_n(xwel_data + wellstate.getRestartTemperatureOffset(), wellstate.temperature().size(), wellstate.temperature().begin());
                std::copy_n(xwel_data + wellstate.getRestartBhpOffset(), wellstate.bhp().size(), wellstate.bhp().begin());
                std::copy_n(xwel_data + wellstate.getRestartPerfPressOffset(), wellstate.perfPress().size(), wellstate.perfPress().begin());
                std::copy_n(xwel_data + wellstate.getRestartPerfRatesOffset(), wellstate.perfRates().size(), wellstate.perfRates().begin());
                std::copy_n(xwel_data + wellstate.getRestartWellRatesOffset(), wellstate.wellRates().size(), wellstate.wellRates().begin());
            }

            ecl_file_close(file_type);
        }
    }



namespace {
    void restoreTemperatureData(const ecl_file_type* file,
                                const EclipseState& eclipse_state,
                                const UnstructuredGrid& grid,
                                SimulatorState& simulator_state) {
        const char* temperature_kw_string = "TEMP";
        ecl_kw_type* temperature_kw       = ecl_file_iget_named_kw(file, temperature_kw_string, 0);

        if (ecl_kw_get_size(temperature_kw) != Opm::UgGridHelpers::numCells(grid)) {
            throw std::runtime_error("Read of restart file: Could not restore temperature data, length of data from file not equal number of cells");
        }

        float* temperature_data = ecl_kw_get_float_ptr(temperature_kw);

        // factor and offset from the temperature values given in the deck to Kelvin
        double scaling = eclipse_state.getDeckUnitSystem()->parse("Temperature")->getSIScaling();
        double offset  = eclipse_state.getDeckUnitSystem()->parse("Temperature")->getSIOffset();

        for (size_t index = 0; index < simulator_state.temperature().size(); ++index) {
            simulator_state.temperature()[index] = unit::convert::from((double)temperature_data[index] - offset, scaling);
        }
    }


    void restorePressureData(const ecl_file_type* file,
                             const EclipseState& eclipse_state,
                             const UnstructuredGrid& grid,
                             SimulatorState& simulator_state) {
        const char* pressure_keyword    = "PRESSURE";
        ecl_kw_type* pressure_kw        = ecl_file_iget_named_kw(file, pressure_keyword, 0);

        if (ecl_kw_get_size(pressure_kw) != Opm::UgGridHelpers::numCells(grid)) {
            throw std::runtime_error("Read of restart file: Could not restore pressure data, length of data from file not equal number of cells");
        }

        float* pressure_data = ecl_kw_get_float_ptr(pressure_kw);
        const double deck_pressure_unit = (eclipse_state.getDeckUnitSystem()->getType() == UnitSystem::UNIT_TYPE_METRIC) ? Opm::unit::barsa : Opm::unit::psia;
        for (size_t index = 0; index < simulator_state.pressure().size(); ++index) {
            simulator_state.pressure()[index] = unit::convert::from((double)pressure_data[index], deck_pressure_unit);
        }
    }
}


    void restoreSOLUTIONData(const std::string& restart_filename,
                             int reportstep,
                             const EclipseState& eclipseState,
                             const UnstructuredGrid& grid,
                             const PhaseUsage& phaseUsage,
                             SimulatorState& simulator_state)
    {
        const char* filename = restart_filename.c_str();
        ecl_file_type* file_type = ecl_file_open(filename, 0);

        if (file_type != NULL) {
            bool block_selected = ecl_file_select_rstblock_report_step(file_type , reportstep);

            if (block_selected) {

                restorePressureData(file_type, eclipseState, grid, simulator_state);
                restoreTemperatureData(file_type, eclipseState, grid, simulator_state);

                int numcells = Opm::UgGridHelpers::numCells(grid);

                float* sgas_data = NULL;
                float* swat_data = NULL;

                if (phaseUsage.phase_used[BlackoilPhases::Aqua]) {
                    ecl_kw_type* swat_kw = ecl_file_iget_named_kw(file_type , "SWAT", 0);
                    swat_data = ecl_kw_get_float_ptr(swat_kw);
                    std::vector<double> swat_datavec(&swat_data[0], &swat_data[numcells]);
                    EclipseIOUtil::addToStripedData(swat_datavec, simulator_state.saturation(), phaseUsage.phase_pos[BlackoilPhases::Aqua], phaseUsage.num_phases);
                }

                if (phaseUsage.phase_used[BlackoilPhases::Vapour]) {
                    ecl_kw_type* sgas_kw = ecl_file_iget_named_kw(file_type , "SGAS", 0);
                    sgas_data = ecl_kw_get_float_ptr(sgas_kw);
                    std::vector<double> sgas_datavec(&sgas_data[0], &sgas_data[numcells]);
                    EclipseIOUtil::addToStripedData(sgas_datavec, simulator_state.saturation(), phaseUsage.phase_pos[BlackoilPhases::Vapour], phaseUsage.num_phases);
                }
            }

            ecl_file_close(file_type);
        }
    }




} // namespace Opm
