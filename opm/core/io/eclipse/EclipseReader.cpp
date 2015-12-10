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
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/io/eclipse/EclipseIOUtil.hpp>

#include <algorithm>

#include <ert/ecl/ecl_file.h>

namespace Opm
{

    void restoreTemperatureData(const ecl_file_type* file,
                              EclipseStateConstPtr eclipse_state,
                              int numcells,
                              SimulatorState& simulator_state) {
        const char* temperature = "TEMP";

        if (ecl_file_has_kw(file , temperature)) {
            ecl_kw_type* temperature_kw = ecl_file_iget_named_kw(file, temperature, 0);

            if (ecl_kw_get_size(temperature_kw) != numcells) {
                throw std::runtime_error("Read of restart file: Could not restore temperature data, length of data from file not equal number of cells");
            }

            float* temperature_data = ecl_kw_get_float_ptr(temperature_kw);

            // factor and offset from the temperature values given in the deck to Kelvin
            double scaling = eclipse_state->getDeckUnitSystem()->parse("Temperature")->getSIScaling();
            double offset  = eclipse_state->getDeckUnitSystem()->parse("Temperature")->getSIOffset();

            for (size_t index = 0; index < simulator_state.temperature().size(); ++index) {
                simulator_state.temperature()[index] = unit::convert::from((double)temperature_data[index] - offset, scaling);
            }
          } else {
              throw std::runtime_error("Read of restart file: File does not contain TEMP data\n");
          }
    }


    void restorePressureData(const ecl_file_type* file,
                             EclipseStateConstPtr eclipse_state,
                             int numcells,
                             SimulatorState& simulator_state) {
        const char* pressure = "PRESSURE";

        if (ecl_file_has_kw(file , pressure)) {

            ecl_kw_type* pressure_kw = ecl_file_iget_named_kw(file, pressure, 0);

            if (ecl_kw_get_size(pressure_kw) != numcells) {
                throw std::runtime_error("Read of restart file: Could not restore pressure data, length of data from file not equal number of cells");
            }

            float* pressure_data = ecl_kw_get_float_ptr(pressure_kw);
            const double deck_pressure_unit = (eclipse_state->getDeckUnitSystem()->getType() == UnitSystem::UNIT_TYPE_METRIC) ? Opm::unit::barsa : Opm::unit::psia;
            for (size_t index = 0; index < simulator_state.pressure().size(); ++index) {
                simulator_state.pressure()[index] = unit::convert::from((double)pressure_data[index], deck_pressure_unit);
            }
        } else {
            throw std::runtime_error("Read of restart file: File does not contain PRESSURE data\n");
        }
    }


    void restoreSaturation(const ecl_file_type* file_type,
                           const PhaseUsage& phaseUsage,
                           int numcells,
                           SimulatorState& simulator_state) {

        float* sgas_data = NULL;
        float* swat_data = NULL;

        if (phaseUsage.phase_used[BlackoilPhases::Aqua]) {
            const char* swat = "SWAT";
            if (ecl_file_has_kw(file_type, swat)) {
                ecl_kw_type* swat_kw = ecl_file_iget_named_kw(file_type , swat, 0);
                swat_data = ecl_kw_get_float_ptr(swat_kw);
                std::vector<double> swat_datavec(&swat_data[0], &swat_data[numcells]);
                EclipseIOUtil::addToStripedData(swat_datavec, simulator_state.saturation(), phaseUsage.phase_pos[BlackoilPhases::Aqua], phaseUsage.num_phases);
            } else {
                std::string error_str = "Restart file is missing SWAT data!\n";
                throw std::runtime_error(error_str);
            }
        }

        if (phaseUsage.phase_used[BlackoilPhases::Vapour]) {
            const char* sgas = "SGAS";
            if (ecl_file_has_kw(file_type, sgas)) {
                ecl_kw_type* sgas_kw = ecl_file_iget_named_kw(file_type , sgas, 0);
                sgas_data = ecl_kw_get_float_ptr(sgas_kw);
                std::vector<double> sgas_datavec(&sgas_data[0], &sgas_data[numcells]);
                EclipseIOUtil::addToStripedData(sgas_datavec, simulator_state.saturation(), phaseUsage.phase_pos[BlackoilPhases::Vapour], phaseUsage.num_phases);
            } else {
                std::string error_str = "Restart file is missing SGAS data!\n";
                throw std::runtime_error(error_str);
            }
        }
    }


    void restoreRSandRV(const ecl_file_type* file_type,
                        SimulationConfigConstPtr sim_config,
                        int numcells,
                        BlackoilState* blackoil_state) {

        if (sim_config->hasDISGAS()) {
            const char* RS = "RS";
            if (ecl_file_has_kw(file_type, RS)) {
                ecl_kw_type* rs_kw = ecl_file_iget_named_kw(file_type, RS, 0);
                float* rs_data = ecl_kw_get_float_ptr(rs_kw);
                std::vector<double> rs_datavec(&rs_data[0], &rs_data[numcells]);
                blackoil_state->gasoilratio().clear();
                blackoil_state->gasoilratio().insert(blackoil_state->gasoilratio().begin(), rs_datavec.begin(), rs_datavec.end());
            } else {
                throw std::runtime_error("Restart file is missing RS data!\n");
            }
        }

        if (sim_config->hasVAPOIL()) {
            const char* RV = "RV";
            if (ecl_file_has_kw(file_type, RV)) {
                ecl_kw_type* rv_kw = ecl_file_iget_named_kw(file_type, RV, 0);
                float* rv_data = ecl_kw_get_float_ptr(rv_kw);
                std::vector<double> rv_datavec(&rv_data[0], &rv_data[numcells]);
                blackoil_state->rv().clear();
                blackoil_state->rv().insert(blackoil_state->rv().begin(), rv_datavec.begin(), rv_datavec.end());
            } else {
                throw std::runtime_error("Restart file is missing RV data!\n");
            }
        }
    }


    void restoreSOLUTION(const std::string& restart_filename,
                         int reportstep,
                         bool unified,
                         EclipseStateConstPtr eclipseState,
                         int numcells,
                         const PhaseUsage& phaseUsage,
                         SimulatorState& simulator_state)
    {
        const char* filename = restart_filename.c_str();
        ecl_file_type* file_type = ecl_file_open(filename, 0);

        if (file_type) {
            bool block_selected = unified ? ecl_file_select_rstblock_report_step(file_type , reportstep) : true;

            if (block_selected) {
                restorePressureData(file_type, eclipseState, numcells, simulator_state);
                restoreTemperatureData(file_type, eclipseState, numcells, simulator_state);
                restoreSaturation(file_type, phaseUsage, numcells, simulator_state);
                BlackoilState* blackoilState = dynamic_cast<BlackoilState*>(&simulator_state);
                if (blackoilState) {
                    SimulationConfigConstPtr sim_config = eclipseState->getSimulationConfig();
                    restoreRSandRV(file_type, sim_config, numcells, blackoilState);
                }
            } else {
                std::string error_str = "Restart file " +  restart_filename + " does not contain data for report step " + std::to_string(reportstep) + "!\n";
                throw std::runtime_error(error_str);
            }
            ecl_file_close(file_type);
        } else {
            std::string error_str = "Restart file " + restart_filename + " not found!\n";
            throw std::runtime_error(error_str);
        }
    }


    void restoreOPM_XWELKeyword(const std::string& restart_filename, int reportstep, bool unified, WellState& wellstate)
    {
        const char * keyword = "OPM_XWEL";
        const char* filename = restart_filename.c_str();
        ecl_file_type* file_type = ecl_file_open(filename, 0);

        if (file_type != NULL) {

            bool block_selected = unified ? ecl_file_select_rstblock_report_step(file_type , reportstep) : true;

            if (block_selected) {
                ecl_kw_type* xwel = ecl_file_iget_named_kw(file_type , keyword, 0);
                const double* xwel_data = ecl_kw_get_double_ptr(xwel);
                std::copy_n(xwel_data + wellstate.getRestartTemperatureOffset(), wellstate.temperature().size(), wellstate.temperature().begin());
                std::copy_n(xwel_data + wellstate.getRestartBhpOffset(), wellstate.bhp().size(), wellstate.bhp().begin());
                std::copy_n(xwel_data + wellstate.getRestartPerfPressOffset(), wellstate.perfPress().size(), wellstate.perfPress().begin());
                std::copy_n(xwel_data + wellstate.getRestartPerfRatesOffset(), wellstate.perfRates().size(), wellstate.perfRates().begin());
                std::copy_n(xwel_data + wellstate.getRestartWellRatesOffset(), wellstate.wellRates().size(), wellstate.wellRates().begin());
            } else {
                std::string error_str = "Restart file " +  restart_filename + " does not contain data for report step " + std::to_string(reportstep) + "!\n";
                throw std::runtime_error(error_str);
            }
            ecl_file_close(file_type);
        } else {
            std::string error_str = "Restart file " + restart_filename + " not found!\n";
            throw std::runtime_error(error_str);
        }
    }



    void init_from_restart_file(EclipseStateConstPtr eclipse_state,
                                int numcells,
                                const PhaseUsage& phase_usage,
                                SimulatorState& simulator_state,
                                WellState& wellstate) {

        InitConfigConstPtr initConfig        = eclipse_state->getInitConfig();
        IOConfigConstPtr ioConfig            = eclipse_state->getIOConfig();
        int restart_step                     = initConfig->getRestartStep();
        const std::string& restart_file_root = initConfig->getRestartRootName();
        bool output                          = false;
        const std::string& restart_file_name = ioConfig->getRestartFileName(restart_file_root, restart_step, output);

        Opm::restoreSOLUTION(restart_file_name, restart_step, ioConfig->getUNIFIN(), eclipse_state, numcells, phase_usage, simulator_state);
        Opm::restoreOPM_XWELKeyword(restart_file_name, restart_step, ioConfig->getUNIFIN(), wellstate);
    }


} // namespace Opm
