 /*
  Copyright 2013 Andreas Lauser
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#include "BlackoilEclipseOutputWriter.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/DataMap.hpp>

#ifdef HAVE_ERT
#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_rst_file.h>
#endif

namespace Opm {
void BlackoilEclipseOutputWriter::writeInitFile(const SimulatorTimer &timer)
{
    tm tm = boost::posix_time::to_tm(timer.currentDateTime());
    startTime_ = mktime(&tm);

#if HAVE_ERT
    writeGridInitFile_(timer);
    writeSummaryHeaderFile_(timer);
#else
    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
#endif // HAVE_ERT
}

void BlackoilEclipseOutputWriter::writeReservoirState(const BlackoilState& reservoirState, const SimulatorTimer& timer)
{
#if HAVE_ERT
    ecl_file_enum file_type = ECL_UNIFIED_RESTART_FILE;  // Alternatively ECL_RESTART_FILE for multiple restart files.
    bool fmt_file           = false;

    char *fileName = ecl_util_alloc_filename(outputDir_.c_str(),
                                             baseName_.c_str(),
                                             /*file_type=*/ECL_UNIFIED_RESTART_FILE,
                                             fmt_file,
                                             timer.currentStepNum());
    int phases = ECL_OIL_PHASE + ECL_GAS_PHASE + ECL_WATER_PHASE;
    double days = Opm::unit::convert::to(timer.currentTime(), Opm::unit::day);
    int nx = grid_.cartdims[0];
    int ny = grid_.cartdims[1];
    int nz = grid_.cartdims[2];
    int nactive = grid_.number_of_cells;
    ecl_rst_file_type* rst_file;

    time_t curTime;
    tm tm = boost::posix_time::to_tm(timer.currentDateTime());
    curTime = mktime(&tm);

    if (timer.currentStepNum() > 0 && file_type == ECL_UNIFIED_RESTART_FILE)
        rst_file = ecl_rst_file_open_append(fileName);
    else
        rst_file = ecl_rst_file_open_write(fileName);

    ecl_rst_file_fwrite_header(rst_file, timer.currentStepNum(), curTime, days, nx, ny, nz, nactive, phases);
    ecl_rst_file_start_solution(rst_file);

    {
        ecl_kw_type* pressure_kw = newEclDoubleKeyword_("PRESSURE", reservoirState.pressure());
        ecl_rst_file_add_kw(rst_file, pressure_kw);
        ecl_kw_free(pressure_kw);
    }

    {
        ecl_kw_type* swat_kw = newEclDoubleKeyword_("SWAT", reservoirState.saturation(), /*offset=*/0, /*stride=*/3);
        ecl_rst_file_add_kw(rst_file, swat_kw);
        ecl_kw_free(swat_kw);
    }

    {
        ecl_kw_type* soil_kw = newEclDoubleKeyword_("SOIL", reservoirState.saturation(), /*offset=*/1, /*stride=*/3);
        ecl_rst_file_add_kw(rst_file, soil_kw);
        ecl_kw_free(soil_kw);
    }

    {
        ecl_kw_type* sgas_kw = newEclDoubleKeyword_("SGAS", reservoirState.saturation(), /*offset=*/2, /*stride=*/3);
        ecl_rst_file_add_kw(rst_file, sgas_kw);
        ecl_kw_free(sgas_kw);
    }

    ecl_rst_file_end_solution(rst_file);
    ecl_rst_file_close(rst_file);
    free(fileName);
#else
    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
#endif // HAVE_ERT
}

void BlackoilEclipseOutputWriter::writeWellState(const WellState& wellState, const SimulatorTimer& timer)
{
#if HAVE_ERT
    tm tm = boost::posix_time::to_tm(timer.currentDateTime());
    time_t curTime = mktime(&tm);

    // create a new timestep for the summary file (at least if the
    // timer was advanced since the last call to writeWellState())
    ecl_sum_tstep_type* tstep=
        ecl_sum_add_tstep(sumWriter_,
                          timer.currentStepNum() + 1,
                          (curTime - startTime_)/(24*60*60));

    int numWells = accumulatedProducedFluids_.size();
    for (int wellIdx = 0; wellIdx < numWells; ++wellIdx) {
        // set the value for the well oil production rate
        double woprValue = wellState.wellRates()[wellIdx*3 + BlackoilPhases::Liquid];
        woprValue *= - 1 * (24 * 60 * 60); // convert kg^3/s of injected fluid to kg^3/d of produced fluid
        woprValue /= 700; // convert from kg/d to m^3/d. TODO: ignore dissolved gas, use real surface density!
        woprValue = std::max(0.0, woprValue);

        int woprIdx = smspec_node_get_params_index(woprSmspec_[wellIdx]);
        ecl_sum_tstep_iset(tstep, woprIdx, woprValue);

        double wwirValue = wellState.wellRates()[wellIdx*3 + BlackoilPhases::Aqua];
        wwirValue *= 1 * (24 * 60 * 60);  // convert kg^3/s to kg^3/d
        wwirValue /= 700; // convert from kg/d to m^3/d. TODO: use real surface density!
        wwirValue = std::max(0.0, wwirValue);

        int wwirIdx = smspec_node_get_params_index(wwirSmspec_[wellIdx]);
        ecl_sum_tstep_iset(tstep, wwirIdx, wwirValue);

        // accumulate injected produced fluids
        for (int phaseIdx = 0; phaseIdx < /*numPhases=*/3; ++phaseIdx) {
            // convert from m^3 of injected fluid to m^3 of produced fluid
            double injectedMass = wellState.wellRates()[wellIdx*3 + phaseIdx];
            injectedMass *= timer.currentStepLength();

            // TODO: use real surface density!
            double rho;
            switch (phaseIdx) {
            case BlackoilPhases::Liquid:
                rho = 700;
                break;
            case BlackoilPhases::Aqua:
                rho = 700;
                break;
            case BlackoilPhases::Vapour:
                rho = 1;
                break;
            }
            double injectedVolume = injectedMass/rho; // convert from kg/d to m^3/d. TODO: ignore dissolved gas

            if (injectedMass < 0)
                accumulatedProducedFluids_[wellIdx][phaseIdx] += -injectedVolume;
            else
                accumulatedInjectedFluids_[wellIdx][phaseIdx] += injectedVolume;

            int woptIdx = smspec_node_get_params_index(woptSmspec_[wellIdx]);
            ecl_sum_tstep_iset(tstep, woptIdx, accumulatedProducedFluids_[wellIdx][BlackoilPhases::Liquid]);

            int wwitIdx = smspec_node_get_params_index(wwitSmspec_[wellIdx]);
            ecl_sum_tstep_iset(tstep, wwitIdx, accumulatedInjectedFluids_[wellIdx][BlackoilPhases::Aqua]);
        }
    }

    ecl_sum_fwrite(sumWriter_);
#else
    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
#endif // HAVE_ERT
}

#if HAVE_ERT
void BlackoilEclipseOutputWriter::writeGridInitFile_(const SimulatorTimer &timer)
{
    int phases = ECL_OIL_PHASE + ECL_GAS_PHASE + ECL_WATER_PHASE;
    bool endian_flip  = true;//ECL_ENDIAN_FLIP;
    bool fmt_file = false;
    ecl_file_enum file_type = ECL_EGRID_FILE;

    ecl_grid_type* ecl_grid = newEclGrid_();
    char* gridFileName = ecl_util_alloc_filename(outputDir_.c_str(), baseName_.c_str(), file_type, fmt_file, timer.currentStepNum());
    fortio_type* fortio;
    ecl_grid_fwrite_EGRID(ecl_grid, gridFileName);
    free(gridFileName);

    char* initFileName = ecl_util_alloc_filename(outputDir_.c_str(), baseName_.c_str(), /*file_type=*/ECL_INIT_FILE, fmt_file, timer.currentStepNum());
    if (!ecl_util_fmt_file(initFileName, &fmt_file)) {
        OPM_THROW(std::runtime_error,
                  "Could not determine formatted/unformatted status of file:" << initFileName << " non-standard name?" << std::endl);
    }
    fortio = fortio_open_writer(initFileName, fmt_file, endian_flip);
    {
        time_t start_date;

        {
            boost::posix_time::ptime start_date_(timer.currentDateTime());

            tm td_tm = boost::posix_time::to_tm(start_date_);
            start_date = mktime(&td_tm);
        }

        ecl_kw_type* poro_kw = newEclDoubleKeyword_(PORO_KW, eclipseParser_.getFloatingPointValue("PORO"));
        ecl_init_file_fwrite_header(fortio, ecl_grid, poro_kw, phases, start_date);
        ecl_kw_free(poro_kw);
    }

    /* This collection of keywords is somewhat arbitrary and random. */
    saveEclKeyword_(fortio, "PERMX", ECL_FLOAT_TYPE);
    saveEclKeyword_(fortio, "PERMY", ECL_FLOAT_TYPE);
    saveEclKeyword_(fortio, "PERMZ", ECL_FLOAT_TYPE);

    fortio_fclose(fortio);
    free(initFileName);
}

void BlackoilEclipseOutputWriter::writeSummaryHeaderFile_(const SimulatorTimer &timer)
{
    std::string caseName;
    if (!outputDir_.empty())
        caseName += outputDir_ + "/";
    caseName += baseName_;

    if (sumWriter_)
        ecl_sum_free(sumWriter_);

    // allocate the data structure for the writer
    sumWriter_ =
        ecl_sum_alloc_writer(caseName.c_str(),
                             /*formattedOutput=*/true,
                             /*unifiedOutput=*/true,
                             /*joinString=*/":",
                             startTime_,
                             grid_.cartdims[0],grid_.cartdims[1],grid_.cartdims[2]);

    // initialize the accumulated masses to zero
    const auto &wellSpecs = eclipseParser_.getWELSPECS().welspecs;
    int numWells = wellSpecs.size();
    accumulatedProducedFluids_.resize(numWells);
    accumulatedInjectedFluids_.resize(numWells);
    for (int wellIdx = 0; wellIdx < numWells; ++wellIdx) {
        for (int phaseIdx = 0; phaseIdx < /*numPhases=*/3; ++phaseIdx) {
            accumulatedProducedFluids_[wellIdx][phaseIdx] = 0;
            accumulatedInjectedFluids_[wellIdx][phaseIdx] = 0;
        }
    }

    woprSmspec_.resize(numWells);
    woptSmspec_.resize(numWells);
    wwirSmspec_.resize(numWells);
    wwitSmspec_.resize(numWells);
    auto wellIt = wellSpecs.begin();
    const auto &wellEndIt = wellSpecs.end();
    for (int wellIdx = 0; wellIt != wellEndIt; ++wellIt, ++wellIdx) {
        // add the variables which ought to be included in the summary
        // file
        woprSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                  /*varName=*/"WOPR",
                                  /*wellGroupName=*/wellIt->name_.c_str(),
                                  /*num=*/0,
                                  /*unit=*/"SM3/DAY",
                                  /*defaultValue=*/0.0);

        woptSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                               /*varName=*/"WOPT",
                                               /*wellGroupName=*/wellIt->name_.c_str(),
                                               /*num=*/0,
                                               /*unit=*/"SM3",
                                               /*defaultValue=*/0.0);

        wwirSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                  /*varName=*/"WWIR",
                                  /*wellGroupName=*/wellIt->name_.c_str(),
                                  /*num=*/0,
                                  /*unit=*/"SM3/DAY",
                                  /*defaultValue=*/0.0);

        wwitSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                               /*varName=*/"WWIT",
                                               /*wellGroupName=*/wellIt->name_.c_str(),
                                               /*num=*/0,
                                               /*unit=*/"SM3",
                                               /*defaultValue=*/0.0);
    }

    ecl_sum_fwrite(sumWriter_);
}

ecl_grid_type* BlackoilEclipseOutputWriter::newEclGrid_()
{
    if (eclipseParser_.hasField("DXV")) {
        // make sure that the DYV and DZV keywords are present if the
        // DXV keyword is used in the deck...
        assert(eclipseParser_.hasField("DYV"));
        assert(eclipseParser_.hasField("DZV"));

        const auto &dxv = eclipseParser_.getFloatingPointValue("DXV");
        const auto &dyv = eclipseParser_.getFloatingPointValue("DYV");
        const auto &dzv = eclipseParser_.getFloatingPointValue("DZV");

        // creating a C array out of std::vector like this is pretty
        // hacky and might even be unportable. having said that, it
        // probably works with all currently known STL
        // implementations...
        return ecl_grid_alloc_dxv_dyv_dzv(dxv.size(), dyv.size(), dzv.size(),
                                          &dxv[0], &dyv[0], &dzv[0],
                                          /*actnum=*/NULL);
    }
    if (eclipseParser_.hasField("ZCORN")) {
        struct grdecl grdecl = eclipseParser_.get_grdecl();

        ecl_kw_type * coord_kw   = newEclDoubleKeyword_(COORD_KW, eclipseParser_.getFloatingPointValue("COORD"));
        ecl_kw_type * zcorn_kw   = newEclDoubleKeyword_(ZCORN_KW, eclipseParser_.getFloatingPointValue("ZCORN"));
        ecl_kw_type * actnum_kw  = newEclIntKeyword_(ACTNUM_KW, eclipseParser_.getIntegerValue("ACTNUM"));
        ecl_kw_type * mapaxes_kw = NULL;

        ecl_grid_type * grid ;
        if (grdecl.mapaxes != NULL)
            mapaxes_kw = newEclDoubleKeyword_(MAPAXES_KW, eclipseParser_.getFloatingPointValue("MAPAXES"));

        grid = ecl_grid_alloc_GRDECL_kw(grdecl.dims[0], grdecl.dims[1], grdecl.dims[2], zcorn_kw, coord_kw, actnum_kw, mapaxes_kw);

        ecl_kw_free(coord_kw);
        ecl_kw_free(zcorn_kw);
        ecl_kw_free(actnum_kw);
        if (mapaxes_kw != NULL)
            ecl_kw_free(mapaxes_kw);

        return grid;
    }
    OPM_THROW(std::runtime_error,
              "Can't create an ERT grid (no supported keywords found in deck)");
}

ecl_kw_type* BlackoilEclipseOutputWriter::newEclIntKeyword_(const std::string& kwName,
                                                            const std::vector<int> &data,
                                                            int offset,
                                                            int stride)
{
    assert(offset >= 0 && offset < data.size());
    assert(stride > 0 && stride < data.size() - offset);

    ecl_kw_type* eclKw =
        ecl_kw_alloc(kwName.c_str(),
                     grid_.number_of_cells, ECL_INT_TYPE);
    for (int i=0; i < grid_.number_of_cells; i++)
        ecl_kw_iset_float(eclKw, i, data[i*stride + offset]);
    return eclKw;
}

ecl_kw_type* BlackoilEclipseOutputWriter::newEclDoubleKeyword_(const std::string& kwName,
                                                               const std::vector<double> &data,
                                                               int offset,
                                                               int stride)
{
    assert(offset >= 0 && offset < data.size());
    assert(stride > 0 && stride < data.size() - offset);

    ecl_kw_type* eclKw =
        ecl_kw_alloc(kwName.c_str(),
                     grid_.number_of_cells, ECL_FLOAT_TYPE);
    for (int i=0; i < grid_.number_of_cells; i++)
        ecl_kw_iset_float(eclKw, i, static_cast<float>(data[i*stride + offset]));
    return eclKw;
}


void BlackoilEclipseOutputWriter::saveEclKeyword_(fortio_type* fortio, const std::string& kw, ecl_type_enum eclType)
{
    ecl_kw_type* eclKw;
    if (eclType == ECL_INT_TYPE)
        eclKw = newEclIntKeyword_(kw, eclipseParser_.getIntegerValue(kw));
    else if (eclType == ECL_FLOAT_TYPE)
        eclKw = newEclDoubleKeyword_(kw, eclipseParser_.getFloatingPointValue(kw));
    else
        OPM_THROW(std::logic_error,
                  "Not implemented: ECL keywords of type " << ECL_FLOAT_TYPE);
    if (eclKw != NULL) {
        ecl_kw_fwrite(eclKw, fortio);
        ecl_kw_free(eclKw);
    }
}

#endif // HAVE_ERT

} // namespace Opm
