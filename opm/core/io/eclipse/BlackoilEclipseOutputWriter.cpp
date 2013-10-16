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

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/DataMap.hpp>

#ifdef HAVE_ERT
#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_rst_file.h>
#endif

namespace Opm {
void BlackoilEclipseOutputWriter::writeInitFile(const SimulatorTimer &timer)
{
#if HAVE_ERT
    int phases        = ECL_OIL_PHASE + ECL_WATER_PHASE;
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
        ecl_kw_type* poro_kw = newEclKeyword_(PORO_KW, ECL_FLOAT_TYPE);
        time_t start_date;

        {
            boost::posix_time::ptime start_date_(timer.currentDateTime());

            tm td_tm = boost::posix_time::to_tm(start_date_);
            start_date = mktime(&td_tm);
        }

        ecl_init_file_fwrite_header(fortio, ecl_grid, poro_kw, phases, start_date);
        ecl_kw_free(poro_kw);
    }

    /* This collection of keywords is somewhat arbitrary and random. */
    saveEclKeyword_(fortio, "PERMX", ECL_FLOAT_TYPE);
    saveEclKeyword_(fortio, "PERMY", ECL_FLOAT_TYPE);
    saveEclKeyword_(fortio, "PERMZ", ECL_FLOAT_TYPE);

    fortio_fclose(fortio);
    free(initFileName);
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
    std::cout << "writing: " << fileName << "\n";
    int phases = ECL_OIL_PHASE + ECL_GAS_PHASE + ECL_WATER_PHASE;
    double days = Opm::unit::convert::to(timer.currentTime(), Opm::unit::day);
    time_t date = 0;
    int nx = grid_.cartdims[0];
    int ny = grid_.cartdims[1];
    int nz = grid_.cartdims[2];
    int nactive = grid_.number_of_cells;
    ecl_rst_file_type* rst_file;

    {
        using namespace boost::posix_time;
        ptime t0(boost::gregorian::date(1970, 1, 1));
        time_duration::sec_type seconds = (timer.currentDateTime() - t0).total_seconds();

        date = time_t(seconds);
    }

    if (timer.currentStepNum() > 0 && file_type == ECL_UNIFIED_RESTART_FILE)
        rst_file = ecl_rst_file_open_append(fileName);
    else
        rst_file = ecl_rst_file_open_write(fileName);

    ecl_rst_file_fwrite_header(rst_file, timer.currentStepNum(), date, days, nx, ny, nz, nactive, phases);
    ecl_rst_file_start_solution(rst_file);

    {
        ecl_kw_type* pressure_kw = eclKeywordWrapper_("PRESSURE", reservoirState.pressure(),
                                                      /*offset=*/0, /*stride=*/1);
        ecl_rst_file_add_kw(rst_file, pressure_kw);
        ecl_kw_free(pressure_kw);
    }

    {
        ecl_kw_type* swat_kw = eclKeywordWrapper_("SWAT", reservoirState.saturation(), /*offset=*/0, /*stride=*/3);
        ecl_rst_file_add_kw(rst_file, swat_kw);
        ecl_kw_free(swat_kw);
    }

    {
        ecl_kw_type* soil_kw = eclKeywordWrapper_("SOIL", reservoirState.saturation(), /*offset=*/1, /*stride=*/3);
        ecl_rst_file_add_kw(rst_file, soil_kw);
        ecl_kw_free(soil_kw);
    }

    {
        ecl_kw_type* sgas_kw = eclKeywordWrapper_("SGAS", reservoirState.saturation(), /*offset=*/2, /*stride=*/3);
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

void BlackoilEclipseOutputWriter::writeWellState(const WellState& wellState)
{
    for (int i = 0; i < wellState.perfPress().size(); ++i)
        std::cout << "Perforation pressure " << i << ": " << wellState.perfPress()[i] << "\n";
    for (int i = 0; i < wellState.perfRates().size(); ++i)
        std::cout << "Well rate " << i << ": " << wellState.perfRates()[i] << "\n";
}

#if HAVE_ERT
ecl_grid_type* BlackoilEclipseOutputWriter::newEclGrid_() {
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

        ecl_kw_type * coord_kw   = newEclKeyword_(COORD_KW, ECL_FLOAT_TYPE);
        ecl_kw_type * zcorn_kw   = newEclKeyword_(ZCORN_KW, ECL_FLOAT_TYPE);
        ecl_kw_type * actnum_kw  = newEclKeyword_(ACTNUM_KW, ECL_INT_TYPE);
        ecl_kw_type * mapaxes_kw = NULL;

        ecl_grid_type * grid ;
        if (grdecl.mapaxes != NULL)
            mapaxes_kw = newEclKeyword_(MAPAXES_KW, ECL_FLOAT_TYPE);

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

ecl_kw_type *BlackoilEclipseOutputWriter::newEclKeyword_(const std::string& keyword, ecl_type_enum ecl_type) const
{
    ecl_kw_type* ecl_kw = NULL;

    if (ecl_type == ECL_INT_TYPE) {
        std::vector<int> data = eclipseParser_.getIntegerValue(keyword);
        ecl_kw = ecl_kw_alloc(keyword.c_str(), data.size(), ecl_type);
        ecl_kw_set_memcpy_data(ecl_kw, &data[0]);
    } else {
        std::vector<double> data = eclipseParser_.getFloatingPointValue(keyword);
        if (ecl_type == ECL_DOUBLE_TYPE) {
            ecl_kw = ecl_kw_alloc(keyword.c_str(), data.size(), ecl_type);
            ecl_kw_set_memcpy_data(ecl_kw, &data[0]);
        } else if (ecl_type == ECL_FLOAT_TYPE) {
            ecl_kw = ecl_kw_alloc(keyword.c_str(), data.size(), ecl_type);
            for (std::vector<double>::size_type i=0; i < data.size(); i++)
                ecl_kw_iset_float(ecl_kw, i, data[i]);
        }
    }
    return ecl_kw;
}

// TODO (?): unify with newEclKeyword_()
ecl_kw_type* BlackoilEclipseOutputWriter::eclKeywordWrapper_(const std::string& kw_name,
                                                             const std::vector<double> &data,
                                                             int offset,
                                                             int stride)
{
    if (stride <= 0)
      OPM_THROW(std::runtime_error, "Vector strides must be positive. Got stride = " << stride);
    if ((stride * std::vector<double>::size_type(grid_.number_of_cells)) != data.size())
      OPM_THROW(std::runtime_error, "Internal mismatch grid_.number_of_cells: " << grid_.number_of_cells << " data size: " << data.size() / stride);

    ecl_kw_type * ecl_kw = ecl_kw_alloc(kw_name.c_str(), grid_.number_of_cells, ECL_FLOAT_TYPE);
    for (int i=0; i < grid_.number_of_cells; i++)
        ecl_kw_iset_float(ecl_kw, i, data[i*stride + offset]);
    return ecl_kw;
}

void BlackoilEclipseOutputWriter::saveEclKeyword_(fortio_type* fortio, const std::string& kw, ecl_type_enum ecl_type)
{
    ecl_kw_type * ecl_kw = newEclKeyword_(kw, ecl_type);
    if (ecl_kw != NULL) {
        ecl_kw_fwrite(ecl_kw, fortio);
        ecl_kw_free(ecl_kw);
    }
}

#endif // HAVE_ERT

} // namespace Opm
