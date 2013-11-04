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

#include <memory>     // unique_ptr
#include <utility>    // move
using namespace Opm;

/// Smart pointer/handle class for ERT opaque types, such as ecl_kw_type*.
///
/// \tparam T Type of handle being wrapper
template <typename T>
struct EclipseHandle : private std::unique_ptr <T, void (*)(T*) throw()> {
    // to save ourselves from some typing we introduce this alias
    typedef std::unique_ptr <T, void (*)(T*) throw()> base;

    /// Construct a new smart handle based on the returned value of
    /// an allocation function and the corresponding destroyer function.
    EclipseHandle (T* t, void (*destroy)(T*))
        : base (t, destroy) { }

    /// Convenience operator that lets us use this type as if
    /// it was a handle directly.
    operator T* () const { return base::get (); }
};

/**
 * Eclipse "keyword" (i.e. named data) for a vector. (This class is
 * different from EclKW in the constructors it provide).
 */
template <typename T>
struct EclipseKeyword : public EclipseHandle <ecl_kw_type> {
    EclipseKeyword (const std::string& name,    /// identification
                    const std::vector<T>& data, /// array holding values
                    const int num,              /// actual number to take
                    const int offset,           /// distance to first
                    const int stride)           /// distance between each

        // allocate handle and put in smart pointer base class
        : EclipseHandle (ecl_kw_alloc (name.c_str(), num, type ()),
                         ecl_kw_free) {

        // range cannot start outside of data set
        assert(offset >= 0 && offset < data.size());

        // don't jump out of the set when trying to
        assert(stride > 0 && stride < data.size() - offset);

        // fill it with values
        for (int i = 0; i < num; ++i) {
            // access from data store
            const float value = data[i * stride + offset];

            // write into memory represented by handle
            ecl_kw_iset_float(*this, i, value);
        }
    }

    /// Convenience constructor that takes the entire array
    EclipseKeyword (const std::string& name,
                    const std::vector<T>& data)
        : EclipseKeyword (name, data, data.size(), 0, 1) { }

    /// Convenience constructor that gets the set of data
    /// from the samely named item in the parser
    EclipseKeyword (const std::string& name,
                    const EclipseGridParser& parser)
        : EclipseKeyword (name, parser.getValue<T> (name)) { }

    /// Constructor for optional fields
    EclipseKeyword (const std::string& name)
        : EclipseHandle (0, ecl_kw_free) {
        static_cast<void> (name);
    }

private:
    /// Map the C++ data type (given by T) to an Eclipse type enum
    static ecl_type_enum type ();
};

// specializations for known keyword types
template <> ecl_type_enum EclipseKeyword<int   >::type () { return ECL_INT_TYPE   ; }
template <> ecl_type_enum EclipseKeyword<double>::type () { return ECL_FLOAT_TYPE; }

/**
 * Extract the current time from a timer object into the C type used by ERT.
 */
time_t current (const SimulatorTimer& timer) {
    tm t = boost::posix_time::to_tm (timer.currentDateTime());
    return mktime(&t);
}

namespace Opm {
void BlackoilEclipseOutputWriter::writeInitFile(const SimulatorTimer &timer)
{
    startTime_ = current (timer);

#if HAVE_ERT
    writeGridInitFile_(timer);
    writeSummaryHeaderFile_(timer);
#else
    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
#endif // HAVE_ERT
}

/**
 * Pointer to memory that holds the name to an Eclipse output file.
 */
struct EclipseFileName : public EclipseHandle <const char> {
    EclipseFileName (const std::string& outputDir,
                     const std::string& baseName,
                     ecl_file_enum type,
                     const SimulatorTimer& timer)

        // filename formatting function returns a pointer to allocated
        // memory that must be released with the free() function
        : EclipseHandle (ecl_util_alloc_filename (outputDir.c_str(),
                                                  baseName.c_str(),
                                                  type,
                                                  false, // formatted?
                                                  timer.currentStepNum ()),
                         EclipseFileName::freestr) { }
private:
    /// Facade which allows us to free a const char*
    static void freestr (const char* ptr) {
        ::free (const_cast<char*>(ptr));
    }
};

void BlackoilEclipseOutputWriter::writeReservoirState(const BlackoilState& reservoirState, const SimulatorTimer& timer)
{
#if HAVE_ERT
    ecl_file_enum file_type = ECL_UNIFIED_RESTART_FILE;  // Alternatively ECL_RESTART_FILE for multiple restart files.

    EclipseFileName fileName (outputDir_,
                              baseName_,
                              ECL_UNIFIED_RESTART_FILE,
                              timer);
    int phases = ECL_OIL_PHASE + ECL_GAS_PHASE + ECL_WATER_PHASE;
    double days = Opm::unit::convert::to(timer.currentTime(), Opm::unit::day);
    int nx = grid_.cartdims[0];
    int ny = grid_.cartdims[1];
    int nz = grid_.cartdims[2];
    int nactive = grid_.number_of_cells;
    ecl_rst_file_type* rst_file;

    time_t curTime = current (timer);

    if (timer.currentStepNum() > 0 && file_type == ECL_UNIFIED_RESTART_FILE)
        rst_file = ecl_rst_file_open_append(fileName);
    else
        rst_file = ecl_rst_file_open_write(fileName);

    ecl_rst_file_fwrite_header(rst_file, timer.currentStepNum(), curTime, days, nx, ny, nz, nactive, phases);
    ecl_rst_file_start_solution(rst_file);

    {
        // convert the pressures from Pascals to bar because eclipse
        // seems to write bars
        std::vector<double> pressureBar(reservoirState.pressure());
        auto it = pressureBar.begin();
        const auto &endIt = pressureBar.end();
        for (; it != endIt; ++it)
            (*it) /= 1e5;

        EclipseKeyword<double> pressure_kw ("PRESSURE", pressureBar);
        ecl_rst_file_add_kw(rst_file, pressure_kw);
    }

    {
        EclipseKeyword<double> swat_kw ("SWAT",
                                         reservoirState.saturation(),
                                         grid_.number_of_cells,
                                         BlackoilPhases::Aqua,
                                         BlackoilPhases::MaxNumPhases);
        ecl_rst_file_add_kw(rst_file, swat_kw);
    }

    {
        EclipseKeyword<double> soil_kw ("SOIL",
                                         reservoirState.saturation(),
                                         grid_.number_of_cells,
                                         BlackoilPhases::Liquid,
                                         BlackoilPhases::MaxNumPhases);
        ecl_rst_file_add_kw(rst_file, soil_kw);
    }

    {
        EclipseKeyword<double> sgas_kw ("SGAS",
                                         reservoirState.saturation(),
                                         grid_.number_of_cells,
                                         BlackoilPhases::Vapour,
                                         BlackoilPhases::MaxNumPhases);
        ecl_rst_file_add_kw(rst_file, sgas_kw);
    }

    ecl_rst_file_end_solution(rst_file);
    ecl_rst_file_close(rst_file);
#else
    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
#endif // HAVE_ERT
}

void BlackoilEclipseOutputWriter::writeWellState(const WellState& wellState, const SimulatorTimer& timer)
{
#if HAVE_ERT
    time_t curTime = current (timer);

    // create a new timestep for the summary file (at least if the
    // timer was advanced since the last call to writeWellState())
    ecl_sum_tstep_type* tstep=
        ecl_sum_add_tstep(sumWriter_,
                          timer.currentStepNum() + 1,
                          (curTime - startTime_)/(24*60*60));

    int numWells = accumulatedProducedFluids_.size();
    for (int wellIdx = 0; wellIdx < numWells; ++wellIdx) {
        // set the value for the well oil production rate. For this,
        // be aware that the rates in the well state are _surface_
        // volume rates...
        double woprValue = wellState.wellRates()[wellIdx*3 + BlackoilPhases::Liquid];
        woprValue *= - 1 * (24 * 60 * 60); // convert m^3/s of injected fluid to m^3/d of produced fluid
        woprValue = std::max(0.0, woprValue);

        int woprIdx = smspec_node_get_params_index(woprSmspec_[wellIdx]);
        ecl_sum_tstep_iset(tstep, woprIdx, woprValue);

        // set the value for the well gas production rate
        double wgprValue = wellState.wellRates()[wellIdx*3 + BlackoilPhases::Vapour];
        wgprValue *= - 1 * (24 * 60 * 60); // convert m^3/s of injected fluid to m^3/d of produced fluid
        wgprValue = std::max(0.0, wgprValue);

        int wgprIdx = smspec_node_get_params_index(wgprSmspec_[wellIdx]);
        ecl_sum_tstep_iset(tstep, wgprIdx, wgprValue);

        // water injection rate
        double wwirValue = wellState.wellRates()[wellIdx*3 + BlackoilPhases::Aqua];
        wwirValue *= 1 * (24 * 60 * 60);  // convert m^3/s to m^3/d
        wwirValue = std::max(0.0, wwirValue);

        int wwirIdx = smspec_node_get_params_index(wwirSmspec_[wellIdx]);
        ecl_sum_tstep_iset(tstep, wwirIdx, wwirValue);

        // gas injection rate
        double wgirValue = wellState.wellRates()[wellIdx*3 + BlackoilPhases::Vapour];
        wgirValue *= - 1 * (24 * 60 * 60); // convert m^3/s of injected fluid to m^3/d of produced fluid
        wgirValue = std::max(0.0, wgirValue);

        int wgirIdx = smspec_node_get_params_index(wgirSmspec_[wellIdx]);
        ecl_sum_tstep_iset(tstep, wgirIdx, wgirValue);

        // accumulate injected produced fluids
        for (int phaseIdx = 0; phaseIdx < /*numPhases=*/3; ++phaseIdx) {
            // accumulate the produced/injected surface volumes
            double injectedVolume = wellState.wellRates()[wellIdx*3 + phaseIdx];
            injectedVolume *= timer.currentStepLength();

            if (injectedVolume < 0)
                accumulatedProducedFluids_[wellIdx][phaseIdx] += -injectedVolume;
            else
                accumulatedInjectedFluids_[wellIdx][phaseIdx] += injectedVolume;

            int woptIdx = smspec_node_get_params_index(woptSmspec_[wellIdx]);
            ecl_sum_tstep_iset(tstep, woptIdx, accumulatedProducedFluids_[wellIdx][BlackoilPhases::Liquid]);

            int wgptIdx = smspec_node_get_params_index(wgptSmspec_[wellIdx]);
            ecl_sum_tstep_iset(tstep, wgptIdx, accumulatedProducedFluids_[wellIdx][BlackoilPhases::Vapour]);

            int wwitIdx = smspec_node_get_params_index(wwitSmspec_[wellIdx]);
            ecl_sum_tstep_iset(tstep, wwitIdx, accumulatedInjectedFluids_[wellIdx][BlackoilPhases::Aqua]);

            int wgitIdx = smspec_node_get_params_index(wgitSmspec_[wellIdx]);
            ecl_sum_tstep_iset(tstep, wgitIdx, accumulatedProducedFluids_[wellIdx][BlackoilPhases::Vapour]);
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

    ecl_grid_type* ecl_grid = newEclGrid_();
    EclipseFileName gridFileName (outputDir_,
                                  baseName_,
                                  ECL_EGRID_FILE,
                                  timer);
    fortio_type* fortio;
    ecl_grid_fwrite_EGRID(ecl_grid, gridFileName);

    EclipseFileName initFileName (outputDir_,
                                  baseName_,
                                  ECL_INIT_FILE,
                                  timer);
    if (!ecl_util_fmt_file(initFileName, &fmt_file)) {
        OPM_THROW(std::runtime_error,
                  "Could not determine formatted/unformatted status of file:" << initFileName << " non-standard name?" << std::endl);
    }
    fortio = fortio_open_writer(initFileName, fmt_file, endian_flip);
    {
        time_t start_date = current (timer);
        EclipseKeyword<double> poro_kw (PORO_KW, eclipseParser_);
        ecl_init_file_fwrite_header(fortio, ecl_grid, poro_kw, phases, start_date);
    }

    /* This collection of keywords is somewhat arbitrary and random. */
    EclipseKeyword<double> permx_kw ("PERMX", eclipseParser_);
    EclipseKeyword<double> permy_kw ("PERMY", eclipseParser_);
    EclipseKeyword<double> permz_kw ("PERMZ", eclipseParser_);
    ecl_kw_fwrite(permx_kw, fortio);
    ecl_kw_fwrite(permy_kw, fortio);
    ecl_kw_fwrite(permz_kw, fortio);

    fortio_fclose(fortio);
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
                             /*formattedOutput=*/false,
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
    wgprSmspec_.resize(numWells);
    wgptSmspec_.resize(numWells);
    wwirSmspec_.resize(numWells);
    wwitSmspec_.resize(numWells);
    wgirSmspec_.resize(numWells);
    wgitSmspec_.resize(numWells);
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

        wgprSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                               /*varName=*/"WGPR",
                                               /*wellGroupName=*/wellIt->name_.c_str(),
                                               /*num=*/0,
                                               /*unit=*/"SM3/DAY",
                                               /*defaultValue=*/0.0);

        wgptSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                               /*varName=*/"WGPT",
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

        wgirSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                               /*varName=*/"WGIR",
                                               /*wellGroupName=*/wellIt->name_.c_str(),
                                               /*num=*/0,
                                               /*unit=*/"SM3/DAY",
                                               /*defaultValue=*/0.0);

        wgitSmspec_[wellIdx] = ecl_sum_add_var(sumWriter_,
                                               /*varName=*/"WGIT",
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

        EclipseKeyword<double> coord_kw   (COORD_KW,  eclipseParser_);
        EclipseKeyword<double> zcorn_kw   (ZCORN_KW,  eclipseParser_);
        EclipseKeyword<double> actnum_kw  (ACTNUM_KW, eclipseParser_);
        EclipseKeyword<double> mapaxes_kw (MAPAXES_KW);

        ecl_grid_type * grid ;
        if (grdecl.mapaxes != NULL)
            mapaxes_kw = std::move (EclipseKeyword<double> (MAPAXES_KW, eclipseParser_));

        grid = ecl_grid_alloc_GRDECL_kw(grdecl.dims[0], grdecl.dims[1], grdecl.dims[2], zcorn_kw, coord_kw, actnum_kw, mapaxes_kw);

        return grid;
    }
    OPM_THROW(std::runtime_error,
              "Can't create an ERT grid (no supported keywords found in deck)");
}

#endif // HAVE_ERT

} // namespace Opm
