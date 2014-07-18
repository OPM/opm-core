/*
  Copyright (c) 2013-2014 Andreas Lauser
  Copyright (c) 2013 SINTEF ICT, Applied Mathematics.
  Copyright (c) 2013 Uni Research AS

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

#include "EclipseWriter.hpp"

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/cpgpreprocess/preprocess.h>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/parameters/Parameter.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/wells.h> // WellType

#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Utility/SpecgridWrapper.hpp>
#include <opm/parser/eclipse/Utility/WelspecsWrapper.hpp>

#include <boost/algorithm/string/case_conv.hpp> // to_upper_copy
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp> // path

#include <ctime>      // mktime
#include <forward_list>
#include <memory>     // unique_ptr
#include <utility>    // move

#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_endian_flip.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_rst_file.h>

// namespace start here since we don't want the ERT headers in it
namespace Opm {
namespace EclipseWriterDetails {

/// Helper function when we don't really want any transformation
/// (The C++ committee removed std::identity because it was "troublesome" (!?!)
static double noConversion(const double& u)
{ return u; }

/// Helper method that can be used in keyword transformation (must carry
/// the barsa argument)
static double toBar(const double& pressure)
{ return Opm::unit::convert::to(pressure, Opm::unit::barsa); }

/// Helper method that can be used in keyword transformation (must carry
/// the milliDarcy argument)
static double toMilliDarcy(const double& permeability)
{ return Opm::unit::convert::to(permeability, Opm::prefix::milli * Opm::unit::darcy); }

/// Names of the saturation property for each phase. The order of these
/// names are critical; they must be the same as the BlackoilPhases enum
static const char* saturationKeywordNames[] = { "SWAT", "SOIL", "SGAS" };

// throw away the data for all non-active cells in an array
void restrictToActiveCells(std::vector<double> &data, const std::vector<int> &actnumData)
{
    assert(actnumData.size() == data.size());

    size_t curActiveIdx = 0;
    for (size_t curIdx = 0; curIdx < data.size(); ++curIdx) {
        if (!actnumData[curIdx])
            continue; // ignore non-active cells

        assert(curActiveIdx <= curIdx);
        data[curActiveIdx] = data[curIdx];
        ++ curActiveIdx;
    }

    data.resize(curActiveIdx);
}

// throw away the data for all non-active cells in an array. (this is
// the variant of the function which takes an UnstructuredGrid object.)
void restrictToActiveCells(std::vector<double> &data,
                           int numCells,
                           const int* compressedToCartesianCellIdx)
{
    if (!compressedToCartesianCellIdx)
        // if there is no active -> global mapping, all cells
        // are considered active
        return;

    // activate those cells that are actually there
    for (int i = 0; i < numCells; ++i) {
        // make sure that global cell indices are always at least as
        // large as the active one and that the global cell indices
        // are in increasing order. the latter might become
        // problematic if cells are extensively re-ordered, but that
        // does not seem to be the case so far
        assert(compressedToCartesianCellIdx[i] >= i);
        assert(i == 0 || compressedToCartesianCellIdx[i - 1] < compressedToCartesianCellIdx[i]);

        data[i] = data[compressedToCartesianCellIdx[i]];
    }
    data.resize(numCells);
}

// convert the units of an array
template <class TransferFunction>
void convertUnit(std::vector<double> &data, TransferFunction &transferFn)
{
    for (size_t curIdx = 0; curIdx < data.size(); ++curIdx) {
        data[curIdx] = transferFn(data[curIdx]);
    }
}

// extract a sub-array of a larger one which represents multiple
// striped ones
void extractFromStripedData(std::vector<double> &data,
                            int offset,
                            int stride)
{
    size_t tmpIdx = 0;
    for (size_t curIdx = offset; curIdx < data.size(); curIdx += stride) {
        assert(tmpIdx <= curIdx);
        data[tmpIdx] = data[curIdx];
        ++tmpIdx;
    }
    // shirk the result
    data.resize(tmpIdx);
}

/// Convert OPM phase usage to ERT bitmask
int ertPhaseMask(const PhaseUsage uses)
{
    return (uses.phase_used[BlackoilPhases::Liquid] ? ECL_OIL_PHASE : 0)
        | (uses.phase_used[BlackoilPhases::Aqua] ? ECL_WATER_PHASE : 0)
        | (uses.phase_used[BlackoilPhases::Vapour] ? ECL_GAS_PHASE : 0);
}

/**
 * Eclipse "keyword" (i.e. named data) for a vector.
 */
template <typename T>
class Keyword : private boost::noncopyable
{
public:
    // Default constructor
    Keyword()
        : ertHandle_(0)
    {}

    /// Initialization from double-precision array.
    Keyword(const std::string& name,
            const std::vector<double>& data)
        : ertHandle_(0)
    { set(name, data); }

    /// Initialization from double-precision array.
    Keyword(const std::string& name,
            const std::vector<int>& data)
        : ertHandle_(0)
    { set(name, data); }

    ~Keyword()
    {
        if (ertHandle_)
            ecl_kw_free(ertHandle_);
    }

    template <class DataElementType>
    void set(const std::string name, const std::vector<DataElementType>& data)
    {
        if(ertHandle_) {
            ecl_kw_free(ertHandle_);
        }

        ertHandle_ = ecl_kw_alloc(name.c_str(),
                                  data.size(),
                                  ertType_());

        // number of elements to take
        const int numEntries = data.size();

        // fill it with values
        T* target = static_cast<T*>(ecl_kw_get_ptr(ertHandle()));
        for (int i = 0; i < numEntries; ++i) {
            target[i] = static_cast<T>(data[i]);
        }
    }

    ecl_kw_type *ertHandle() const
    { return ertHandle_; }

private:
    static ecl_type_enum ertType_()
    {
        if (std::is_same<T, float>::value)
        { return ECL_FLOAT_TYPE; }
        if (std::is_same<T, double>::value)
        { return ECL_DOUBLE_TYPE; }
        if (std::is_same<T, int>::value)
        { return ECL_INT_TYPE; }

        OPM_THROW(std::logic_error,
                  "Unhandled type for data elements in EclipseWriterDetails::Keyword");
    }

    ecl_kw_type *ertHandle_;
};

/**
 * Pointer to memory that holds the name to an Eclipse output file.
 */
class FileName : private boost::noncopyable
{
public:
    FileName(const std::string& outputDir,
             const std::string& baseName,
             ecl_file_enum type,
             int reportStepIdx)
    {
        ertHandle_ = ecl_util_alloc_filename(outputDir.c_str(),
                                             baseName.c_str(),
                                             type,
                                             false, // formatted?
                                             reportStepIdx);
    }

    ~FileName()
    { std::free(ertHandle_); }

    const char *ertHandle() const
    { return ertHandle_; }

private:
    char *ertHandle_;
};

class Restart : private boost::noncopyable
{
public:
    Restart(const std::string& outputDir,
            const std::string& baseName,
            int reportStepIdx)
    {
        restartFileName_ = ecl_util_alloc_filename(outputDir.c_str(),
                                                   baseName.c_str(),
                                                   /*type=*/ECL_UNIFIED_RESTART_FILE,
                                                   false, // use formatted instead of binary output?
                                                   reportStepIdx);

        if (reportStepIdx == 0) {
            restartFileHandle_ = ecl_rst_file_open_write(restartFileName_);
        }
        else {
            restartFileHandle_ = ecl_rst_file_open_append(restartFileName_);
        }
    }

    ~Restart()
    {
        free(restartFileName_);
        ecl_rst_file_close(restartFileHandle_);
    }

    void writeHeader(const SimulatorTimer& timer,
                     int reportStepIdx,
                     int numCells,
                     int nx,
                     int ny,
                     int nz,
                     const int *compressedToCartesianCellIdx,
                     const PhaseUsage uses)
    {
        ecl_rst_file_fwrite_header(restartFileHandle_,
                                   reportStepIdx,
                                   timer.currentPosixTime(),
                                   Opm::unit::convert::to(timer.simulationTimeElapsed(),
                                                          Opm::unit::day),
                                   nx, ny, nz,
                                   numCells,
                                   ertPhaseMask(uses));
    }

    ecl_rst_file_type *ertHandle() const
    { return restartFileHandle_; }

private:
    char *restartFileName_;
    ecl_rst_file_type *restartFileHandle_;
};

/**
 * The Solution class wraps the actions that must be done to the restart file while
 * writing solution variables; it is not a handle on its own.
 */
class Solution : private boost::noncopyable
{
public:
    Solution(Restart& restartHandle)
        : restartHandle_(&restartHandle)
    {  ecl_rst_file_start_solution(restartHandle_->ertHandle()); }

    ~Solution()
    { ecl_rst_file_end_solution(restartHandle_->ertHandle()); }

    template <typename T>
    void add(const Keyword<T>& kw)
    { ecl_rst_file_add_kw(restartHandle_->ertHandle(), kw.ertHandle()); }

    ecl_rst_file_type *ertHandle() const
    { return restartHandle_->ertHandle(); }

private:
    Restart* restartHandle_;
};

/// Supported well types. Enumeration doesn't let us get all the members,
/// so we must have an explicit array.
static WellType WELL_TYPES[] = { INJECTOR, PRODUCER };

class WellReport;

class Summary : private boost::noncopyable
{
public:
    Summary(const std::string& outputDir,
            const std::string& baseName,
            const SimulatorTimer& timer,
            int nx,
            int ny,
            int nz)
    {
        boost::filesystem::path casePath(outputDir);
        casePath /= boost::to_upper_copy(baseName);

        ertHandle_ = ecl_sum_alloc_writer(casePath.string().c_str(),
                                          false, /* formatted   */
                                          true,  /* unified     */
                                          ":",    /* join string */
                                          timer.simulationTimeElapsed(),
                                          nx,
                                          ny,
                                          nz);
    }

    ~Summary()
    { ecl_sum_free(ertHandle_); }

    typedef std::unique_ptr <WellReport> var_t;
    typedef std::vector <var_t> vars_t;

    Summary& addWell(var_t var)
    {
        vars_.push_back(std::move(var));
        return *this;
    }

    // no inline implementation of these two methods since they depend
    // on the classes defined in the following.

    // add rate variables for each of the well in the input file
    void addAllWells(Opm::EclipseStateConstPtr eclipseState,
                     const PhaseUsage& uses);
    void writeTimeStep(int reportStepIdx,
                       const SimulatorTimer& timer,
                       const WellState& wellState);

    ecl_sum_type *ertHandle() const
    { return ertHandle_; }

private:
    ecl_sum_type *ertHandle_;

    vars_t vars_;
};

class SummaryTimeStep : private boost::noncopyable
{
public:
    SummaryTimeStep(Summary& summaryHandle,
                    int reportStepIdx,
                    const SimulatorTimer &timer)
    {
        ertHandle_ = ecl_sum_add_tstep(summaryHandle.ertHandle(),
                                       reportStepIdx,
                                       Opm::unit::convert::to(timer.simulationTimeElapsed(),
                                                              Opm::unit::day));
    }

    // no destructor in this class as ERT takes care of freeing the
    // handle as part of freeing the solution handle!

    ecl_sum_tstep_type *ertHandle() const
    { return ertHandle_; };

private:
    ecl_sum_tstep_type *ertHandle_;
};


/**
 * Initialization file which contains static properties (such as
 * porosity and permeability) for the simulation field.
 */
class Init : private boost::noncopyable
{
public:
    Init(const std::string& outputDir,
         const std::string& baseName,
         int reportStepIdx)
        : egridFileName_(outputDir,
                         baseName,
                         ECL_EGRID_FILE,
                         reportStepIdx)
    {
        FileName initFileName(outputDir,
                              baseName,
                              ECL_INIT_FILE,
                              reportStepIdx);

        bool isFormatted;
        if (!ecl_util_fmt_file(initFileName.ertHandle(), &isFormatted)) {
            OPM_THROW(std::runtime_error,
                      "Could not determine formatted/unformatted status of file:" << initFileName.ertHandle() << " non-standard name?" << std::endl);
        }

        ertHandle_ = fortio_open_writer(initFileName.ertHandle(),
                                        isFormatted,
                                        ECL_ENDIAN_FLIP);
    }

    ~Init()
    { fortio_fclose(ertHandle_); }

    void writeHeader(int numCells,
                     const int* compressedToCartesianCellIdx,
                     const SimulatorTimer& timer,
                     Opm::EclipseStateConstPtr eclipseState,
                     const PhaseUsage uses)
    {
        auto dataField = eclipseState->getDoubleGridProperty("PORO")->getData();
        restrictToActiveCells(dataField, numCells, compressedToCartesianCellIdx);

        auto eclGrid = eclipseState->getEclipseGridCopy();

        // update the ACTNUM array using the processed cornerpoint grid
        std::vector<int> actnumData(eclGrid->getCartesianSize(), 1);
        if (compressedToCartesianCellIdx) {
            std::fill(actnumData.begin(), actnumData.end(), 0);
            for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
                int cartesianCellIdx = compressedToCartesianCellIdx[cellIdx];
                actnumData[cartesianCellIdx] = 1;
            }
        }
        eclGrid->resetACTNUM(&actnumData[0]);

        // finally, write the grid to disk
        eclGrid->fwriteEGRID(egridFileName_.ertHandle());

        Keyword<float> poro_kw("PORO", dataField);
        ecl_init_file_fwrite_header(ertHandle(),
                                    eclGrid->c_ptr(),
                                    poro_kw.ertHandle(),
                                    ertPhaseMask(uses),
                                    timer.currentPosixTime());
    }

    void writeKeyword(const std::string& keywordName, const std::vector<double> &data)
    {
        Keyword <float> kw(keywordName, data);
        ecl_kw_fwrite(kw.ertHandle(), ertHandle());
    }

    fortio_type *ertHandle() const
    { return ertHandle_; }

private:
    fortio_type *ertHandle_;
    FileName egridFileName_;
};

/**
 * Summary variable that reports a characteristics of a well.
 */
class WellReport : private boost::noncopyable
{
protected:
    // this is only needed to allow derived classes to hide their copy
    // constructor
    WellReport()
        : index_(0)
        , sign_(0)
    {}

    WellReport(const Summary& summary,    /* section to add to  */
               Opm::EclipseStateConstPtr eclipseState,/* well names         */
               int whichWell,                    /* index of well line */
               PhaseUsage uses,                  /* phases present     */
               BlackoilPhases::PhaseIndex phase, /* oil, water or gas  */
               WellType type,                    /* prod. or inj.      */
               char aggregation,                 /* rate or total      */
               std::string unit)
        // save these for when we update the value in a timestep
        : index_(whichWell * uses.num_phases + uses.phase_pos[phase])

        // producers can be seen as negative injectors
        , sign_(type == INJECTOR ? +1. : -1.)
    {
        ertHandle_ = ecl_sum_add_var(summary.ertHandle(),
                                     varName_(phase,
                                              type,
                                              aggregation).c_str(),
                                     wellName_(eclipseState, whichWell).c_str(),
                                     /*num=*/ 0,
                                     unit.c_str(),
                                     /*defaultValue=*/ 0.);
    }

public:
    /// Allows us to pass this type to ecl_sum_tstep_iset
    operator int()
    { return smspec_node_get_params_index(ertHandle()); }

    /// Update the monitor according to the new state of the well, and
    /// get the reported value. Note: Only call this once for each timestep.
    virtual double update(const SimulatorTimer& timer,
                          const WellState& wellState) = 0;

    smspec_node_type *ertHandle() const
    { return ertHandle_; }

private:
    smspec_node_type *ertHandle_;

    /// index into a (flattened) wells*phases matrix
    const int index_;

    /// natural sign of the rate
    const double sign_;

    /// Get the name associated with this well
    std::string wellName_(Opm::EclipseStateConstPtr eclipseState,
                          int whichWell)
    {
        return eclipseState->getSchedule()->getWells()[whichWell]->name();
    }

    /// Compose the name of the summary variable, e.g. "WOPR" for
    /// well oil production rate.
    std::string varName_(BlackoilPhases::PhaseIndex phase,
                         WellType type,
                         char aggregation)
    {
        std::string name;
        name += 'W'; // well
        if (aggregation == 'B') {
            name += "BHP";
        } else {
            switch (phase) {
            case BlackoilPhases::Aqua:   name += 'W'; break; /* water */
            case BlackoilPhases::Vapour: name += 'G'; break; /* gas */
            case BlackoilPhases::Liquid: name += 'O'; break; /* oil */
            default:
                OPM_THROW(std::runtime_error,
                          "Unknown phase used in blackoil reporting");
            }
            switch (type) {
            case WellType::INJECTOR: name += 'I'; break;
            case WellType::PRODUCER: name += 'P'; break;
            default:
                OPM_THROW(std::runtime_error,
                          "Unknown well type used in blackoil reporting");
            }
            name += aggregation; /* rate ('R') or total ('T') */
        }
        return name;
    }

protected:
    double rate(const WellState& wellState)
    {
        // convert m^3/s of injected fluid to m^3/d of produced fluid
        const double convFactor = Opm::unit::convert::to(1., Opm::unit::day);
        double value = 0;
        if (wellState.wellRates().size() > 0) {
            assert(int(wellState.wellRates().size()) > index_);
            value = sign_ * wellState.wellRates()[index_] * convFactor;
        }
        return value;
    }

    double bhp(const WellState& wellstate)
    {
        if (wellstate.bhp().size() > 0) {
            // Note that 'index_' is used here even though it is meant
            // to give a (well,phase) pair.
            const int num_phases = wellstate.wellRates().size() / wellstate.bhp().size();

            return wellstate.bhp()[index_/num_phases];
        }
        return 0.0;
    }
};

/// Monitors the rate given by a well.
class WellRate : public WellReport
{
public:
    WellRate(const Summary& summary,
             Opm::EclipseStateConstPtr eclipseState,
             int whichWell,
             PhaseUsage uses,
             BlackoilPhases::PhaseIndex phase,
             WellType type)
        : WellReport(summary,
                     eclipseState,
                     whichWell,
                     uses,
                     phase,
                     type,
                     'R',
                     "SM3/DAY" /* surf. cub. m. per day */)
    { }

    virtual double update(const SimulatorTimer& /*timer*/,
                          const WellState& wellState)
    {
        // TODO: Why only positive rates?
        return std::max (0., rate (wellState));
    }
};

/// Monitors the total production in a well.
class WellTotal : public WellReport
{
public:
    WellTotal(const Summary& summary,
              Opm::EclipseStateConstPtr eclipseState,
              int whichWell,
              PhaseUsage uses,
              BlackoilPhases::PhaseIndex phase,
              WellType type)
        : WellReport(summary,
                     eclipseState,
                     whichWell,
                     uses,
                     phase,
                     type,
                     'T',
                     "SM3" /* surface cubic meter */ )
          // nothing produced when the reporting starts
        , total_(0.)
    { }

    virtual double update(const SimulatorTimer& timer,
                          const WellState& wellState)
    {
        if (timer.currentStepNum() == 0) {
            // We are at the initial state.
            // No step has been taken yet.
            return 0.0;
        }
        // TODO: Is the rate average for the timestep, or is in
        // instantaneous (in which case trapezoidal or Simpson integration
        // would probably be better)
        const double intg = timer.stepLengthTaken() * rate(wellState);
        // add this timesteps production to the total
        total_ += intg;
        // report the new production total
        return total_;
    }

private:
    /// Aggregated value of the course of the simulation
    double total_;
};

/// Monitors the bottom hole pressure in a well.
class WellBhp : public WellReport
{
public:
    WellBhp(const Summary& summary,
            Opm::EclipseStateConstPtr eclipseState,
            int whichWell,
            PhaseUsage uses,
            BlackoilPhases::PhaseIndex phase,
            WellType type)
        : WellReport(summary,
                     eclipseState,
                     whichWell,
                     uses,
                     phase,
                     type,
                     'B',
                     "Pascal")
    { }

    virtual double update(const SimulatorTimer& /*timer*/,
                          const WellState& wellState)
    {
        return bhp(wellState);
    }
};

// no inline implementation of this since it depends on the
// WellReport type being completed first
void Summary::writeTimeStep(int reportStepIdx,
                            const SimulatorTimer& timer,
                            const WellState& wellState)
{
    // internal view; do not move this code out of Summary!
    SummaryTimeStep tstep(*this, reportStepIdx, timer);
    // write all the variables
    for (vars_t::iterator v = vars_.begin(); v != vars_.end(); ++v) {
        const double value = (*v)->update(timer, wellState);
        ecl_sum_tstep_iset(tstep.ertHandle(), *(*v).get(), value);
    }

    // write the summary file to disk
    ecl_sum_fwrite(ertHandle());
}

void Summary::addAllWells(Opm::EclipseStateConstPtr eclipseState,
                          const PhaseUsage& uses)
{
    // TODO: Only create report variables that are requested with keywords
    // (e.g. "WOPR") in the input files, and only for those wells that are
    // mentioned in those keywords
    const int numWells = eclipseState->getSchedule()->numWells();
    for (int phaseIdx = 0; phaseIdx != BlackoilPhases::MaxNumPhases; ++phaseIdx) {
        const BlackoilPhases::PhaseIndex ertPhaseIdx =
            static_cast <BlackoilPhases::PhaseIndex>(phaseIdx);
        // don't bother with reporting for phases that aren't there
        if (!uses.phase_used[phaseIdx]) {
            continue;
        }
        for (size_t typeIndex = 0;
             typeIndex < sizeof(WELL_TYPES) / sizeof(WELL_TYPES[0]);
             ++typeIndex) {
            const WellType type = WELL_TYPES[typeIndex];
            for (int whichWell = 0; whichWell != numWells; ++whichWell) {
                // W{O,G,W}{I,P}R
                addWell(std::unique_ptr <WellReport>(
                            new WellRate(*this,
                                         eclipseState,
                                         whichWell,
                                         uses,
                                         ertPhaseIdx,
                                         type)));
                // W{O,G,W}{I,P}T
                addWell(std::unique_ptr <WellReport>(
                            new WellTotal(*this,
                                          eclipseState,
                                          whichWell,
                                          uses,
                                          ertPhaseIdx,
                                          type)));
            }
        }
    }

    // Add BHP monitors
    for (int whichWell = 0; whichWell != numWells; ++whichWell) {
        // In the call below: uses, phase and the well type arguments
        // are not used, except to set up an index that stores the
        // well indirectly. For details see the implementation of the
        // WellReport constructor, and the method
        // WellReport::bhp().
        BlackoilPhases::PhaseIndex ertPhaseIdx = BlackoilPhases::Liquid;
        if (!uses.phase_used[BlackoilPhases::Liquid]) {
            ertPhaseIdx = BlackoilPhases::Vapour;
        }
        addWell(std::unique_ptr <WellReport>(
                    new WellBhp(*this,
                                eclipseState,
                                whichWell,
                                uses,
                                ertPhaseIdx,
                                WELL_TYPES[0])));
    }
}
} // end namespace EclipseWriterDetails

void EclipseWriter::writeInit(const SimulatorTimer &timer)
{
    // if we don't want to write anything, this method becomes a
    // no-op...
    if (!enableOutput_) {
        return;
    }

    reportStepIdx_ = 0;

    EclipseWriterDetails::Init fortio(outputDir_, baseName_, /*stepIdx=*/0);
    fortio.writeHeader(numCells_,
                       compressedToCartesianCellIdx_,
                       timer,
                       eclipseState_,
                       phaseUsage_);

    if (eclipseState_->hasDoubleGridProperty("PERMX")) {
        auto data = eclipseState_->getDoubleGridProperty("PERMX")->getData();
        EclipseWriterDetails::convertUnit(data, EclipseWriterDetails::toMilliDarcy);
        fortio.writeKeyword("PERMX", data);
    }
    if (eclipseState_->hasDoubleGridProperty("PERMY")) {
        auto data = eclipseState_->getDoubleGridProperty("PERMY")->getData();
        EclipseWriterDetails::convertUnit(data, EclipseWriterDetails::toMilliDarcy);
        fortio.writeKeyword("PERMY", data);
    }
    if (eclipseState_->hasDoubleGridProperty("PERMZ")) {
        auto data = eclipseState_->getDoubleGridProperty("PERMZ")->getData();
        EclipseWriterDetails::convertUnit(data, EclipseWriterDetails::toMilliDarcy);
        fortio.writeKeyword("PERMZ", data);
    }

    /* Create summary object (could not do it at construction time,
       since it requires knowledge of the start time). */
    auto eclGrid = eclipseState_->getEclipseGrid();
    summary_.reset(new EclipseWriterDetails::Summary(outputDir_,
                                                     baseName_,
                                                     timer,
                                                     eclGrid->getNX(),
                                                     eclGrid->getNY(),
                                                     eclGrid->getNZ()));
    summary_->addAllWells(eclipseState_, phaseUsage_);
}

void EclipseWriter::writeTimeStep(const SimulatorTimer& timer,
                                  const SimulatorState& reservoirState,
                                  const WellState& wellState)
{
    // if we don't want to write anything, this method becomes a
    // no-op...
    if (!enableOutput_) {
        return;
    }

    // respected the output_interval parameter
    if (reportStepIdx_ % outputInterval_ != 0) {
        return;
    }

    // start writing to files
    EclipseWriterDetails::Restart restartHandle(outputDir_, baseName_, reportStepIdx_);
    restartHandle.writeHeader(timer,
                              reportStepIdx_,
                              numCells_,
                              cartesianSize_[0],
                              cartesianSize_[1],
                              cartesianSize_[2],
                              compressedToCartesianCellIdx_,
                              phaseUsage_);
    EclipseWriterDetails::Solution sol(restartHandle);

    // write out the pressure of the reference phase (whatever
    // phase that is...). this is not the most performant solution
    // thinkable, but this is also not in the most performance
    // critical code path!
    std::vector<double> tmp = reservoirState.pressure();
    EclipseWriterDetails::convertUnit(tmp, EclipseWriterDetails::toBar);

    sol.add(EclipseWriterDetails::Keyword<float>("PRESSURE", tmp));

    for (int phase = 0; phase != BlackoilPhases::MaxNumPhases; ++phase) {
        // Eclipse never writes the oil saturation, so all post-processors
        // must calculate this from the other saturations anyway
        if (phase == BlackoilPhases::PhaseIndex::Liquid) {
            continue;
        }
        if (phaseUsage_.phase_used[phase]) {
            tmp = reservoirState.saturation();
            EclipseWriterDetails::extractFromStripedData(tmp,
                                                         /*offset=*/phaseUsage_.phase_pos[phase],
                                                         /*stride=*/phaseUsage_.num_phases);
            sol.add(EclipseWriterDetails::Keyword<float>(EclipseWriterDetails::saturationKeywordNames[phase], tmp));
        }
    }

    /* Summary variables (well reporting) */
    // TODO: instead of writing the header (smspec) every time, it should
    // only be written when there is a change in the well configuration
    // (first timestep, in practice), and reused later. but how to do this
    // without keeping the complete summary in memory (which will then
    // accumulate all the timesteps)?
    //
    // Note: The answer to the question above is still not settled, but now we do keep
    // the complete summary in memory, as a member variable in the EclipseWriter class,
    // instead of creating a temporary EclipseWriterDetails::Summary in this function
    // every time it is called.  This has been changed so that the final summary file
    // will contain data from the whole simulation, instead of just the last step.
    summary_->writeTimeStep(reportStepIdx_, timer, wellState);

    ++reportStepIdx_;
}


EclipseWriter::EclipseWriter(const parameter::ParameterGroup& params,
                             Opm::EclipseStateConstPtr eclipseState,
                             const Opm::PhaseUsage &phaseUsage,
                             int numCells,
                             const int* compressedToCartesianCellIdx)
    : eclipseState_(eclipseState)
    , numCells_(numCells)
    , compressedToCartesianCellIdx_(compressedToCartesianCellIdx)
    , phaseUsage_(phaseUsage)
{
    const auto eclGrid = eclipseState->getEclipseGrid();
    cartesianSize_[0] = eclGrid->getNX();
    cartesianSize_[1] = eclGrid->getNY();
    cartesianSize_[2] = eclGrid->getNZ();

    init(params);
}

void EclipseWriter::init(const parameter::ParameterGroup& params)
{
    // get the base name from the name of the deck
    using boost::filesystem::path;
    path deckPath(params.get <std::string>("deck_filename"));
    if (boost::to_upper_copy(path(deckPath.extension()).string()) == ".DATA") {
        baseName_ = path(deckPath.stem()).string();
    }
    else {
        baseName_ = path(deckPath.filename()).string();
    }

    // make uppercase of everything (or otherwise we'll get uppercase
    // of some of the files (.SMSPEC, .UNSMRY) and not others
    baseName_ = boost::to_upper_copy(baseName_);

    // retrieve the value of the "output" parameter
    enableOutput_ = params.getDefault<bool>("output", /*defaultValue=*/true);

    // retrieve the interval at which something should get written to
    // disk (once every N timesteps)
    outputInterval_ = params.getDefault<int>("output_interval", /*defaultValue=*/1);

    // store in current directory if not explicitly set
    outputDir_ = params.getDefault<std::string>("output_dir", ".");

    // set the index of the first time step written to 0...
    reportStepIdx_ = 0;

    if (enableOutput_) {
        // make sure that the output directory exists, if not try to create it
        if (!boost::filesystem::exists(outputDir_)) {
            std::cout << "Trying to create directory \"" << outputDir_ << "\" for the simulation output\n";
            boost::filesystem::create_directories(outputDir_);
        }

        if (!boost::filesystem::is_directory(outputDir_)) {
            OPM_THROW(std::runtime_error,
                      "The path specified as output directory '" << outputDir_
                      << "' is not a directory");
        }
    }
}

// default destructor is OK, just need to be defined
EclipseWriter::~EclipseWriter()
{ }

} // namespace Opm
