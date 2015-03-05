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
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/cpgpreprocess/preprocess.h>
#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/simulator/SimulatorTimerInterface.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/io/eclipse/EclipseWriteRFTHandler.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/parameters/Parameter.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/wells.h> // WellType

#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Utility/SpecgridWrapper.hpp>
#include <opm/parser/eclipse/Utility/WelspecsWrapper.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

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
#include <ert/ecl_well/well_const.h>
#include <ert/ecl/ecl_rsthead.h>

// namespace start here since we don't want the ERT headers in it
namespace Opm {
namespace EclipseWriterDetails {
/// Names of the saturation property for each phase. The order of these
/// names are critical; they must be the same as the BlackoilPhases enum
static const char* saturationKeywordNames[] = { "SWAT", "SOIL", "SGAS" };

// throw away the data for all non-active cells and reorder to the Cartesian logic of
// eclipse
void restrictAndReorderToActiveCells(std::vector<double> &data,
                                     int numCells,
                                     const int* compressedToCartesianCellIdx)
{
    if (!compressedToCartesianCellIdx)
        // if there is no active -> global mapping, all cells
        // are considered active
        return;

    std::vector<double> eclData;
    eclData.reserve( numCells );

    // activate those cells that are actually there
    for (int i = 0; i < numCells; ++i) {
        eclData.push_back( data[ compressedToCartesianCellIdx[i] ] );
    }
    data.swap( eclData );
}

// convert the units of an array
void convertFromSiTo(std::vector<double> &siValues, double toSiConversionFactor)
{
    for (size_t curIdx = 0; curIdx < siValues.size(); ++curIdx) {
        siValues[curIdx] = unit::convert::to(siValues[curIdx], toSiConversionFactor);
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

    Keyword(const std::string& name,
            const std::vector<const char*>& data)
        : ertHandle_(0)
    {set(name, data); }

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

    void set(const std::string name, const std::vector<const char *>& data)
    {
      if(ertHandle_) {
          ecl_kw_free(ertHandle_);
      }


      ertHandle_ = ecl_kw_alloc(name.c_str(),
                                data.size(),
                                ertType_());

      // number of elements to take
      const int numEntries = data.size();
      for (int i = 0; i < numEntries; ++i) {
          ecl_kw_iset_char_ptr( ertHandle_, i, data[i]);
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
        if (std::is_same<T, const char *>::value)
        { return ECL_CHAR_TYPE; }


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
             int writeStepIdx)
    {
        ertHandle_ = ecl_util_alloc_filename(outputDir.c_str(),
                                             baseName.c_str(),
                                             type,
                                             false, // formatted?
                                             writeStepIdx);
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
    static const int NIWELZ = 11; //Number of data elements per well in IWEL array in restart file
    static const int NZWELZ = 3;  //Number of 8-character words per well in ZWEL array restart file
    static const int NICONZ = 15; //Number of data elements per completion in ICON array restart file

    /**
     * The constants NIWELZ and NZWELZ referes to the number of elements per
     * well that we write to the IWEL and ZWEL eclipse restart file data
     * arrays. The constant NICONZ refers to the number of elements per
     * completion in the eclipse restart file ICON data array.These numbers are
     * written to the INTEHEAD header.

     * The elements are added in the method addRestartFileIwelData(...) and and
     * addRestartFileIconData(...), respectively.  We write as many elements
     * that we need to be able to view the restart file in Resinsight.  The
     * restart file will not be possible to restart from with Eclipse, we write
     * to little information to be able to do this.
     *
     * Observe that all of these values are our "current-best-guess" for how
     * many numbers are needed; there might very well be third party
     * applications out there which have a hard expectation for these values.
    */


    Restart(const std::string& outputDir,
            const std::string& baseName,
            int writeStepIdx)
    {
        restartFileName_ = ecl_util_alloc_filename(outputDir.c_str(),
                                                   baseName.c_str(),
                                                   /*type=*/ECL_UNIFIED_RESTART_FILE,
                                                   false, // use formatted instead of binary output?
                                                   writeStepIdx);

        if (writeStepIdx == 0) {
            restartFileHandle_ = ecl_rst_file_open_write(restartFileName_);
        }
        else {
            restartFileHandle_ = ecl_rst_file_open_append(restartFileName_);
        }
    }

    template <typename T>
    void add_kw(const Keyword<T>& kw)
    { ecl_rst_file_add_kw(restartFileHandle_, kw.ertHandle()); }


    void addRestartFileIwelData(std::vector<int>& iwel_data, size_t currentStep , WellConstPtr well, size_t offset) const {
        CompletionSetConstPtr completions = well->getCompletions( currentStep );

        iwel_data[ offset + IWEL_HEADI_ITEM ] = well->getHeadI() + 1;
        iwel_data[ offset + IWEL_HEADJ_ITEM ] = well->getHeadJ() + 1;
        iwel_data[ offset + IWEL_CONNECTIONS_ITEM ] = completions->size();
        iwel_data[ offset + IWEL_GROUP_ITEM ] = 1;

        {
            WellType welltype = well->isProducer(currentStep) ? PRODUCER : INJECTOR;
            int ert_welltype = EclipseWriter::eclipseWellTypeMask(welltype, well->getInjectionProperties(currentStep).injectorType);
            iwel_data[ offset + IWEL_TYPE_ITEM ] = ert_welltype;
        }

        iwel_data[ offset + IWEL_STATUS_ITEM ] = EclipseWriter::eclipseWellStatusMask(well->getStatus(currentStep));
    }




    void addRestartFileIconData(std::vector<int>& icon_data, CompletionSetConstPtr completions , size_t wellICONOffset) const {
        for (size_t i = 0; i < completions->size(); ++i) {
            CompletionConstPtr completion = completions->get(i);
            size_t iconOffset = wellICONOffset + i * Opm::EclipseWriterDetails::Restart::NICONZ;
            icon_data[ iconOffset + ICON_IC_ITEM] = 1;

            icon_data[ iconOffset + ICON_I_ITEM] = completion->getI() + 1;
            icon_data[ iconOffset + ICON_J_ITEM] = completion->getJ() + 1;
            icon_data[ iconOffset + ICON_K_ITEM] = completion->getK() + 1;

            {
                WellCompletion::StateEnum completion_state = completion->getState();
                if (completion_state == WellCompletion::StateEnum::OPEN) {
                    icon_data[ iconOffset + ICON_STATUS_ITEM ] = 1;
                } else {
                    icon_data[ iconOffset + ICON_STATUS_ITEM ] = 0;
                }
            }
            icon_data[ iconOffset + ICON_DIRECTION_ITEM] = (int)completion->getDirection();
        }
    }


    ~Restart()
    {
        free(restartFileName_);
        ecl_rst_file_close(restartFileHandle_);
    }

    void writeHeader(const SimulatorTimerInterface& /*timer*/,
                     int writeStepIdx,
                     ecl_rsthead_type * rsthead_data)
    {

      ecl_rst_file_fwrite_header(restartFileHandle_,
                                 writeStepIdx,
                                 rsthead_data);

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
            const SimulatorTimerInterface& timer,
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

    typedef std::unique_ptr<WellReport> SummaryReportVar;
    typedef std::vector<SummaryReportVar> SummaryReportVarCollection;

    Summary& addWell(SummaryReportVar var)
    {
        summaryReportVars_.push_back(std::move(var));
        return *this;
    }

    // no inline implementation of these two methods since they depend
    // on the classes defined in the following.

    // add rate variables for each of the well in the input file
    void addAllWells(Opm::EclipseStateConstPtr eclipseState,
                     const PhaseUsage& uses);
    void writeTimeStep(int writeStepIdx,
                       const SimulatorTimerInterface& timer,
                       const WellState& wellState);

    ecl_sum_type *ertHandle() const
    { return ertHandle_; }

private:
    ecl_sum_type *ertHandle_;

    Opm::EclipseStateConstPtr eclipseState_;
    SummaryReportVarCollection summaryReportVars_;
};

class SummaryTimeStep : private boost::noncopyable
{
public:
    SummaryTimeStep(Summary& summaryHandle,
                    int writeStepIdx,
                    const SimulatorTimerInterface &timer)
    {
        ertHandle_ = ecl_sum_add_tstep(summaryHandle.ertHandle(),
                                       writeStepIdx,
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
         int writeStepIdx)
        : egridFileName_(outputDir,
                         baseName,
                         ECL_EGRID_FILE,
                         writeStepIdx)
    {
        FileName initFileName(outputDir,
                              baseName,
                              ECL_INIT_FILE,
                              writeStepIdx);

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
                     const SimulatorTimerInterface& timer,
                     Opm::EclipseStateConstPtr eclipseState,
                     const PhaseUsage uses)
    {
        auto dataField = eclipseState->getDoubleGridProperty("PORO")->getData();
        restrictAndReorderToActiveCells(dataField, numCells, compressedToCartesianCellIdx);

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
    WellReport(const Summary& summary,    /* section to add to  */
               Opm::EclipseStateConstPtr eclipseState,
               Opm::WellConstPtr& well,
               PhaseUsage uses,                  /* phases present     */
               BlackoilPhases::PhaseIndex phase, /* oil, water or gas  */
               WellType type,                    /* prod. or inj.      */
               char aggregation,                 /* rate or total      */
               std::string unit)
        // save these for when we update the value in a timestep
        : eclipseState_(eclipseState)
        , well_(well)
        , phaseUses_(uses)
        , phaseIdx_(phase)
    {
        // producers can be seen as negative injectors
        if (type == INJECTOR)
            sign_ = +1.0;
        else
            sign_ = -1.0;
        ertHandle_ = ecl_sum_add_var(summary.ertHandle(),
                                     varName_(phase,
                                              type,
                                              aggregation).c_str(),
                                     well_->name().c_str(),
                                     /*num=*/ 0,
                                     unit.c_str(),
                                     /*defaultValue=*/ 0.);
    }

public:
    /// Retrieve the value which the monitor is supposed to write to the summary file
    /// according to the state of the well.
    virtual double retrieveValue(const int writeStepIdx,
                                 const SimulatorTimerInterface& timer,
                                 const WellState& wellState,
                                 const std::map<std::string, int>& nameToIdxMap) = 0;

    smspec_node_type *ertHandle() const
    { return ertHandle_; }

protected:
    void updateTimeStepWellIndex_(const std::map<std::string, int>& nameToIdxMap)
    {
        const std::string& wellName = well_->name();

        const auto wellIdxIt = nameToIdxMap.find(wellName);
        if (wellIdxIt == nameToIdxMap.end()) {
            timeStepWellIdx_ = -1;
            flatIdx_ = -1;
            return;
        }

        timeStepWellIdx_ = wellIdxIt->second;
        flatIdx_ = timeStepWellIdx_*phaseUses_.num_phases + phaseUses_.phase_pos[phaseIdx_];
    }

    // return m^3/s of injected or produced fluid
    double rate(const WellState& wellState)
    {
        double value = 0;
        if (wellState.wellRates().size() > 0) {
            assert(int(wellState.wellRates().size()) > flatIdx_);
            value = sign_ * wellState.wellRates()[flatIdx_];
        }
        return value;
    }

    double bhp(const WellState& wellState)
    {
        if (wellState.bhp().size() > 0) {
            // Note that 'flatIdx_' is used here even though it is meant
            // to give a (well,phase) pair.
            const int numPhases = wellState.wellRates().size() / wellState.bhp().size();

            return wellState.bhp()[flatIdx_/numPhases];
        }
        return 0.0;
    }

    /// Get the index associated a well name
    int wellIndex_(Opm::EclipseStateConstPtr eclipseState)
    {
        const Opm::ScheduleConstPtr schedule = eclipseState->getSchedule();

        const std::string& wellName = well_->name();
        const auto& wells = schedule->getWells();
        for (size_t wellIdx = 0; wellIdx < wells.size(); ++wellIdx) {
            if (wells[wellIdx]->name() == wellName) {
                return wellIdx;
            }
        }

        OPM_THROW(std::runtime_error,
                  "Well '" << wellName << "' is not present in deck");
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

    smspec_node_type *ertHandle_;

    Opm::EclipseStateConstPtr eclipseState_;
    Opm::WellConstPtr well_;

    PhaseUsage phaseUses_;
    BlackoilPhases::PhaseIndex phaseIdx_;

    int timeStepWellIdx_;

    /// index into a (flattened) wellsOfTimeStep*phases matrix
    int flatIdx_;

    /// natural sign of the rate
    double sign_;
};

/// Monitors the rate given by a well.
class WellRate : public WellReport
{
public:
    WellRate(const Summary& summary,
             Opm::EclipseStateConstPtr eclipseState,
             Opm::WellConstPtr well,
             PhaseUsage uses,
             BlackoilPhases::PhaseIndex phase,
             WellType type)
        : WellReport(summary,
                     eclipseState,
                     well,
                     uses,
                     phase,
                     type,
                     'R',
                     "SM3/DAY" /* surf. cub. m. per day */)
    { }

    virtual double retrieveValue(const int /* writeStepIdx */,
                                 const SimulatorTimerInterface& timer,
                                 const WellState& wellState,
                                 const std::map<std::string, int>& wellNameToIdxMap)
    {
        // find the index for the quantity in the wellState
        this->updateTimeStepWellIndex_(wellNameToIdxMap);
        if (this->flatIdx_ < 0) {
            // well not active in current time step
            return 0.0;
        }

        if (well_->getStatus(timer.reportStepNum()) == WellCommon::SHUT) {
            // well is shut in the current time step
            return 0.0;
        }

        // TODO: Why only positive rates?
        using namespace Opm::unit;
        return convert::to(std::max(0., rate(wellState)),
                           cubic(meter)/day);
    }
};

/// Monitors the total production in a well.
class WellTotal : public WellReport
{
public:
    WellTotal(const Summary& summary,
              Opm::EclipseStateConstPtr eclipseState,
              Opm::WellConstPtr well,
              PhaseUsage uses,
              BlackoilPhases::PhaseIndex phase,
              WellType type)
        : WellReport(summary,
                     eclipseState,
                     well,
                     uses,
                     phase,
                     type,
                     'T',
                     "SM3" /* surface cubic meter */ )
          // nothing produced when the reporting starts
        , total_(0.)
    { }

    virtual double retrieveValue(const int writeStepIdx,
                                 const SimulatorTimerInterface& timer,
                                 const WellState& wellState,
                                 const std::map<std::string, int>& wellNameToIdxMap)
    {
        if (writeStepIdx == 0) {
            // We are at the initial state.
            // No step has been taken yet.
            return 0.0;
        }

        if (well_->getStatus(timer.reportStepNum()) == WellCommon::SHUT) {
            // well is shut in the current time step
            return 0.0;
        }

        // find the index for the quantity in the wellState
        this->updateTimeStepWellIndex_(wellNameToIdxMap);
        if (this->flatIdx_ < 0) {
            // well not active in current time step
            return 0.0;
        }

        // due to using an Euler method as time integration scheme, the well rate is the
        // average for the time step. For more complicated time stepping schemes, the
        // integral of the rate is not simply multiplying two numbers...
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
            Opm::WellConstPtr well,
            PhaseUsage uses,
            BlackoilPhases::PhaseIndex phase,
            WellType type)
        : WellReport(summary,
                     eclipseState,
                     well,
                     uses,
                     phase,
                     type,
                     'B',
                     "Pascal")
    { }

    virtual double retrieveValue(const int /* writeStepIdx */,
                                 const SimulatorTimerInterface& timer,
                                 const WellState& wellState,
                                 const std::map<std::string, int>& wellNameToIdxMap)
    {
        // find the index for the quantity in the wellState
        this->updateTimeStepWellIndex_(wellNameToIdxMap);
        if (this->flatIdx_ < 0) {
            // well not active in current time step
            return 0.0;
        }
        if (well_->getStatus(timer.reportStepNum()) == WellCommon::SHUT) {
            // well is shut in the current time step
            return 0.0;
        }

        return bhp(wellState);
    }
};

// no inline implementation of this since it depends on the
// WellReport type being completed first
void Summary::writeTimeStep(int writeStepIdx,
                            const SimulatorTimerInterface& timer,
                            const WellState& wellState)
{
    // create a name -> well index map
    const Opm::ScheduleConstPtr schedule = eclipseState_->getSchedule();
    const auto& timeStepWells = schedule->getWells(timer.reportStepNum());
    std::map<std::string, int> wellNameToIdxMap;
    int openWellIdx = 0;
    for (size_t tsWellIdx = 0; tsWellIdx < timeStepWells.size(); ++tsWellIdx) {
        if (timeStepWells[tsWellIdx]->getStatus(timer.reportStepNum()) != WellCommon::SHUT ) {
            wellNameToIdxMap[timeStepWells[tsWellIdx]->name()] = openWellIdx;
            openWellIdx++;
        }
    }

    // internal view; do not move this code out of Summary!
    SummaryTimeStep tstep(*this, writeStepIdx, timer);
    // write all the variables
    for (auto varIt = summaryReportVars_.begin(); varIt != summaryReportVars_.end(); ++varIt) {
        ecl_sum_tstep_iset(tstep.ertHandle(),
                           smspec_node_get_params_index((*varIt)->ertHandle()),
                           (*varIt)->retrieveValue(writeStepIdx, timer, wellState, wellNameToIdxMap));
    }

    // write the summary file to disk
    ecl_sum_fwrite(ertHandle());
}

void Summary::addAllWells(Opm::EclipseStateConstPtr eclipseState,
                          const PhaseUsage& uses)
{
    eclipseState_ = eclipseState;
    // TODO: Only create report variables that are requested with keywords
    // (e.g. "WOPR") in the input files, and only for those wells that are
    // mentioned in those keywords
    Opm::ScheduleConstPtr schedule = eclipseState->getSchedule();
    const auto& wells = schedule->getWells();
    const int numWells = schedule->numWells();
    for (int phaseIdx = 0; phaseIdx != BlackoilPhases::MaxNumPhases; ++phaseIdx) {
        const BlackoilPhases::PhaseIndex ertPhaseIdx =
            static_cast <BlackoilPhases::PhaseIndex>(phaseIdx);
        // don't bother with reporting for phases that aren't there
        if (!uses.phase_used[phaseIdx]) {
            continue;
        }
        size_t numWellTypes = sizeof(WELL_TYPES) / sizeof(WELL_TYPES[0]);
        for (size_t wellTypeIdx = 0; wellTypeIdx < numWellTypes; ++wellTypeIdx) {
            const WellType wellType = WELL_TYPES[wellTypeIdx];
            for (int wellIdx = 0; wellIdx != numWells; ++wellIdx) {
                // W{O,G,W}{I,P}R
                addWell(std::unique_ptr <WellReport>(
                            new WellRate(*this,
                                         eclipseState,
                                         wells[wellIdx],
                                         uses,
                                         ertPhaseIdx,
                                         wellType)));
                // W{O,G,W}{I,P}T
                addWell(std::unique_ptr <WellReport>(
                            new WellTotal(*this,
                                          eclipseState,
                                          wells[wellIdx],
                                          uses,
                                          ertPhaseIdx,
                                          wellType)));
            }
        }
    }

    // Add BHP monitors
    for (int wellIdx = 0; wellIdx != numWells; ++wellIdx) {
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
                                wells[wellIdx],
                                uses,
                                ertPhaseIdx,
                                WELL_TYPES[0])));
    }
}
} // end namespace EclipseWriterDetails



/**
 * Convert opm-core WellType and InjectorType to eclipse welltype
 */
int EclipseWriter::eclipseWellTypeMask(WellType wellType, WellInjector::TypeEnum injectorType)
{
  int ert_well_type = IWEL_UNDOCUMENTED_ZERO;

  if (PRODUCER == wellType) {
      ert_well_type = IWEL_PRODUCER;
  } else if (INJECTOR == wellType) {
      switch (injectorType) {
        case WellInjector::WATER:
          ert_well_type = IWEL_WATER_INJECTOR;
          break;
        case WellInjector::GAS:
          ert_well_type = IWEL_GAS_INJECTOR;
          break;
        case WellInjector::OIL :
          ert_well_type = IWEL_OIL_INJECTOR;
          break;
        default:
          ert_well_type = IWEL_UNDOCUMENTED_ZERO;
      }
  }

  return ert_well_type;
}


/**
 * Convert opm-core WellStatus to eclipse format: > 0 open, <= 0 shut
 */
int EclipseWriter::eclipseWellStatusMask(WellCommon::StatusEnum wellStatus)
{
  int well_status = 0;

  if (wellStatus == WellCommon::OPEN) {
    well_status = 1;
  }
  return well_status;
}



/**
 * Convert opm-core UnitType to eclipse format: ert_ecl_unit_enum
 */
ert_ecl_unit_enum EclipseWriter::convertUnitTypeErtEclUnitEnum(UnitSystem::UnitType unit)
{
    ert_ecl_unit_enum ecl_type;
    switch (unit) {
      case(UnitSystem::UNIT_TYPE_METRIC):
          ecl_type = ERT_ECL_METRIC_UNITS;
          break;
      case(UnitSystem::UNIT_TYPE_FIELD)          :
          ecl_type = ERT_ECL_FIELD_UNITS;
          break;
      case(UnitSystem::UNIT_TYPE_LAB):
          ecl_type = ERT_ECL_LAB_UNITS;
          break;
      default:
          break;
    };

    return ecl_type;
}


void EclipseWriter::writeInit(const SimulatorTimerInterface &timer)
{
    // if we don't want to write anything, this method becomes a
    // no-op...
    if (!enableOutput_) {
        return;
    }

    writeStepIdx_ = 0;

    EclipseWriterDetails::Init fortio(outputDir_, baseName_, /*stepIdx=*/0);
    fortio.writeHeader(numCells_,
                       compressedToCartesianCellIdx_,
                       timer,
                       eclipseState_,
                       phaseUsage_);

    if (eclipseState_->hasDoubleGridProperty("PERMX")) {
        auto data = eclipseState_->getDoubleGridProperty("PERMX")->getData();
        EclipseWriterDetails::convertFromSiTo(data, Opm::prefix::milli * Opm::unit::darcy);
        fortio.writeKeyword("PERMX", data);
    }
    if (eclipseState_->hasDoubleGridProperty("PERMY")) {
        auto data = eclipseState_->getDoubleGridProperty("PERMY")->getData();
        EclipseWriterDetails::convertFromSiTo(data, Opm::prefix::milli * Opm::unit::darcy);
        fortio.writeKeyword("PERMY", data);
    }
    if (eclipseState_->hasDoubleGridProperty("PERMZ")) {
        auto data = eclipseState_->getDoubleGridProperty("PERMZ")->getData();
        EclipseWriterDetails::convertFromSiTo(data, Opm::prefix::milli * Opm::unit::darcy);
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

// implementation of the writeTimeStep method
void EclipseWriter::writeTimeStep(const SimulatorTimerInterface& timer,
                                  const SimulatorState& reservoirState,
                                  const WellState& wellState)
{
    // if we don't want to write anything, this method becomes a
    // no-op...
    if (!enableOutput_) {
        return;
    }

    // respected the output_interval parameter
    if (writeStepIdx_ % outputInterval_ != 0) {
        return;
    }

    const size_t ncwmax                 = eclipseState_->getSchedule()->getMaxNumCompletionsForWells(timer.reportStepNum());
    const size_t numWells               = eclipseState_->getSchedule()->numWells(timer.reportStepNum());
    std::vector<WellConstPtr> wells_ptr = eclipseState_->getSchedule()->getWells(timer.reportStepNum());

    std::vector<const char*> zwell_data( numWells * Opm::EclipseWriterDetails::Restart::NZWELZ , "");
    std::vector<int>         iwell_data( numWells * Opm::EclipseWriterDetails::Restart::NIWELZ , 0 );
    std::vector<int>         icon_data( numWells * ncwmax * Opm::EclipseWriterDetails::Restart::NICONZ , 0 );

    EclipseWriterDetails::Restart restartHandle(outputDir_, baseName_, writeStepIdx_);

    for (size_t iwell = 0; iwell < wells_ptr.size(); ++iwell) {
        WellConstPtr well = wells_ptr[iwell];
        {
            size_t wellIwelOffset = Opm::EclipseWriterDetails::Restart::NIWELZ * iwell;
            restartHandle.addRestartFileIwelData(iwell_data, timer.reportStepNum(), well , wellIwelOffset);
        }
        {
            size_t wellIconOffset = ncwmax * Opm::EclipseWriterDetails::Restart::NICONZ * iwell;
            restartHandle.addRestartFileIconData(icon_data,  well->getCompletions( timer.reportStepNum() ), wellIconOffset);
        }
        zwell_data[ iwell * Opm::EclipseWriterDetails::Restart::NZWELZ ] = well->name().c_str();
    }

    {
        ecl_rsthead_type rsthead_data = {};
        rsthead_data.sim_time   = timer.currentPosixTime();
        rsthead_data.nactive    = numCells_;
        rsthead_data.nx         = cartesianSize_[0];
        rsthead_data.ny         = cartesianSize_[1];
        rsthead_data.nz         = cartesianSize_[2];
        rsthead_data.nwells     = numWells;
        rsthead_data.niwelz     = EclipseWriterDetails::Restart::NIWELZ;
        rsthead_data.nzwelz     = EclipseWriterDetails::Restart::NZWELZ;
        rsthead_data.niconz     = EclipseWriterDetails::Restart::NICONZ;
        rsthead_data.ncwmax     = ncwmax;
        rsthead_data.phase_sum  = Opm::EclipseWriterDetails::ertPhaseMask(phaseUsage_);
        rsthead_data.sim_days   = Opm::unit::convert::to(timer.simulationTimeElapsed(), Opm::unit::day); //data for doubhead

        restartHandle.writeHeader(timer,
                                  writeStepIdx_,
                                  &rsthead_data);
    }


    restartHandle.add_kw(EclipseWriterDetails::Keyword<int>(IWEL_KW, iwell_data));
    restartHandle.add_kw(EclipseWriterDetails::Keyword<const char *>(ZWEL_KW, zwell_data));
    restartHandle.add_kw(EclipseWriterDetails::Keyword<int>(ICON_KW, icon_data));

    EclipseWriterDetails::Solution sol(restartHandle);

    // write out the pressure of the reference phase (whatever phase that is...). this is
    // not the most performant solution thinkable, but this is also not in the most
    // performance critical code path!
    //
    // Also, we want to use the same units as the deck for pressure output, i.e. we have
    // to mutliate our nice SI pressures by the inverse of the conversion factor of deck
    // to SI pressure units...
    std::vector<double> pressure = reservoirState.pressure();
    EclipseWriterDetails::convertFromSiTo(pressure, deckToSiPressure_);
    EclipseWriterDetails::restrictAndReorderToActiveCells(pressure, gridToEclipseIdx_.size(), gridToEclipseIdx_.data());

    sol.add(EclipseWriterDetails::Keyword<float>("PRESSURE", pressure));

    std::vector<double> saturation_water;
    std::vector<double> saturation_gas;


    if (phaseUsage_.phase_used[BlackoilPhases::Aqua]) {
        saturation_water = reservoirState.saturation();
        EclipseWriterDetails::extractFromStripedData(saturation_water,
                                                     /*offset=*/phaseUsage_.phase_pos[BlackoilPhases::Aqua],
                                                     /*stride=*/phaseUsage_.num_phases);
        EclipseWriterDetails::restrictAndReorderToActiveCells(saturation_water, gridToEclipseIdx_.size(), gridToEclipseIdx_.data());
        sol.add(EclipseWriterDetails::Keyword<float>(EclipseWriterDetails::saturationKeywordNames[BlackoilPhases::PhaseIndex::Aqua], saturation_water));
    }


    if (phaseUsage_.phase_used[BlackoilPhases::Vapour]) {
        saturation_gas = reservoirState.saturation();
        EclipseWriterDetails::extractFromStripedData(saturation_gas,
                                                     /*offset=*/phaseUsage_.phase_pos[BlackoilPhases::Vapour],
                                                     /*stride=*/phaseUsage_.num_phases);
        EclipseWriterDetails::restrictAndReorderToActiveCells(saturation_gas, gridToEclipseIdx_.size(), gridToEclipseIdx_.data());
        sol.add(EclipseWriterDetails::Keyword<float>(EclipseWriterDetails::saturationKeywordNames[BlackoilPhases::PhaseIndex::Vapour], saturation_gas));
    }



    //Write RFT data for current timestep to RFT file
    std::shared_ptr<EclipseWriterDetails::EclipseWriteRFTHandler> eclipseWriteRFTHandler = std::make_shared<EclipseWriterDetails::EclipseWriteRFTHandler>(
                                                                                                                      compressedToCartesianCellIdx_,
                                                                                                                      numCells_,
                                                                                                                      eclipseState_->getEclipseGrid()->getCartesianSize());


    char * rft_filename = ecl_util_alloc_filename(outputDir_.c_str(),
                                                  baseName_.c_str(),
                                                  ECL_RFT_FILE,
                                                  false,
                                                  0);

    std::shared_ptr<const UnitSystem> unitsystem = eclipseState_->getDeckUnitSystem();
    ert_ecl_unit_enum ecl_unit = convertUnitTypeErtEclUnitEnum(unitsystem->getType());

    std::vector<WellConstPtr> wells = eclipseState_->getSchedule()->getWells(timer.currentStepNum());


    eclipseWriteRFTHandler->writeTimeStep(rft_filename,
                                          ecl_unit,
                                          timer,
                                          wells,
                                          eclipseState_->getEclipseGrid(),
                                          pressure,
                                          saturation_water,
                                          saturation_gas);



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
    summary_->writeTimeStep(writeStepIdx_, timer, wellState);

    ++writeStepIdx_;
}


EclipseWriter::EclipseWriter(const parameter::ParameterGroup& params,
                             Opm::EclipseStateConstPtr eclipseState,
                             const Opm::PhaseUsage &phaseUsage,
                             int numCells,
                             const int* compressedToCartesianCellIdx)
    : eclipseState_(eclipseState)
    , numCells_(numCells)
    , compressedToCartesianCellIdx_(compressedToCartesianCellIdx)
    , gridToEclipseIdx_(numCells, int(-1) )
    , phaseUsage_(phaseUsage)
{
    const auto eclGrid = eclipseState->getEclipseGrid();
    cartesianSize_[0] = eclGrid->getNX();
    cartesianSize_[1] = eclGrid->getNY();
    cartesianSize_[2] = eclGrid->getNZ();

    if( compressedToCartesianCellIdx ) {
        // if compressedToCartesianCellIdx available then
        // compute mapping to eclipse order
        std::map< int , int > indexMap;
        for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            int cartesianCellIdx = compressedToCartesianCellIdx[cellIdx];
            indexMap[ cartesianCellIdx ] = cellIdx;
        }

        int idx = 0;
        for( auto it = indexMap.begin(), end = indexMap.end(); it != end; ++it ) {
            gridToEclipseIdx_[ idx++ ] = (*it).second;
        }
    }
    else {
        // if not compressedToCartesianCellIdx was given use identity
        for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            gridToEclipseIdx_[ cellIdx ] = cellIdx;
        }
    }


    // factor from the pressure values given in the deck to Pascals
    deckToSiPressure_ =
        eclipseState->getDeckUnitSystem()->parse("Pressure")->getSIScaling();

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
    writeStepIdx_ = 0;

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
