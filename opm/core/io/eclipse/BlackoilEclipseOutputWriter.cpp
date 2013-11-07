 /*
  Copyright (c) 2013 Andreas Lauser
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
#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "BlackoilEclipseOutputWriter.hpp"

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/parameters/Parameter.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/wells.h> // WellType

#include <boost/algorithm/string/case_conv.hpp> // to_upper_copy
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp> // path

#include <memory>     // unique_ptr
#include <utility>    // move

using namespace Opm;
using namespace Opm::parameter;

#ifdef HAVE_ERT
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

// namespace start here since we don't want the ERT headers in it;
// this means we must also do it in the #else section at the bottom
namespace Opm {

namespace internal {

/// Smart pointer/handle class for ERT opaque types, such as ecl_kw_type*.
///
/// \tparam T Type of handle being wrapper
template <typename T>
struct EclipseHandle {
    /// Instead of inheriting std::unique_ptr and letting the compiler
    /// provide a default implementation which calls the base class, we
    /// define the move constructor and assignment operator ourselves
    /// and call and aggregated pointer, because of bugs in GCC 4.4
    EclipseHandle <T> (EclipseHandle <T>&& rhs)
        : h_ (std::move (rhs.h_)) { }

    EclipseHandle <T>& operator= (EclipseHandle <T>&& rhs) {
        h_ = std::move (rhs.h_);
        return *this;
    }

    /// Prevent GCC 4.4 from the urge of generating a copy constructor
    EclipseHandle (const EclipseHandle&) = delete;
    EclipseHandle <T>& operator= (const EclipseHandle <T>&) = delete;

    /// Construct a new smart handle based on the returned value of
    /// an allocation function and the corresponding destroyer function.
    EclipseHandle <T> (T* t, void (*destroy)(T*))
        : h_ (t, destroy) { }

    /// Convenience operator that lets us use this type as if
    /// it was a handle directly.
    operator T* () const { return h_.get (); }

private:
    std::unique_ptr <T, void (*)(T*) throw()> h_; // handle
};

/**
 * Eclipse "keyword" (i.e. named data) for a vector. (This class is
 * different from EclKW in the constructors it provide).
 */
template <typename T>
struct EclipseKeyword : public EclipseHandle <ecl_kw_type> {
    EclipseKeyword (const std::string& name,    /// identification
                    const std::vector<T>& data, /// array holding values
                    const int offset = 0,      /// distance to first
                    const int stride = 1)      /// distance between each

        // allocate handle and put in smart pointer base class
        : EclipseHandle <ecl_kw_type> (
              ecl_kw_alloc (name.c_str(), data.size (), type ()),
              ecl_kw_free) {
        copyData (data, offset, stride);
    }

    /// Convenience constructor that gets the set of data
    /// from the samely named item in the parser
    EclipseKeyword (const std::string& name,
                    const EclipseGridParser& parser)
        // allocate handle and put in smart pointer base class
        : EclipseHandle <ecl_kw_type> (
              ecl_kw_alloc (name.c_str(),
                            parser.getValue <T> (name).size (),
                            type ()),
              ecl_kw_free) {
        copyData (parser.getValue <T> (name), 0, 1);
    }

    /// Constructor for optional fields
    EclipseKeyword (const std::string& name)
        : EclipseHandle <ecl_kw_type> (0, ecl_kw_free) {
        static_cast<void> (name);
    }

    // GCC 4.4 doesn't generate these constructors for us; provide the
    // default implementation explicitly here instead
    EclipseKeyword (EclipseKeyword&& rhs)
        : EclipseHandle <ecl_kw_type> (std::move (rhs)) { }
    EclipseKeyword& operator= (EclipseKeyword&& rhs) {
        EclipseHandle <ecl_kw_type>::operator= (std::move(rhs));
        return *this;
    }
    EclipseKeyword (const EclipseKeyword&) = delete;
    EclipseKeyword& operator= (const EclipseKeyword&) = delete;

private:
    /// Map the C++ data type (given by T) to an Eclipse type enum
    static ecl_type_enum type ();

    /// Helper function that is the meat of the constructor
    void copyData (const std::vector <T>& data,
                    const int offset,
                    const int stride) {
        // number of elements we can possibly take from the vector
        const int num = data.size ();

        // range cannot start outside of data set
        assert(offset >= 0 && offset < num);

        // don't jump out of the set when trying to
        assert(stride > 0 && stride < num - offset);

        // fill it with values
        for (int i = 0; i < num; ++i) {
            // access from data store
            const float value = data[i * stride + offset];

            // write into memory represented by handle
            ecl_kw_iset_float(*this, i, value);
        }
    }
};

// specializations for known keyword types
template <> ecl_type_enum EclipseKeyword<int   >::type () { return ECL_INT_TYPE   ; }
template <> ecl_type_enum EclipseKeyword<double>::type () { return ECL_FLOAT_TYPE; }

/**
 * Extract the current time from a timer object into the C type used by ERT.
 */
static time_t current (const SimulatorTimer& timer) {
    tm t = boost::posix_time::to_tm (timer.currentDateTime());
    return mktime(&t);
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
        : EclipseHandle <const char> (
              ecl_util_alloc_filename (outputDir.c_str(),
                                       baseName.c_str(),
                                       type,
                                       false, // formatted?
                                       timer.currentStepNum ()),
              freestr) { }
private:
    /// Facade which allows us to free a const char*
    static void freestr (const char* ptr) {
        ::free (const_cast<char*>(ptr));
    }
};

struct EclipseRestart : public EclipseHandle <ecl_rst_file_type> {
    EclipseRestart (const std::string& outputDir,
                    const std::string& baseName,
                    const SimulatorTimer& timer)
        // notice the poor man's polymorphism of the allocation function
        : EclipseHandle <ecl_rst_file_type> (
              (timer.currentStepNum () > 0 ? ecl_rst_file_open_append
                                           : ecl_rst_file_open_write)(
                  EclipseFileName (outputDir,
                                   baseName,
                                   ECL_UNIFIED_RESTART_FILE,
                                   timer)),
              ecl_rst_file_close) { }

    void writeHeader (const SimulatorTimer& timer,
                      const int phases,
                      const EclipseGridParser parser,
                      const int num_active_cells) {
        const std::vector<int> dim = parser.getSPECGRID ().dimensions;
        ecl_rst_file_fwrite_header (*this,
                                    timer.currentStepNum (),
                                    current (timer),
                                    Opm::unit::convert::to (timer.currentTime (),
                                                            Opm::unit::day),
                                    dim[0],
                                    dim[1],
                                    dim[2],
                                    num_active_cells,
                                    phases);
    }
};

/**
 * The EclipseSolution class wraps the actions that must be done to the
 * restart file while writing solution variables; it is not a handle on
 * its own.
 */
struct EclipseSolution : public EclipseHandle <ecl_rst_file_type> {
    EclipseSolution (EclipseRestart& rst_file)
        : EclipseHandle <ecl_rst_file_type> (start_solution (rst_file),
                                             ecl_rst_file_end_solution) { }

    template <typename T>
    void add (const EclipseKeyword<T>& kw) {
        ecl_rst_file_add_kw (*this, kw);
    }

private:
    /// Helper method to call function *and* return the handle
    static ecl_rst_file_type* start_solution (EclipseRestart& rst_file) {
        ecl_rst_file_start_solution (rst_file);
        return rst_file;
    }
};

/**
 * Representation of an Eclipse grid.
 */
struct EclipseGrid : public EclipseHandle <ecl_grid_type> {
    /// Create a grid based on the keywords available in input file
    static EclipseGrid make (const EclipseGridParser& parser) {
        if (parser.hasField("DXV")) {
            // make sure that the DYV and DZV keywords are present if the
            // DXV keyword is used in the deck...
            assert(parser.hasField("DYV"));
            assert(parser.hasField("DZV"));

            const auto& dxv = parser.getFloatingPointValue("DXV");
            const auto& dyv = parser.getFloatingPointValue("DYV");
            const auto& dzv = parser.getFloatingPointValue("DZV");

            return EclipseGrid (dxv, dyv, dzv);
        }
        else if (parser.hasField("ZCORN")) {
            struct grdecl g = parser.get_grdecl ();

            EclipseKeyword<double> coord_kw   (COORD_KW,  parser);
            EclipseKeyword<double> zcorn_kw   (ZCORN_KW,  parser);
            EclipseKeyword<double> actnum_kw  (ACTNUM_KW, parser);
            EclipseKeyword<double> mapaxes_kw (MAPAXES_KW);

            if (g.mapaxes) {
                mapaxes_kw = std::move (EclipseKeyword<double> (MAPAXES_KW, parser));
            }

            return EclipseGrid (g.dims, zcorn_kw, coord_kw, actnum_kw, mapaxes_kw);
        }
        else {
            OPM_THROW(std::runtime_error,
                  "Can't create an ERT grid (no supported keywords found in deck)");
        }
    }

    /**
     * Save the grid in an .EGRID file.
     */
    void write (const std::string& outputDir,
                 const std::string& baseName,
                 const SimulatorTimer& timer) {
        ecl_grid_fwrite_EGRID (*this,
                               EclipseFileName (outputDir,
                                                baseName,
                                                ECL_EGRID_FILE,
                                                timer));
    }

    // GCC 4.4 doesn't generate these constructors for us; provide the
    // default implementation explicitly here instead
    EclipseGrid (EclipseGrid&& rhs)
        : EclipseHandle <ecl_grid_type> (std::move (rhs)) { }
    EclipseGrid& operator= (EclipseGrid&& rhs) {
        EclipseHandle <ecl_grid_type>::operator= (std::move(rhs));
        return *this;
    }
    EclipseGrid (const EclipseGrid&) = delete;
    EclipseGrid& operator= (const EclipseGrid&) = delete;

private:
    // each of these cases could have been their respective subclass,
    // but there is not any polymorphism on each of these grid types
    // once we have the handle

    // setup smart pointer for Cartesian grid
    EclipseGrid (const std::vector<double>& dxv,
                 const std::vector<double>& dyv,
                 const std::vector<double>& dzv)
        : EclipseHandle <ecl_grid_type> (
              ecl_grid_alloc_dxv_dyv_dzv (dxv.size (),
                                          dyv.size (),
                                          dzv.size (),
                                          &dxv[0],
                                          &dyv[0],
                                          &dzv[0],
                                          NULL),
              ecl_grid_free) { }

    // setup smart pointer for cornerpoint grid
    EclipseGrid (const int dims[],
                 const EclipseKeyword<double>& zcorn,
                 const EclipseKeyword<double>& coord,
                 const EclipseKeyword<double>& actnum,
                 const EclipseKeyword<double>& mapaxes)
        : EclipseHandle <ecl_grid_type> (
              ecl_grid_alloc_GRDECL_kw(dims[0],
                                       dims[1],
                                       dims[2],
                                       zcorn,
                                       coord,
                                       actnum,
                                       mapaxes),
              ecl_grid_free) { }
};

/**
 * Initialization file which contains static properties (such as
 * porosity and permeability) for the simulation field.
 */
struct EclipseInit : public EclipseHandle <fortio_type> {
    // contrary to the grid, the location of the file goes here because
    // there is only one construction method but several write methods
    // (but we need to do a bit of logic before we can call the actual
    // constructor, so we'll have to do with a static wrapper)
    static EclipseInit make (const std::string& outputDir,
                              const std::string& baseName,
                              const SimulatorTimer& timer) {
        EclipseFileName initFileName (outputDir,
                                      baseName,
                                      ECL_INIT_FILE,
                                      timer);
        bool fmt_file;
        if (!ecl_util_fmt_file(initFileName, &fmt_file)) {
            OPM_THROW(std::runtime_error,
                      "Could not determine formatted/unformatted status of file:" << initFileName << " non-standard name?" << std::endl);
        }
        return EclipseInit (initFileName, fmt_file);
    }

    void writeHeader (const EclipseGrid& grid,
                       const SimulatorTimer& timer,
                       const EclipseGridParser& parser,
                       const int phases) {
        EclipseKeyword<double> poro (PORO_KW, parser);
        ecl_init_file_fwrite_header (*this,
                                     grid,
                                     poro,
                                     phases,
                                     current (timer));
    }

    template <typename T>
    void writeKeyword (const std::string& keyword,
                        const EclipseGridParser& parser) {
        EclipseKeyword <T> kw (keyword, parser);
        ecl_kw_fwrite (kw, *this);
    }

    // GCC 4.4 doesn't generate these constructors for us; provide the
    // default implementation explicitly here instead
    EclipseInit (EclipseInit&& rhs)
        : EclipseHandle <fortio_type> (std::move (rhs)) { }
    EclipseInit& operator= (EclipseInit& rhs) {
        EclipseHandle <fortio_type>::operator= (std::move(rhs));
        return *this;
    }
    EclipseInit (const EclipseInit&) = delete;
    EclipseInit& operator= (const EclipseInit&) = delete;
private:
    EclipseInit (const EclipseFileName& fname, const bool formatted)
        : EclipseHandle <fortio_type> (
              fortio_open_writer (fname, formatted, ECL_ENDIAN_FLIP),
              fortio_fclose) { }
};

// forward decl. of mutually dependent type
struct EclipseWellReport;

struct EclipseSummary : public EclipseHandle <ecl_sum_type> {
    EclipseSummary (const std::string& outputDir,
                    const std::string& baseName,
                    const SimulatorTimer& timer,
                    const EclipseGridParser parser)
        : EclipseHandle <ecl_sum_type> (
              alloc_writer (outputDir, baseName, timer, parser),
              ecl_sum_free) { }

    typedef std::unique_ptr <EclipseWellReport> var_t;
    typedef std::vector <var_t> vars_t;

    EclipseSummary& add (var_t var) {
        vars_.push_back (std::move (var));
        return *this;
    }

    // no inline implementation of this since it depends on the
    // EclipseWellReport type being completed first
    void writeTimeStep (const SimulatorTimer& timer,
                         const WellState& wellState);

private:
    vars_t vars_;

    /// Make sure a new timestep is flushed to the summary section
    struct EclipseTimeStep : public EclipseHandle <ecl_sum_tstep_type> {
        EclipseTimeStep (const EclipseSummary& sum,
                         const SimulatorTimer& timer)
            : EclipseHandle <ecl_sum_tstep_type> (
                  ecl_sum_add_tstep (sum,
                                     timer.currentStepNum () + 1,
                                     // currentTime is always relative to start
                                     Opm::unit::convert::to (timer.currentTime (),
                                                             Opm::unit::day)),
                  ecl_sum_tstep_free)
            , sum_ (sum) { }

        ~EclipseTimeStep () {
            ecl_sum_fwrite (sum_);
        }
    private:
        ecl_sum_type* sum_;
    };

    /// Helper routine that lets us use local variables to hold
    /// intermediate results while filling out the allocations function's
    /// argument list.
    static ecl_sum_type* alloc_writer (const std::string& outputDir,
                                        const std::string& baseName,
                                        const SimulatorTimer& timer,
                                        const EclipseGridParser& parser) {
        boost::filesystem::path casePath (outputDir);
        casePath /= boost::to_upper_copy (baseName);
        const std::vector<int>& dim = parser.getSPECGRID().dimensions;
        return ecl_sum_alloc_writer (casePath.string ().c_str (),
                                     false, /* formatted   */
                                     true,  /* unified     */
                                     ":",    /* join string */
                                     current (timer),
                                     dim[0],
                                     dim[1],
                                     dim[2]);
    }
};


/**
 * Summary variable that reports a characteristics of a well.
 */
struct EclipseWellReport : public EclipseHandle <smspec_node_type> {
protected:
    EclipseWellReport (const EclipseSummary& summary,    /* section to add to  */
                       const EclipseGridParser& parser,  /* well names         */
                       int whichWell,                    /* index of well line */
                       BlackoilPhases::PhaseIndex phase, /* oil, water or gas  */
                       WellType type,                    /* prod. or inj.      */
                       char aggregation,                 /* rate or total      */
                       std::string unit)
        : EclipseHandle <smspec_node_type> (
              ecl_sum_add_var (summary,
                               varName (phase,
                                        type,
                                        aggregation).c_str (),
                               wellName (parser, whichWell).c_str (),
                               /* num = */ 0,
                               unit.c_str(),
                               /* defaultValue = */ 0.),
              smspec_node_free)
        // save these for when we update the value in a timestep
        , index_ (whichWell * BlackoilPhases::MaxNumPhases + phase)

        // injectors can be seen as negative producers
        , sign_ (type == INJECTOR ? -1. : +1.) { }

public:
    /// Allows us to pass this type to ecl_sum_tstep_iset
    operator int () {
        return smspec_node_get_params_index (*this);
    }

    /// Update the monitor according to the new state of the well, and
    /// get the reported value. Note: Only call this once for each timestep.
    virtual double update (const SimulatorTimer& timer,
                             const WellState& wellState) = 0;

private:
    /// index into a (flattened) wells*phases matrix
    const int index_;

    /// natural sign of the rate
    const double sign_;

    /// Get the name associated with this well
    std::string wellName (const EclipseGridParser& parser,
                          int whichWell) {
        return parser.getWELSPECS().welspecs[whichWell].name_;
    }

    /// Compose the name of the summary variable, e.g. "WOPR" for
    /// well oil production rate.
    std::string varName (BlackoilPhases::PhaseIndex phase,
                         WellType type,
                         char aggregation) {
        std::string name;
        name += 'W'; // well
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
        return name;
    }
protected:
    double rate (const WellState& wellState) {
        // convert m^3/s of injected fluid to m^3/d of produced fluid
        const double convFactor = Opm::unit::convert::to (1., Opm::unit::day);
        // TODO: Correct to flip sign to get positive values?
        const double value = sign_ * wellState.wellRates () [index_] * convFactor;
        return value;
    }
};

/// Monitors the rate given by a well.
struct EclipseWellRate : public EclipseWellReport {
    EclipseWellRate (const EclipseSummary& summary,
                     const EclipseGridParser& parser,
                     int whichWell,
                     BlackoilPhases::PhaseIndex phase,
                     WellType type)
        : EclipseWellReport (summary,
                             parser,
                             whichWell,
                             phase,
                             type,
                             'R',
                             "SM3/DAY" /* surf. cub. m. per day */ ) { }
    virtual double update (const SimulatorTimer& timer,
                             const WellState& wellState) {
        // TODO: Why only positive rates?
        return std::max (0., rate (wellState));
    }
};

/// Monitors the total production in a well.
struct EclipseWellTotal : public EclipseWellReport {
    EclipseWellTotal (const EclipseSummary& summary,
                      const EclipseGridParser& parser,
                      int whichWell,
                      BlackoilPhases::PhaseIndex phase,
                      WellType type)
        : EclipseWellReport (summary,
                             parser,
                             whichWell,
                             phase,
                             type,
                             'T',
                             "SM3" /* surface cubic meter */ )

        // nothing produced when the reporting starts
        , total_ (0.) { }

    virtual double update (const SimulatorTimer& timer,
                             const WellState& wellState) {
        // TODO: Is the rate average for the timestep, or is in
        // instantaneous (in which case trapezoidal or Simpson integration
        // would probably be better)
        const double intg = timer.currentStepLength () * rate (wellState);
        // add this timesteps production to the total
        total_ += intg;
        // report the new production total
        return total_;
    }

private:
    /// Aggregated value of the course of the simulation
    double total_;
};

inline void
EclipseSummary::writeTimeStep (const SimulatorTimer& timer,
                               const WellState& wellState) {
    EclipseTimeStep tstep (*this, timer);
    // write all the variables
    for (vars_t::iterator v = vars_.begin(); v != vars_.end(); ++v) {
        const double value = (*v)->update (timer, wellState);
        ecl_sum_tstep_iset(tstep, *(*v).get (), value);
    }
}

/// Supported well types. Enumeration doesn't let us get all the members,
/// so we must have an explicit array.
static WellType WELL_TYPES[] = { INJECTOR, PRODUCER };

/// Helper method that can be used in std::transform (must curry the barsa
/// argument)
static double pasToBar (double pressureInPascal) {
    return Opm::unit::convert::to (pressureInPascal, Opm::unit::barsa);
}

} // namespace Opm::internal

using namespace Opm::internal;

void BlackoilEclipseOutputWriter::writeInit(const SimulatorTimer &timer) {
    /* Grid files */
    EclipseGrid ecl_grid = EclipseGrid::make (eclipseParser_);
    ecl_grid.write (outputDir_, baseName_, timer);

    EclipseInit fortio = EclipseInit::make (outputDir_, baseName_, timer);
    fortio.writeHeader (ecl_grid,
                        timer,
                        eclipseParser_,
                        ECL_OIL_PHASE | ECL_GAS_PHASE | ECL_WATER_PHASE);

    fortio.writeKeyword<double> ("PERMX", eclipseParser_);
    fortio.writeKeyword<double> ("PERMY", eclipseParser_);
    fortio.writeKeyword<double> ("PERMZ", eclipseParser_);

    /* Summary files */
    sum_ = std::move (std::unique_ptr <EclipseSummary> (
                          new EclipseSummary (outputDir_,
                                               baseName_,
                                               timer,
                                               eclipseParser_)));

    // TODO: Only create report variables that are requested with keywords
    // (e.g. "WOPR") in the input files, and only for those wells that are
    // mentioned in those keywords
    const int numWells = eclipseParser_.getWELSPECS().welspecs.size();
    for (int whichWell = 0; whichWell != numWells; ++whichWell) {
        for (int phaseCounter = 0;
             phaseCounter != BlackoilPhases::MaxNumPhases;
             ++phaseCounter) {
            const BlackoilPhases::PhaseIndex phase =
                    static_cast <BlackoilPhases::PhaseIndex> (phaseCounter);
            for (size_t typeIndex = 0;
                 typeIndex < sizeof (WELL_TYPES) / sizeof (WELL_TYPES[0]);
                 ++typeIndex) {
                const WellType type = WELL_TYPES[typeIndex];
                // W{O,G,W}{I,P}R
                sum_->add (std::unique_ptr <EclipseWellReport> (
                              new EclipseWellRate (*sum_,
                                                    eclipseParser_,
                                                    whichWell,
                                                    phase,
                                                    type)));
                // W{O,G,W}{I,P}T
                sum_->add (std::unique_ptr <EclipseWellReport> (
                              new EclipseWellTotal (*sum_,
                                                     eclipseParser_,
                                                     whichWell,
                                                     phase,
                                                     type)));
            }
        }
    }

    // flush after all variables are allocated
    ecl_sum_fwrite(*sum_);
}

void BlackoilEclipseOutputWriter::writeTimeStep(
        const SimulatorTimer& timer,
        const BlackoilState& reservoirState,
        const WellState& wellState) {
    // convert the pressures from Pascals to bar because eclipse
    // seems to write bars
    const std::vector<double>& pas = reservoirState.pressure ();
    std::vector<double> bar (pas.size (), 0.);
    std::transform (pas.begin(), pas.end(), bar.begin(), pasToBar);

    // start writing to files
    EclipseRestart rst (outputDir_,
                        baseName_,
                        timer);
    rst.writeHeader (timer,
                     ECL_OIL_PHASE | ECL_GAS_PHASE | ECL_WATER_PHASE,
                     eclipseParser_,
                     pas.size ());
    EclipseSolution sol (rst);

    // write pressure and saturation fields (same as DataMap holds)
    sol.add (EclipseKeyword<double> ("PRESSURE", bar));

    sol.add (EclipseKeyword<double> ("SWAT",
                                      reservoirState.saturation(),
                                      BlackoilPhases::Aqua,
                                      BlackoilPhases::MaxNumPhases));

    sol.add (EclipseKeyword<double> ("SOIL",
                                      reservoirState.saturation(),
                                      BlackoilPhases::Liquid,
                                      BlackoilPhases::MaxNumPhases));

    sol.add (EclipseKeyword<double> ("SGAS",
                                      reservoirState.saturation(),
                                      BlackoilPhases::Vapour,
                                      BlackoilPhases::MaxNumPhases));

    /* Summary variables (well reporting) */
    sum_->writeTimeStep (timer, wellState);
}

#else
namespace Opm {

void BlackoilEclipseOutputWriter::writeInit(const SimulatorTimer &timer) {
    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
}

void BlackoilEclipseOutputWriter::writeTimeStep(
        const SimulatorTimer& timer,
        const BlackoilState& reservoirState,
        const WellState& wellState) {
    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
}

#endif // HAVE_ERT

BlackoilEclipseOutputWriter::BlackoilEclipseOutputWriter (
        const ParameterGroup& params,
        const EclipseGridParser& parser)
    : eclipseParser_ (parser) {

    // get the base name from the name of the deck
    boost::filesystem::path deck (params.get <std::string> ("deck_filename"));
    if (boost::to_upper_copy (deck.extension ().string ()) == ".DATA") {
        baseName_ = deck.stem ().string ();
    }
    else {
        baseName_ = deck.filename ().string ();
    }

    // store in current directory if not explicitly set
    if (params.has ("output_dir")) {
        outputDir_ = params.get <std::string> ("output_dir");
    }
}

} // namespace Opm
