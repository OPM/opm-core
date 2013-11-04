/*
  Copyright 2013 Andreas Lauser

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

#ifndef OPM_BLACKOIL_ECLIPSE_OUTPUT_WRITER_HPP
#define OPM_BLACKOIL_ECLIPSE_OUTPUT_WRITER_HPP

#include <opm/core/grid.h>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <list>
#include <string>
#include <array>
#include <vector>

#ifdef HAVE_ERT
#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_util.h>
#endif

#include <string>

namespace Opm {
/*!
 * \brief A class to write the reservoir state and the well state of a
 *        blackoil simulation to disk using the Eclipse binary format.
 *
 * This class only writes files if the 'write_output' parameter is set
 * to 1. It needs the ERT libraries to write to disk, so if the
 * 'write_output' parameter is set but ERT is not available, all
 * methods throw a std::runtime_error.
 */
class BlackoilEclipseOutputWriter
{
public:
    /*!
     * \brief Sets the common attributes required to write eclipse
     * binary files using ERT.
     */
    BlackoilEclipseOutputWriter(const EclipseGridParser& eclipseParser,
                                const UnstructuredGrid& grid,
                                const std::string &outputDir,
                                const std::string &baseName)
    : eclipseParser_(eclipseParser)
    , grid_(grid)
    , outputDir_(outputDir)
    , baseName_(baseName)
    {
#if HAVE_ERT
        sumWriter_ = 0;
#endif
    }

    ~BlackoilEclipseOutputWriter()
    {
#if HAVE_ERT
        if (sumWriter_) {
            // clean after ourselfs
            ecl_sum_free(sumWriter_);
        }
#endif
    }

    /*!
     * \brief Write the static eclipse data (grid, PVT curves, etc) to disk
     */
    void writeInitFile(const SimulatorTimer &timer);

    /*!
     * \brief Write a blackoil reservoir state to disk for later inspection with
     *        visualization tools like ResInsight
     *
     * \param[in] reservoirState The thermodynamic state of the reservoir
     */
    void writeReservoirState(const BlackoilState& reservoirState,
                             const SimulatorTimer& timer);

    /*!
     * \brief Write a well state to disk for later inspection with
     *        visualization tools
     *
     * \param[in] wellState The production/injection data for all wells
     */
    void writeWellState(const WellState& wellState, const SimulatorTimer& timer);

private:
    const EclipseGridParser& eclipseParser_;
    const UnstructuredGrid& grid_;
    std::string outputDir_;
    std::string baseName_;

    time_t startTime_;

#if HAVE_ERT
    void writeSummaryHeaderFile_(const SimulatorTimer &timer);
    void writeGridInitFile_(const SimulatorTimer &timer);

    ecl_grid_type* newEclGrid_();

    // keyword handles per well each
    std::vector<smspec_node_type*> woprSmspec_;
    std::vector<smspec_node_type*> woptSmspec_;
    std::vector<smspec_node_type*> wgprSmspec_;
    std::vector<smspec_node_type*> wgptSmspec_;
    std::vector<smspec_node_type*> wwirSmspec_;
    std::vector<smspec_node_type*> wwitSmspec_;
    std::vector<smspec_node_type*> wgirSmspec_;
    std::vector<smspec_node_type*> wgitSmspec_;

    ecl_sum_type* sumWriter_;

    std::vector<std::array<double, /*numPhases=*/3> > accumulatedProducedFluids_;
    std::vector<std::array<double, /*numPhases=*/3> > accumulatedInjectedFluids_;
#endif
};
} // namespace Opm


#endif // OPM_BLACKOIL_ECLIPSE_OUTPUT_WRITER_HPP
