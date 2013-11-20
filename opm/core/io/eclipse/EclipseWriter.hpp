/*
  Copyright (c) 2013 Andreas Lauser
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

#ifndef OPM_ECLIPSE_WRITER_HPP
#define OPM_ECLIPSE_WRITER_HPP

#include <opm/core/io/OutputWriter.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <string>
#include <memory>  // std::unique_ptr

struct UnstructuredGrid;

namespace Opm {

// forward declarations
class BlackoilState;
class EclipseGridParser;
class SimulatorTimer;
class WellState;

namespace parameter { class ParameterGroup; }

/*!
 * Internal class. Forward-declared here since it is part of the writer.
 */
namespace internal { struct EclipseSummary; }

/*!
 * \brief A class to write the reservoir state and the well state of a
 *        blackoil simulation to disk using the Eclipse binary format.
 *
 * This class only writes files if the 'write_output' parameter is set
 * to 1. It needs the ERT libraries to write to disk, so if the
 * 'write_output' parameter is set but ERT is not available, all
 * methods throw a std::runtime_error.
 */
class EclipseWriter : public OutputWriter
{
public:
    /*!
     * \brief Sets the common attributes required to write eclipse
     * binary files using ERT.
     */
    EclipseWriter(const parameter::ParameterGroup& params,
                  std::shared_ptr <const EclipseGridParser> parser,
                  std::shared_ptr <const UnstructuredGrid> grid);

    /*!
     * \brief Write the static eclipse data (grid, PVT curves, etc) to disk
     */
    virtual void writeInit(const SimulatorTimer &timer);

    /*!
     * \brief Write a blackoil reservoir state to disk for later inspection with
     *        visualization tools like ResInsight
     *
     * \param[in] reservoirState The thermodynamic state of the reservoir
     * \param[in] wellState The production/injection data for all wells
     */
    virtual void writeTimeStep(const SimulatorTimer& timer,
                               const BlackoilState& reservoirState,
                               const WellState& wellState);

private:
    std::shared_ptr <const EclipseGridParser> parser_;
    std::shared_ptr <const UnstructuredGrid> grid_;
    std::string outputDir_;
    std::string baseName_;
    PhaseUsage uses_;           // active phases in the input deck

    std::unique_ptr <internal::EclipseSummary> sum_;
};
} // namespace Opm


#endif // OPM_ECLIPSE_WRITER_HPP
