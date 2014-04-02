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

#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <string>
#include <vector>
#include <memory>  // std::unique_ptr

struct UnstructuredGrid;
struct EclipseSummary;

namespace Opm {

// forward declarations
class SimulatorState;
class SimulatorTimer;
class WellState;

namespace parameter { class ParameterGroup; }

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
     *        binary files using ERT.
     */
    EclipseWriter(const parameter::ParameterGroup& params,
                  Opm::DeckConstPtr newParserDeck,
                  std::shared_ptr <const UnstructuredGrid> grid);

    /**
     * We need a destructor in the compilation unit to avoid the
     * EclipseSummary being a complete type here.
     */
    virtual ~EclipseWriter ();

    /**
     * Write the static eclipse data (grid, PVT curves, etc) to disk.
     */
    virtual void writeInit(const SimulatorTimer &timer);

    /*!
     * \brief Write a reservoir state and summary information to disk.
     *
     *
     * The reservoir state can be inspected with visualization tools like
     * ResInsight.
     *
     * The summary information can then be visualized using tools from
     * ERT or ECLIPSE. Note that calling this method is only
     * meaningful after the first time step has been completed.
     *
     * \param[in] reservoirState The thermodynamic state of the reservoir
     * \param[in] wellState The production/injection data for all wells
     */
    virtual void writeTimeStep(const SimulatorTimer& timer,
                               const SimulatorState& reservoirState,
                               const WellState& wellState);

private:
    Opm::DeckConstPtr newParserDeck_;
    std::shared_ptr <const UnstructuredGrid> grid_;
    bool enableOutput_;
    int outputInterval_;
    int outputTimeStepIdx_;
    std::string outputDir_;
    std::string baseName_;
    PhaseUsage uses_;           // active phases in the input deck
    std::shared_ptr <EclipseSummary> summary_;

    void activeToGlobalCellData_(std::vector<double> &globalCellsBuf,
                                 const std::vector<double> &activeCellsBuf,
                                 const std::vector<double> &inactiveCellsBuf) const;

    /// Write solution field variables (pressure and saturation)
    void writeSolution_(const SimulatorTimer& timer,
                        const SimulatorState& reservoirState);
};
} // namespace Opm


#endif // OPM_ECLIPSE_WRITER_HPP
