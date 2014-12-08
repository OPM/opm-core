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
#include <opm/core/wells.h> // WellType

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <string>
#include <vector>
#include <array>
#include <memory>

struct UnstructuredGrid;

namespace Opm {

// forward declarations
namespace EclipseWriterDetails {
class Summary;
}

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
                  Opm::EclipseStateConstPtr eclipseState,
                  const Opm::PhaseUsage &phaseUsage,
                  int numCells,
                  const int* compressedToCartesianCellIdx);

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

    static int eclipseWellTypeMask(WellType wellType, WellInjector::TypeEnum injectorType);
    static int eclipseWellStatusMask(WellCommon::StatusEnum wellStatus);

private:
    Opm::EclipseStateConstPtr eclipseState_;
    int numCells_;
    std::array<int, 3> cartesianSize_;
    const int* compressedToCartesianCellIdx_;
    std::vector< int > gridToEclipseIdx_;
    double deckToSiPressure_;
    bool enableOutput_;
    int outputInterval_;
    int reportStepIdx_;
    std::string outputDir_;
    std::string baseName_;
    PhaseUsage phaseUsage_; // active phases in the input deck
    std::shared_ptr<EclipseWriterDetails::Summary> summary_;

    void init(const parameter::ParameterGroup& params);
};

typedef std::shared_ptr<EclipseWriter> EclipseWriterPtr;
typedef std::shared_ptr<const EclipseWriter> EclipseWriterConstPtr;

} // namespace Opm


#endif // OPM_ECLIPSE_WRITER_HPP
