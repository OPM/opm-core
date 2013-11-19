/*
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

#ifndef OPM_OUTPUT_WRITER_HPP
#define OPM_OUTPUT_WRITER_HPP

#include <memory>  // unique_ptr, shared_ptr

struct UnstructuredGrid;

namespace Opm {

// forward declaration
class BlackoilState;
class EclipseGridParser;
namespace parameter { class ParameterGroup; }
class SimulatorTimer;
class WellState;

/*!
 * Interface for writing non-compositional (blackoil, two-phase) simulation
 * state to files.
 *
 * Use the create() function to setup a chain of writer based on the
 * configuration values, e.g.
 *
 * \example
 * \code{.cpp}
 *  ParameterGroup params (argc, argv, false);
 *  auto parser = std::make_shared <EclipseGridParser> (
 *                      params.get <string> ("deck_filename"));
 *
 *  std::unique_ptr <OutputWriter> writer =
 *          OutputWriter::create (params, parser);
 *
 *  // before the first timestep
 *  writer->writeInit (timer);
 *
 *  // after each timestep
 *  writer->writeTimeStep (timer, state, wellState);
 *
 * \endcode
 */
class OutputWriter {
public:
    /*!
     * \brief Write the static eclipse data (grid, PVT curves, etc) to disk
     */
    virtual void writeInit(const SimulatorTimer &timer) = 0;

    /*!
     * \brief Write a blackoil reservoir state to disk for later inspection with
     *        visualization tools like ResInsight
     *
     * \param[in] reservoirState The thermodynamic state of the reservoir
     * \param[in] wellState The production/injection data for all wells
     */
    virtual void writeTimeStep(const SimulatorTimer& timer,
                                 const BlackoilState& reservoirState,
                                 const WellState& wellState) = 0;

    /*!
     * Create a suitable set of output formats based on configuration.
     *
     * @param params Configuration properties. This function will setup a
     *               multiplexer of applicable output formats based on the
     *               desired configuration values.
     *
     * @param parser Input deck used to set up the simulation. The lifetime
     *               of this object must exceed the lifetime of the writer
     *               that is returned.
     *
     * @return       Pointer to a multiplexer to all applicable output formats.
     *
     * @see Opm::share_obj
     */
    static std::unique_ptr <OutputWriter>
    create (const parameter::ParameterGroup& params,
            std::shared_ptr <EclipseGridParser> parser,
            std::shared_ptr <UnstructuredGrid> grid);
};

} // namespace Opm

#endif /* OPM_OUTPUT_WRITER_HPP */
