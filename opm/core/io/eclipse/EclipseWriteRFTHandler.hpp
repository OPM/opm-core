/*
  Copyright 2015 Statoil ASA.

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

#ifndef OPM_ECLIPSE_WRITE_RFT_HANDLER_HPP
#define OPM_ECLIPSE_WRITE_RFT_HANDLER_HPP

#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/SimulatorState.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <ert/ecl/ecl_rft_node.h>
#include <ert/ecl/ecl_util.h>


namespace Opm {
    class EclipseGrid;
    class Well;

namespace EclipseWriterDetails {


    class EclipseWriteRFTHandler {

    public:
    EclipseWriteRFTHandler(const int * compressedToCartesianCellIdx, size_t numCells, size_t cartesianSize);


    void writeTimeStep(const std::string& filename,
                       const ert_ecl_unit_enum ecl_unit,
                       const SimulatorTimerInterface& simulatorTimer,
                       std::vector<std::shared_ptr< const Well >>& wells,
                       std::shared_ptr< const EclipseGrid > eclipseGrid,
                       std::vector<double>& pressure,
                       std::vector<double>& swat,
                       std::vector<double>& sgas);



    private:

    ecl_rft_node_type * createEclRFTNode(std::shared_ptr< const Well > well,
                                         const SimulatorTimerInterface& simulatorTimer,
                                         std::shared_ptr< const EclipseGrid > eclipseGrid,
                                         const std::vector<double>& pressure,
                                         const std::vector<double>& swat,
                                         const std::vector<double>& sgas);

    void initGlobalToActiveIndex(const int * compressedToCartesianCellIdx, size_t numCells, size_t cartesianSize);

    std::vector<int> globalToActiveIndex_;

    };




}//namespace EclipseWriterDetails
}//namespace Opm


#endif //OPM_ECLIPSE_WRITE_RFT_HANDLER_HPP
