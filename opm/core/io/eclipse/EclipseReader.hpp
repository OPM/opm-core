#ifndef ECLIPSEREADER_HPP
#define ECLIPSEREADER_HPP

#include <string>

#include <opm/core/simulator/WellState.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>



namespace Opm
{
///
/// \brief init_from_restart_file
///     Reading from the restart file, information stored under the OPM_XWEL keyword and SOLUTION data is in this method filled into
///     an instance of a wellstate object and a SimulatorState object.
/// \param grid
///     UnstructuredGrid reference
/// \param pu
///     PhaseUsage reference
/// \param simulator_state
///     An instance of a SimulatorState object
/// \param wellstate
///     An instance of a WellState object, with correct size for each of the 5 contained std::vector<double> objects
///

    class SimulationDataContainer;

    void init_from_restart_file(EclipseStateConstPtr eclipse_state,
                                int numcells,
                                const PhaseUsage& pu,
                                SimulationDataContainer& simulator_state,
                                WellState& wellstate);


}

#endif // ECLIPSEREADER_HPP
