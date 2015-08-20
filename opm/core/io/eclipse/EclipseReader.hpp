#ifndef ECLIPSEREADER_HPP
#define ECLIPSEREADER_HPP

#include <string>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

namespace Opm
{
///
/// \brief restoreOPM_XWELKeyword
///     Reading from the restart file, information stored under the OPM_XWEL keyword is in this method filled into
///     an instance of a wellstate object.
/// \param restart_filename
///     The filename of the restart file.
/// \param reportstep
///     The report step to restart from.
/// \param wellstate
///     An instance of a WellState object, with correct size for each of the 5 contained std::vector<double> objects.
///
    void restoreOPM_XWELKeyword(const std::string& restart_filename, int report_step, WellState& wellState);
    void restoreSOLUTIONData(const std::string& restart_filename,
                             int report_step,
                             const EclipseState &eclipseState,
                             const UnstructuredGrid& grid,
                             const PhaseUsage& phaseUsage,
                             SimulatorState& simulator_state);


}

#endif // ECLIPSEREADER_HPP
