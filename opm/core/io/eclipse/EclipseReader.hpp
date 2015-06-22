#ifndef ECLIPSEREADER_HPP
#define ECLIPSEREADER_HPP

#include <string>
#include <opm/core/simulator/WellState.hpp>

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
}

#endif // ECLIPSEREADER_HPP
