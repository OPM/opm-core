#ifndef OPM_SIMULATORINCOMPTWOPHASE_HEADER_INCLUDED
#error Do not include SimulatorIncompTwophase directly!
#endif

namespace Opm {

template <typename T, void (T::*callback)()>
inline void SimulatorIncompTwophase::connect_timestep (T& t) {
	connect_timestep_impl (boost::function0<void> (std::bind (callback, t)));
}

} /* namespace Opm */
