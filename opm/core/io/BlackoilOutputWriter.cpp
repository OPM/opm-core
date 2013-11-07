#include "BlackoilOutputWriter.hpp"

#include <opm/core/io/eclipse/BlackoilEclipseOutputWriter.hpp>
#include <opm/core/utility/parameters/Parameter.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <map>
#include <memory> // unique_ptr
#include <vector>

using namespace std;
using namespace Opm;
using namespace Opm::parameter;

namespace {

/// Multiplexer over a list of output writers
struct BlackoilMultiOutputWriter : public BlackoilOutputWriter {
    /// Shorthand for a list of owned output writers
    typedef vector <unique_ptr <BlackoilOutputWriter> > writers_t;
    typedef writers_t::iterator it_t;
    typedef unique_ptr <writers_t> ptr_t;

    /// Adopt a list of writers
    BlackoilMultiOutputWriter (ptr_t writers)
        : writers_ (std::move (writers)) { }

    /// Forward the call to all writers
    virtual void writeInit(const SimulatorTimer &timer) {
        for (it_t it = writers_->begin (); it != writers_->end (); ++it) {
            (*it)->writeInit (timer);
        }
    }

    virtual void writeTimeStep(const SimulatorTimer& timer,
                                 const BlackoilState& reservoirState,
                                 const WellState& wellState) {
        for (it_t it = writers_->begin (); it != writers_->end(); ++it) {
            (*it)->writeTimeStep (timer, reservoirState, wellState);
        }
    }

private:
    ptr_t writers_;
};

/// Psuedo-constructor, can appear in template
template <typename Format> unique_ptr <BlackoilOutputWriter>
create (const ParameterGroup& params, const EclipseGridParser& parser) {
    return unique_ptr <BlackoilOutputWriter> (new Format (params, parser));
}

/// Map between keyword in configuration and the corresponding
/// constructor function (type) that should be called when detected.
/// The writer must have a constructor which takes params and parser.
///
/// If you want to add more possible writer formats, just add them
/// to the list below!
typedef map <const char*,
              unique_ptr <BlackoilOutputWriter> (*)(const ParameterGroup&,
                                                    const EclipseGridParser&)> map_t;
map_t FORMATS = {
    { "output_ecl", &create <BlackoilEclipseOutputWriter> },
};

} // anonymous namespace

unique_ptr <BlackoilOutputWriter>
BlackoilOutputWriter::create (const ParameterGroup& params,
                              const EclipseGridParser& parser) {
    // allocate a list which will be filled with writers. this list
    // is initially empty (no output).
    BlackoilMultiOutputWriter::ptr_t list (
                new BlackoilMultiOutputWriter::writers_t ());

    // loop through the map and see if we can find the key that is
    // specified there
    typedef map_t::iterator map_it_t;
    for (map_it_t it = FORMATS.begin (); it != FORMATS.end(); ++it) {
        // keyword which would indicate that this format should be used
        const std::string name (it->first);

        // invoke the constructor for the type if we found the keyword
        // and put the pointer to this writer onto the list
        if (params.getDefault <bool> (name, false)) {
            list->push_back (it->second (params, parser));
        }
    }

    // create a multiplexer from the list of formats we found
    return unique_ptr <BlackoilOutputWriter> (
                new BlackoilMultiOutputWriter (std::move (list)));
}
