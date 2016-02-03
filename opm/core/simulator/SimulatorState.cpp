#include <opm/common/ErrorMacros.hpp>
#include <opm/core/simulator/SimulatorState.hpp>

#include <cmath>
#include <cassert>

using namespace Opm;

bool
SimulatorState::equals (const SimulatorState& other,
                        double epsilon) const {
    bool equal = (num_phases_ == other.num_phases_);

    // if we use &=, then all the tests will be run regardless
    equal = equal && vectorApproxEqual( pressure() , other.pressure() , epsilon);
    equal = equal && vectorApproxEqual( temperature() , other.temperature() , epsilon);
    equal = equal && vectorApproxEqual( facepressure() , other.facepressure() , epsilon);
    equal = equal && vectorApproxEqual( faceflux() , other.faceflux() , epsilon);
    equal = equal && vectorApproxEqual( saturation() , other.saturation() , epsilon);

    return equal;
}

bool
SimulatorState::vectorApproxEqual(const std::vector<double>& v1,
                                  const std::vector<double>& v2,
                                  double epsilon) {
    if (v1.size() != v2.size()) {
        return false;
    }

    for (size_t i = 0; i < v1.size(); i++) {
        const double diff = std::abs(v1[i] - v2[i]);
        const double scale = std::abs(v1[i]) + std::abs(v2[i]);
        if (diff > epsilon * scale) {
            return false;
        }
    }

    return true;
}


void
SimulatorState::init(int number_of_cells, int number_of_faces, int num_phases)
{
    num_cells_  = number_of_cells;
    num_faces_  = number_of_faces;
    num_phases_ = num_phases;

    // clear memory
    cellData_ = std::vector< std::vector<double> > ();
    faceData_ = std::vector< std::vector<double> > ();

    int id;
    id = registerCellData("PRESSURE", 1, 0.0 );
    if( id != pressureId_ )
        OPM_THROW(std::logic_error,"ids in SimulatorState do not match");
    assert( pressureId_ == id );
    id = registerCellData("TEMPERATURE", 1, 273.15 + 20 );
    assert( temperatureId_ == id );
    id = registerCellData("SATURATION", num_phases_, 0.0 );
    assert( saturationId_ == id );

    for (int cell = 0; cell < number_of_cells; ++cell) {
        // Defaulting the second saturation to 1.0.
        // This will usually be oil in a water-oil case,
        // gas in an oil-gas case.
        // For proper initialization, one should not rely on this,
        // but use available phase information instead.
        saturation()[num_phases_*cell + 1] = 1.0;
    }

    id = registerFaceData("FACEPRESSURE", 1, 0.0 );
    assert( facePressureId_ == id );
    id = registerFaceData("FACEFLUX", 1, 0.0 );
    assert( faceFluxId_ == id );
}

size_t
SimulatorState::registerCellData( const std::string& name, const int components, const double initialValue )
{
    // check if init has been called
    const size_t pos = cellData_.size();
    cellDataNames_.emplace_back( name );
    cellData_.emplace_back( num_cells_ * components, initialValue );
    return pos;
}

size_t
SimulatorState::registerFaceData( const std::string& name, const int components, const double initialValue )
{
    // check if init has been called
    const size_t pos = faceData_.size();
    faceDataNames_.emplace_back( name );
    faceData_.emplace_back( num_faces_ * components, initialValue );
    return pos ;
}

