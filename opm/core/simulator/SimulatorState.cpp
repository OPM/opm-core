#include <opm/common/ErrorMacros.hpp>
#include <opm/core/simulator/SimulatorState.hpp>

#include <cmath>
#include <cassert>

using namespace Opm;


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

template <typename Props> void
SimulatorState::setFirstSat(const std::vector<int>& cells,
                            const Props& props,
                            ExtremalSat es)
{
    if (cells.empty()) {
        return;
    }
    int n = cells.size();
    std::vector<double> smin(num_phases_*n);
    std::vector<double> smax(num_phases_*n);
    props.satRange(n, &cells[0], &smin[0], &smax[0]);
    std::vector< double >& sat = saturation();
    const double* svals = (es == MinSat) ? &smin[0] : &smax[0];
    for (int ci = 0; ci < n; ++ci) {
        const int cell = cells[ci];
        sat[num_phases_*cell] = svals[num_phases_*ci];
        sat[num_phases_*cell + 1] = 1.0 - sat[num_phases_*cell];
    }
}

// template instantiations for all known (to this library) subclasses
// of SimulatorState that will call this method. notice that there are
// no empty angle brackets after "template" -- that would have been
// specialization instead
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>

template void
SimulatorState::setFirstSat <IncompPropertiesInterface> (
        const std::vector<int> &cells,
        const IncompPropertiesInterface &props,
        ExtremalSat es);

template void
SimulatorState::setFirstSat <BlackoilPropertiesInterface> (
        const std::vector<int> &cells,
        const BlackoilPropertiesInterface &props,
        ExtremalSat es);
