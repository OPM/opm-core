#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/grid.h>

#include <cmath>

using namespace Opm;

bool
SimulatorState::equals (const SimulatorState& other,
                        double epsilon) const {
    bool equal = (num_phases_ == other.num_phases_);

    // if we use &=, then all the tests will be run regardless
    equal = equal && vectorApproxEqual( pressure() , other.pressure() , epsilon);
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
SimulatorState::init(const UnstructuredGrid& g, int num_phases)
{
    init(g.number_of_cells, g.number_of_faces, num_phases);
}
void
SimulatorState::init(int number_of_cells, int number_of_faces, int num_phases)
{
    num_phases_ = num_phases;
    press_.resize(number_of_cells, 0.0);
    fpress_.resize(number_of_faces, 0.0);
    flux_.resize(number_of_faces, 0.0);
    sat_.resize(num_phases_ * number_of_cells, 0.0);
    for (int cell = 0; cell < number_of_cells; ++cell) {
        // Defaulting the second saturation to 1.0.
        // This will usually be oil in a water-oil case,
        // gas in an oil-gas case.
        // For proper initialization, one should not rely on this,
        // but use available phase information instead.
        sat_[num_phases_*cell + 1] = 1.0;
    }
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
    const double* svals = (es == MinSat) ? &smin[0] : &smax[0];
    for (int ci = 0; ci < n; ++ci) {
        const int cell = cells[ci];
        sat_[num_phases_*cell] = svals[num_phases_*ci];
        sat_[num_phases_*cell + 1] = 1.0 - sat_[num_phases_*cell];
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
