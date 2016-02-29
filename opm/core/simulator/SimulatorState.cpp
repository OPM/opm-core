#include <opm/common/ErrorMacros.hpp>
#include <opm/common/util/numeric/cmp.hpp>
#include <opm/core/simulator/SimulatorState.hpp>

#include <algorithm>
#include <cmath>
#include <cassert>

using namespace Opm;

bool
SimulatorState::equals (const SimulatorState& other,
                        double epsilon) const {
    bool equal = (num_phases_ == other.num_phases_);

    // if we use &=, then all the tests will be run regardless
    equal = equal && cmp::vector_equal( pressure() , other.pressure() , cmp::default_abs_epsilon , epsilon);
    equal = equal && cmp::vector_equal( temperature() , other.temperature() , cmp::default_abs_epsilon , epsilon);
    equal = equal && cmp::vector_equal( facepressure() , other.facepressure() , cmp::default_abs_epsilon , epsilon);
    equal = equal && cmp::vector_equal( faceflux() , other.faceflux() , cmp::default_abs_epsilon , epsilon);
    equal = equal && cmp::vector_equal( saturation() , other.saturation() , cmp::default_abs_epsilon , epsilon);

    return equal;
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

void SimulatorState::setCellDataComponent( const std::string& name , size_t component , const std::vector<int>& cells , const std::vector<double>& values) {
  const auto iter = std::find( cellDataNames_.begin() , cellDataNames_.end() , name);
  int id = iter - cellDataNames_.begin();
  auto& data = cellData_[id];
  if (component >= size_t(num_phases_))
    throw std::invalid_argument("Invalid component");

  if (cells.size() != values.size())
    throw std::invalid_argument("size mismatch between cells and values");

  /* This is currently quite broken; the setCellDataComponent
     method assumes that the number of components in the field
     we are currently focusing on has num_phases components in
     total. This restriction should be lifted by allowing a per
     field number of components.
  */
  if (data.size() != size_t(num_phases_ * num_cells_))
    throw std::invalid_argument("Can currently only be used on fields with num_components == num_phases (i.e. saturation...) ");

  for (size_t i = 0; i < cells.size(); i++) {
    if (cells[i] < num_cells_) {
      auto field_index = cells[i] * num_phases_ + component;
      auto value = values[i];

      data[field_index] = value;
    } else {
      throw std::invalid_argument("Invalid cell number");
    }
  }
}


std::vector<double>& SimulatorState::getCellData( const std::string& name )  {
    const auto iter = std::find( cellDataNames_.begin() , cellDataNames_.end() , name);
    int id = iter - cellDataNames_.begin();
    auto& data = cellData_[id];
    return data;
}


const std::vector<double>& SimulatorState::getCellData( const std::string& name )  const {
    const auto iter = std::find( cellDataNames_.begin() , cellDataNames_.end() , name);
    int id = iter - cellDataNames_.begin();
    const auto& data = cellData_[id];
    return data;
}

