/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

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



#include <opm/core/grid.h>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/writeECLData.hpp>
#include <opm/core/utility/Units.hpp>

#include <vector>

#include <ecl_grid.h>
#include <ecl_util.h>
#include <ecl_rst_file.h>



namespace Opm
{
  
  static ecl_kw_type * ecl_kw_wrapper( const UnstructuredGrid& grid,
                                       const std::string& kw_name , 
                                       const std::vector<double> * data , 
                                       int offset , 
                                       int stride ) {

    ecl_kw_type * ecl_kw = ecl_kw_alloc( kw_name.c_str() , data->size() / stride , ECL_FLOAT_TYPE );
    if (grid.global_cell == NULL) {
      for (int i=0; i < grid.number_of_cells; i++) 
        ecl_kw_iset_float( ecl_kw , i , (*data)[i*stride + offset]);
    } else {
      for (int i=0; i < grid.number_of_cells; i++) 
        ecl_kw_iset_float( ecl_kw , grid.global_cell[i] , (*data)[i*stride + offset]);
    }
    return ecl_kw;
  }



  void writeECLData(const UnstructuredGrid& grid,
                    const DataMap& data,
                    const SimulatorTimer& simtimer,
                    const std::string& output_dir,
                    const std::string& base_name) {
    
    int step                = simtimer.currentStepNum();
    ecl_file_enum file_type = ECL_UNIFIED_RESTART_FILE;
    bool fmt_file           = true; 
    char * filename         = ecl_util_alloc_filename(output_dir.c_str() , base_name.c_str() , file_type , fmt_file , step );
    int phases              = ECL_OIL_PHASE + ECL_WATER_PHASE;
    double days             = Opm::unit::convert::to(simtimer.currentTime(), Opm::unit::day);
    time_t date             = 0;
    int nx                  = grid.cartdims[0];
    int ny                  = grid.cartdims[1];
    int nz                  = grid.cartdims[2];
    int nactive             = grid.number_of_cells;
    ecl_rst_file_type * rst_file;
    
    if (step > 0 && file_type == ECL_UNIFIED_RESTART_FILE)
      rst_file = ecl_rst_file_open_append( filename );
    else
      rst_file = ecl_rst_file_open_write( filename );
    
    ecl_rst_file_fwrite_header( rst_file , step , date , days , nx , ny , nz , nactive , phases );
    ecl_rst_file_start_solution( rst_file );

    {
      DataMap::const_iterator i = data.find("pressure");
      if (i != data.end()) {
        ecl_kw_type * pressure_kw = ecl_kw_wrapper( grid , "PRESSURE" , i->second , 0 , 1);
        ecl_rst_file_add_kw( rst_file , pressure_kw );
        ecl_kw_free( pressure_kw );
      }
    }
    
    {
      DataMap::const_iterator i = data.find("saturation");
      if (i != data.end()) {
        ecl_kw_type * swat_kw = ecl_kw_wrapper( grid , "SWAT" , i->second , 0 , 2);
        ecl_rst_file_add_kw( rst_file , swat_kw );
        ecl_kw_free( swat_kw );
      }
    }

    ecl_rst_file_end_solution( rst_file );
    ecl_rst_file_close( rst_file );
    free(filename);
  } 
}

