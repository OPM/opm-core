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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/core/io/eclipse/writeECLData.hpp>

#include <opm/core/grid.h>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <vector>

#ifdef HAVE_ERT // This one goes almost to the bottom of the file

#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_rst_file.h>

#include <cstdlib>
#include <ctime>
#include <memory>

namespace {
  std::shared_ptr<ecl_kw_type>
  ecl_kw_wrapper( const UnstructuredGrid& grid,
                  const std::string& kw_name ,
                  const std::vector<double> * data ,
                  int offset ,
                  int stride ) {

    if (stride <= 0)
      OPM_THROW(std::runtime_error, "Vector strides must be positive. Got stride = " << stride);
    if ((stride * std::vector<double>::size_type(grid.number_of_cells)) != data->size())
      OPM_THROW(std::runtime_error, "Internal mismatch grid.number_of_cells: " << grid.number_of_cells << " data size: " << data->size() / stride);
    {
      std::shared_ptr<ecl_kw_type>
          ecl_kw(ecl_kw_alloc( kw_name.c_str() , grid.number_of_cells , ECL_FLOAT_TYPE ),
                 ecl_kw_free);
      for (int i=0; i < grid.number_of_cells; i++)
        ecl_kw_iset_float( ecl_kw.get() , i , (*data)[i*stride + offset]);
      return ecl_kw;
    }
  }

  template <ecl_file_enum FileType  = ECL_UNIFIED_RESTART_FILE,
            bool          Formatted = false>
  class RestartFile {
  public:
    RestartFile(const UnstructuredGrid&         grid,
                const int                       current_step,
                const double                    current_time,
                const boost::posix_time::ptime& current_date_time,
                const std::string&              output_dir,
                const std::string&              base_name)
      : nactive_(grid.number_of_cells)
      , rst_file_(0)
    {
      const int phases  = ECL_OIL_PHASE + ECL_WATER_PHASE;
      const double days = Opm::unit::convert::to(current_time, Opm::unit::day);
      const int nx      = grid.cartdims[0];
      const int ny      = grid.cartdims[1];
      const int nz      = grid.cartdims[2];

      std::time_t date  = 0;
      {
        using namespace boost::posix_time;
        ptime t0( boost::gregorian::date(1970 , 1 ,1) );
        time_duration::sec_type seconds = (current_date_time - t0).total_seconds();

        date = std::time_t( seconds );
      }

      {
        std::shared_ptr<char>
          filename(ecl_util_alloc_filename(output_dir.c_str() ,
                                           base_name.c_str()  ,
                                           FileType, Formatted,
                                           current_step       ),
                   std::free);

        if (current_step > 0 && FileType == ECL_UNIFIED_RESTART_FILE)
          rst_file_ = ecl_rst_file_open_append( filename.get() );
        else
          rst_file_ = ecl_rst_file_open_write( filename.get() );
      }

      ecl_rst_file_fwrite_header( rst_file_ , current_step , date , days ,
                                  nx , ny , nz , nactive_ , phases );
      ecl_rst_file_start_solution( rst_file_ );
    }

    ~RestartFile()
    {
      ecl_rst_file_end_solution( rst_file_ );
      ecl_rst_file_close( rst_file_ );
    }

  private:
    const int          nactive_;
    ecl_rst_file_type* rst_file_;
  };
} // Anonymous namespace


namespace Opm
{
  /*
    This function will write the data solution data in the DataMap
    @data as an ECLIPSE restart file, in addition to the solution
    fields the ECLIPSE restart file will have a minimum (hopefully
    sufficient) amount of header information.

    The ECLIPSE restart files come in two varietes; unified restart
    files which have all the report steps lumped together in one large
    chunk and non-unified restart files which are one file for each
    report step. In addition the files can be either formatted
    (i.e. ASCII) or unformatted (i.e. binary).

    The writeECLData() function has two hardcoded settings:
    'file_type' and 'fmt_file' which regulate which type of files the
    should be created. The extension of the files follow a convention:

      Unified, formatted    : .FUNRST
      Unified, unformatted  : .UNRST
      Multiple, formatted   : .Fnnnn
      Multiple, unformatted : .Xnnnn

    For the multiple files the 'nnnn' part is the report number,
    formatted with '%04d' format specifier. The writeECLData()
    function will use the ecl_util_alloc_filename() function to create
    an ECLIPSE filename according to this conventions.
  */

  void writeECLData(const UnstructuredGrid& grid,
                    const DataMap& data,
                    const int current_step,
                    const double current_time,
                    const boost::posix_time::ptime& current_date_time,
                    const std::string& output_dir,
                    const std::string& base_name) {

    ecl_file_enum file_type = ECL_UNIFIED_RESTART_FILE;  // Alternatively ECL_RESTART_FILE for multiple restart files.
    bool fmt_file           = false;

    std::shared_ptr<char> filename(ecl_util_alloc_filename(output_dir.c_str() ,
                                                           base_name.c_str()  ,
                                                           file_type, fmt_file,
                                                           current_step       ),
                                   std::free);
    int phases              = ECL_OIL_PHASE + ECL_WATER_PHASE;
    double days             = Opm::unit::convert::to(current_time, Opm::unit::day);
    time_t date             = 0;
    int nx                  = grid.cartdims[0];
    int ny                  = grid.cartdims[1];
    int nz                  = grid.cartdims[2];
    int nactive             = grid.number_of_cells;
    ecl_rst_file_type * rst_file;

    {
      using namespace boost::posix_time;
      ptime t0( boost::gregorian::date(1970 , 1 ,1) );
      time_duration::sec_type seconds = (current_date_time - t0).total_seconds();

      date = time_t( seconds );
    }

    if (current_step > 0 && file_type == ECL_UNIFIED_RESTART_FILE)
      rst_file = ecl_rst_file_open_append( filename.get() );
    else
      rst_file = ecl_rst_file_open_write( filename.get() );

    ecl_rst_file_fwrite_header( rst_file , current_step , date , days , nx , ny , nz , nactive , phases );
    ecl_rst_file_start_solution( rst_file );

    {
      DataMap::const_iterator i = data.find("pressure");
      if (i != data.end()) {
        std::shared_ptr<ecl_kw_type>
            pressure_kw = ecl_kw_wrapper( grid , "PRESSURE" , i->second , 0 , 1);
        ecl_rst_file_add_kw( rst_file , pressure_kw.get() );
      }
    }

    {
      DataMap::const_iterator i = data.find("saturation");
      if (i != data.end()) {
        if (int(i->second->size()) != 2 * grid.number_of_cells) {
          OPM_THROW(std::runtime_error, "writeECLData() requires saturation field to have two phases.");
        }
        std::shared_ptr<ecl_kw_type>
            swat_kw = ecl_kw_wrapper( grid , "SWAT" , i->second , 0 , 2);
        ecl_rst_file_add_kw( rst_file , swat_kw.get() );
      }
    }

    ecl_rst_file_end_solution( rst_file );
    ecl_rst_file_close( rst_file );
  }
}

#else // that is, we have not defined HAVE_ERT

namespace Opm
{

    void writeECLData(const UnstructuredGrid&,
                      const DataMap&,
                      const int,
                      const double,
                      const boost::posix_time::ptime&,
                      const std::string&,
                      const std::string&)
    {
        OPM_THROW(std::runtime_error, "Cannot call writeECLData() without ERT library support. Reconfigure opm-core with ERT support and recompile.");
    }
}

#endif
