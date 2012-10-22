#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/initState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/SimpleFluid2p.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

#include <opm/core/utility/ColumnExtract.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/transport/GravityColumnSolver.hpp>

#include <opm/core/transport/reorder/TransportModelTwophase.hpp>

#include <boost/filesystem/convenience.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <cassert>
#include <cstddef>
//#include <cstdio.h>

#include <algorithm>
#include <tr1/array>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <vector>
#include <numeric>

#ifdef HAVE_ERT
#include <opm/core/utility/writeECLData.hpp>

#include <util.h>

#include <ecl_util.h>
#include <ecl_kw.h>
#include <ecl_endian_flip.h>
#include <fortio.h>

#endif


using namespace Opm;



static void skipKeyword( std::ifstream& is , std::ofstream& os) {
  std::ios::pos_type start_pos = is.tellg();
  std::ios::pos_type end_pos;
  size_t length;
  {
    std::string keyword;
    std::cout << "Starting in " << start_pos;
    EclipseGridParser::readKeyword( is , keyword );
    while (true) {
      
      end_pos = is.tellg();
      if (EclipseGridParser::readKeyword( is , keyword )) 
        break;
      else
        is >> ignoreLine;
      
      if (!is.good()) {
        is.clear();
        is.seekg( 0 , std::ios::end );
        end_pos = is.tellg();
        break;
      }
    }
    
    length = end_pos - start_pos;
  }

  {
    char * buffer = new char[length];
    {
      is.seekg( start_pos );
      is.read( buffer , length );
    }
    os.write( buffer , length );
    delete[] buffer;
  }
  std::cout << "Stopping in " << end_pos << std::endl;
}


static void convertKeyword( const char * inputFile , std::ifstream& is , FieldType fieldType , std::ofstream& os ) {
  const ecl_type_enum outputFloatType = ECL_DOUBLE_TYPE;
  ecl_type_enum ecl_type;
  ecl_kw_type * ecl_kw;

  if (fieldType == Integer)
    ecl_type = ECL_INT_TYPE;
  else
    ecl_type = outputFloatType;

  {
    FILE * cstream = util_fopen( inputFile , "r");
    fseek( cstream , is.tellg() , SEEK_SET);
    ecl_kw = ecl_kw_fscanf_alloc_current_grdecl( cstream , ecl_type );
    {
      std::ios::pos_type pos = ftell( cstream );
      is.seekg( pos , std::ios::beg );
    }
    util_fclose( cstream );
    
    {
      fortio_type * fortio = fortio_open_writer( ecl_kw_get_header( ecl_kw ) , false , ECL_ENDIAN_FLIP );
      ecl_kw_fwrite( ecl_kw , fortio );
      std::cout << "Writing binary file: " << ecl_kw_get_header( ecl_kw ) << std::endl;
      fortio_fclose( fortio );
    }
  }
}




int
main(int argc, char** argv)
{
  if (argc != 2)
    THROW("Need the name of ECLIPSE file on command line");

  std::cout << "Reading file " << argv[1] << "\n";
  {
    std::ifstream is(argv[1]);
    std::ofstream os;
    std::string outputFile(argv[1]);
    std::string keyword;


    {
      std::string basename;
      std::string extension;
      
      size_t ext_pos = outputFile.rfind(".");
      if (ext_pos == std::string::npos) {
        basename = outputFile.substr();
        extension = "";
      } else {
        basename = outputFile.substr(0,ext_pos);
        extension = outputFile.substr(ext_pos);
      }
      
      outputFile = basename + "_import" + extension;
    }
    os.open( outputFile.c_str() );
    std::cout << "Writing to file: " << outputFile << "\n";
    
    while(is.good()) {
      is >> ignoreWhitespace;
      { 
        std::ios::pos_type start_pos = is.tellg();
        if (EclipseGridParser::readKeyword( is , keyword )) {
          FieldType fieldType = EclipseGridParser::classifyKeyword( keyword );
          std::cout << "Have read keyword: " << keyword << " Type : " << fieldType << "\n";
          
          is.seekg( start_pos );
          if (fieldType == Integer || fieldType == FloatingPoint) 
            convertKeyword( argv[1] , is , fieldType , os );
          else
            skipKeyword( is , os );
        } else
          is >> ignoreLine;
      }
      
    }
    
    os.close();
    is.close();
  }
}
