#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/eclipse/EclipseGridParser.hpp>

#include <boost/filesystem/convenience.hpp>

#ifdef HAVE_ERT
#include <opm/core/utility/writeECLData.hpp>
#include <util.h>
#include <ecl_util.h>
#include <ecl_kw.h>
#include <ecl_endian_flip.h>
#include <fortio.h>
#endif


using namespace Opm;


static void skipKeyword( std::ifstream& is) {
  std::string keyword;
  EclipseGridParser::readKeyword( is , keyword );
  std::cout << "Skipping: " << keyword << "Pos: " << is.tellg() << std::endl;
  while (true) {
    std::ios::pos_type pos = is.tellg();

    if (EclipseGridParser::readKeyword( is , keyword )) {
      is.seekg( pos );  // Repos to start of keyword for next read.
      break;
    } else
      is >> ignoreLine;
    
    if (!is.good()) {
      is.clear();
      is.seekg( 0 , std::ios::end );
      break;
    }
  }
}


static void copyKeyword( std::ifstream& is , std::ofstream& os) {
  std::ios::pos_type start_pos = is.tellg();
  skipKeyword( is );
  {
    std::ios::pos_type end_pos = is.tellg();
    long length = end_pos - start_pos;
    
    {
      char * buffer = new char[length];
      {
        is.seekg( start_pos );
        is.read( buffer , length );
      }
      os.write( buffer , length );
      delete[] buffer;
    }
  }
}





static void convertKeyword( const std::string& inputFile , const std::string& outputPath , std::ifstream& is , FieldType fieldType , std::ofstream& os ) {
  const ecl_type_enum outputFloatType = ECL_DOUBLE_TYPE;
  ecl_type_enum ecl_type;
  ecl_kw_type * ecl_kw;

  if (fieldType == Integer)
    ecl_type = ECL_INT_TYPE;
  else
    ecl_type = outputFloatType;

  std::cout << "InputPos: " << is.tellg() << std::endl;
  {
    FILE * cstream = util_fopen( inputFile.c_str() , "r");
    fseek( cstream , is.tellg() , SEEK_SET);
    ecl_kw = ecl_kw_fscanf_alloc_current_grdecl( cstream , ecl_type );
    {
      std::ios::pos_type pos = ftell( cstream );
      is.seekg( pos , std::ios::beg );
    }
    util_fclose( cstream );
    
    {
      std::string outputFile = outputPath + "/" + ecl_kw_get_header( ecl_kw );
      fortio_type * fortio = fortio_open_writer( outputFile.c_str() , false , ECL_ENDIAN_FLIP );
      ecl_kw_fwrite( ecl_kw , fortio );
      std::cout << "Writing binary file: " << ecl_kw_get_header( ecl_kw ) << std::endl;
      fortio_fclose( fortio );

      os << "IMPORT" << std::endl << "  '" << outputFile << "'  /" << std::endl << std::endl;
    }
  }
  std::cout << "Exit: InputPos: " << is.tellg() << std::endl;
}


bool parseFile(const std::string& inputFile, std::string& outputFile) {
  bool updateFile = false;
  std::cout << "Reading file " << inputFile << "\n";
  {
    std::ifstream is(inputFile.c_str());
    std::ofstream os;
    std::string keyword;
    std::string path;
    {
      boost::filesystem::path inputPath(inputFile);
      path = inputPath.parent_path().string();
    }

    {
      std::string basename;
      std::string extension;
      
      outputFile = inputFile;
      size_t ext_pos = inputFile.rfind(".");
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
          std::cout << std::endl;
          std::cout << "Have read keyword: " << keyword << " Type : " << fieldType << "\n";
          
          switch (fieldType) {
          case(Integer):
            {
              is.seekg( start_pos );
              std::cout << "Converting keyword: " << keyword << std::endl;
              convertKeyword( inputFile , path , is , fieldType , os );
              updateFile = true;
              break;
            }
          case(FloatingPoint):
            {
              is.seekg( start_pos );
              std::cout << "Converting keyword: " << keyword << std::endl;
              convertKeyword( inputFile , path , is , fieldType , os );
              updateFile = true;
              break;
            }
          case(Include):
            {
              std::string includeFile = readString(is);
              if (!path.empty()) {
                includeFile = path + '/' + includeFile;
              }
              std::cout << "Continue to includeFile: " << includeFile << std::endl;
              {
                bool updateInclude = parseFile( includeFile , outputFile );
                if (updateInclude) {
                  is.seekg( start_pos );
                  skipKeyword( is );
                  os << "INCLUDE" << std::endl << "   '" << outputFile << "'  /" << std::endl << std::endl;
                }
                updateFile |= updateInclude;
              }
              break;
            }
          default:
            {
              std::cout << "Calling copy: " << keyword << " fieldType:" << fieldType << std::endl;
              is.seekg( start_pos );
              copyKeyword( is , os);
              std::cout << "After loop pos: " << is.tellg() << std::endl;
              break;
            }
          }
        } else
          is >> ignoreLine;  // Not at a valid keyword
      }
    }
    
    os.close();
    is.close();
    if (!updateFile)
      remove( outputFile.c_str() );
  }
  return updateFile;
}



int
main(int argc, char** argv)
{
  if (argc != 2)
    THROW("Need the name of ECLIPSE file on command line");
  {
    std::string outputFile;
    if (parseFile(argv[1] , outputFile))
      std::cout << "Have written: " << outputFile << std::endl;
  }
}
