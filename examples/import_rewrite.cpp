#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/io/eclipse/EclipseGridParser.hpp>

#include <boost/filesystem/convenience.hpp>

#ifdef HAVE_ERT
#include <opm/core/io/eclipse/writeECLData.hpp>
#include <ert/util/util.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_endian_flip.h>
#include <ert/ecl/fortio.h>
#endif

/*
  Small utility to read through an ECLIPSE input deck and replace
  occurences of (large) numerical fields like COORD and ZCORN with
  IMPORT statements of a binary version of the same keyword. The
  program will follow INCLUDE statements.

  Usage:  import_rewrite eclipse_case.data
*/


/*
  The numerical keywords in the ECLIPSE datafile like e.g. PORO and
  COORD are not annoted with type information; however when read and
  written in binary form they are of type float. If the updated
  datafile should be used with ECLIPSE these float values must be
  exported as float; this is achieved by setting the outputFloatType
  variable to ECL_FLOAT_TYPE.

  In OPM all numerical fields are treated as double, hence if the OPM
  EclipseParser meets an import of a float keyword it will be
  converted to double. If the output from this little utility will
  only be used from OPM the output can be saved as double directly by
  setting the outputFloatType to ECL_DOUBLE_TYPE.
*/
const ecl_type_enum outputFloatType = ECL_DOUBLE_TYPE;


/*
  Only keywords which have more >= minImportSize elements are
  converted to binary form. This is to avoid conversion of short
  keywords like MAPAXES and TABDIMS.  
*/
const int minImportSize = 10;


using namespace Opm;


static void skipKeyword( std::ifstream& is) {
  std::string keyword;
  EclipseGridParser::readKeyword( is , keyword );
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



static ecl_kw_type * loadFromcstdio( const std::string& filename , std::ios::pos_type& offset , ecl_type_enum ecl_type) {
  ecl_kw_type * ecl_kw;

  FILE * cstream = util_fopen( filename.c_str() , "r");
  fseek( cstream , offset , SEEK_SET);
  ecl_kw = ecl_kw_fscanf_alloc_current_grdecl( cstream , ecl_type );
  offset = ftell( cstream );
  fclose( cstream );

  return ecl_kw;
}



static bool convertKeyword( const std::string& inputFile , const std::string& outputPath , std::string& outputFile , std::ifstream& is , FieldType fieldType , std::ofstream& os ) {
  bool convert = true;
  ecl_type_enum ecl_type;

  if (fieldType == Integer)
    ecl_type = ECL_INT_TYPE;
  else
    ecl_type = outputFloatType;

  {
    std::ios::pos_type inputPos = is.tellg();
    ecl_kw_type * ecl_kw = loadFromcstdio( inputFile , inputPos , ecl_type );
    
    if (ecl_kw_get_size( ecl_kw ) >= minImportSize) {
      {
        if (outputPath.empty())
          outputFile = ecl_kw_get_header( ecl_kw );
        else
          outputFile = outputPath + "/" + ecl_kw_get_header( ecl_kw );
        {
          fortio_type * fortio = fortio_open_writer( outputFile.c_str() , false , ECL_ENDIAN_FLIP );
          ecl_kw_fwrite( ecl_kw , fortio );
          fortio_fclose( fortio );
        }
        
        os << "IMPORT" << std::endl << "  '" << outputFile << "'  /" << std::endl << std::endl;
      }
      is.seekg( inputPos );
    } else {
      copyKeyword( is , os );
      convert = false;
    }
    

    ecl_kw_free( ecl_kw );
  }
  return convert;
}





static bool parseFile(const std::string& inputFile, std::string& outputFile, const std::string& indent = "") {
  bool updateFile = false;
  std::cout << indent << "Parsing " << inputFile << "\n";
  {
    std::ifstream is(inputFile.c_str());  
    if (is) {
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
      
      while(is.good()) {
        is >> ignoreWhitespace;
        { 
          std::ios::pos_type start_pos = is.tellg();
          if (EclipseGridParser::readKeyword( is , keyword )) {
            FieldType fieldType = EclipseGridParser::classifyKeyword( keyword );
            switch (fieldType) {
            case(Integer):
            case(FloatingPoint):
              {
                std::string keywordFile;
                is.seekg( start_pos );
                if (convertKeyword( inputFile , path , keywordFile , is , fieldType , os )) {
                  std::cout << indent  + "   " << "Writing binary file: " << keywordFile << std::endl;
                  updateFile = true;
                }
                break;
              }
            case(Include):
              {
                std::string includeFile = readString(is);
                if (!path.empty()) {
                  includeFile = path + '/' + includeFile;
                }
                {
                  std::string __outputFile;
                  bool updateInclude = parseFile( includeFile , __outputFile , indent + "   ");
                  if (updateInclude) {
                    is.seekg( start_pos );
                    skipKeyword( is );
                    os << "INCLUDE" << std::endl << "   '" << __outputFile << "'  /" << std::endl << std::endl;
                  }
                  updateFile |= updateInclude;
                }
                break;
              }
            default:
              {
                is.seekg( start_pos );
                copyKeyword( is , os);
                break;
              }
            }
          } else
            is >> ignoreLine;  // Not at a valid keyword
        }
      }
      
      os.close();
      is.close();
      if (updateFile)
        std::cout << indent << "Written updated include file: " << outputFile << std::endl;
      else
        remove( outputFile.c_str() );
    } else
      std::cerr << indent << "** WARNING: Failed to open include file: " << inputFile << " for reading **" << std::endl;
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
    parseFile(argv[1] , outputFile);
  }
}
