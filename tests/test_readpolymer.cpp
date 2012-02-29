//===========================================================================
//                                                                           
// File: test_readpolymer.cpp                                                       
//                                                                           
// Created: Thu Jan 12 15:18:46 2012                                         
//                                                                           
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
//===========================================================================

#include "../config.h"

#include <opm/core/eclipse/EclipseGridParser.hpp>

// Test program for reading Eclipse Polymer keywords.

int main(int argc, char** argv)
{
    using namespace std;

    std::string ecl_filename;
    if (argc == 2) {
	ecl_filename = argv[1];
    } else {
	std::cout << "\nUsage: test_readpolymer  filename.grdecl\n";
	exit( 1 );
    }


    bool convert_to_SI = true;
    Opm::EclipseGridParser parser(ecl_filename, convert_to_SI);

    std::cout << "\n  Polymer fields\n\n";

    if (parser.hasField("PLYVISC")) {
	parser.getPLYVISC().write(std::cout);
    }
    if (parser.hasField("PLYROCK")) {
	parser.getPLYROCK().write(std::cout);
    }
    if (parser.hasField("PLYADS")) { 
	parser.getPLYADS().write(std::cout);
    }
    if (parser.hasField("TLMIXPAR")) {
	parser.getTLMIXPAR().write(std::cout);
    }
    if (parser.hasField("PLYMAX")) { 
	parser.getPLYMAX().write(std::cout);
    }
    if (parser.hasField("WPOLYMER")) {
	parser.getWPOLYMER().write(std::cout);
    }

}


