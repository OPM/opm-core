#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <opm/core/vag_format/vag.hpp>
//#include "../config.h"
/* test reading of vag grid format */
int main(int argc, char** argv)
{
    using namespace std;
    using namespace OPM;
    std::string filename;
    if (argc == 2) {
	filename = argv[1];
    } else {
	std::cout << "\nUsage: test_read_vag  filename\n";
	exit( 1 );
    }
    ifstream is(filename.c_str());//"/home/hnil/heim/SVN/simmatlab/projects/clastic/utils/unstructuredgrids/data/3x3_w_layered-vag.dat");
    //ifstream is("/home/hnil/heim/SVN/simmatlab/projects/clastic/utils/unstructuredgrids/data/test.txt");
    //std::ofstream is("");
    VAG vag_grid;
    readVagGrid(is,vag_grid);
}
