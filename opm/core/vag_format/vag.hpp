#ifndef OPM_VAG_HEADER_INCLUDED
#define OPM_VAG_HEADER_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
namespace OPM
{
    /* Struct to hold maping from the natural number less than pos.size()-1 to
       a set of integers. value(pos(i):pos(i+1)-1) hold the integers corresponding to i.
       pos(end)-1==value.size();
    */               
    struct PosStruct{
	std::vector<int> pos;
        std::vector<int>  value;
    };
    /* Structure to represent the unstructured vag grid format
     */
    struct VAG{
	int number_of_vertices;
	int number_of_volumes;
	int number_of_faces;
	int number_of_edges;
	std::vector<double> vertices;
	PosStruct volumes_to_faces;
	PosStruct volumes_to_vertices;
	PosStruct faces_to_edges;
	PosStruct faces_to_vertices;
	std::vector<int> edges;
	std::vector<int> faces_to_volumes;
	std::vector<double> material;
    };
    /* Function the vag grid format and make a vag_grid struct. This structure
       is intended to be converted to a grid*/
    void readVagGrid(std::istream& is,OPM::VAG& vag_grid);
    /*
    void writeVagFormat(std::ostream& os){
	using namespace std;
	os << "File in the Vag grid format" << endl;
    };
    */
    template <typename T>
    void readVector(std::istream& is,std::vector<T>& vec){
	using namespace std;
	for(int i=0;i< vec.size();++i){
	    is >> vec[i];
	}
    }
       
    //PosStruct readPosStruct(std::istream& is,int n){
    void readPosStruct(std::istream& is,int n,PosStruct& pos_struct);	
}
	
#endif // OPE    
