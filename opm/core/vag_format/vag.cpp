#include <opm/core/vag_format/vag.hpp>
#include <iostream>
#include <fstream>
#include <string>
namespace OPM
{
    void readPosStruct(std::istream& is,int n,PosStruct& pos_struct){
	using namespace std;
	//PosStruct pos_struct;
	pos_struct.pos.resize(n+1);
	pos_struct.pos[0]=0;
	for(int i=0;i< n;++i){
	    int number;
	    is >> number ;
	    //cout <<number << endl;
	    pos_struct.pos[i+1]=pos_struct.pos[i]+number;
	    for(int j=0;j< number;++j){
		int value;
		is >> value;
		//	cout << value << " ";
		pos_struct.value.push_back(value);
	    }
	    //cout << endl;
	}
	if(!(pos_struct.value.size()==pos_struct.pos[n])){
	    cerr << "Failed to read pos structure" << endl;
	    cerr << "pos_struct.value.size()" << pos_struct.value.size() << endl;
	    cerr << "pos_struct.pos[n+1]" << pos_struct.pos[n] << endl;
	}
    };
    void readVagGrid(std::istream& is,OPM::VAG& vag_grid){
	using namespace std;
	using namespace OPM;
	while (!is.eof()) {
	    string keyword;
	    is >> keyword;	
	    //cout << keyword<< endl;
	    if(keyword == "Number"){
		string stmp;
		is >> stmp;
		if(stmp == "of"){
		    string entity;
		    is >> entity;
		    getline(is,stmp);
		    int number;
		    is >> number;
		    if(entity=="vertices"){
			vag_grid.number_of_vertices=number;
		    }else if((entity=="volumes") || (entity=="control")){
			vag_grid.number_of_volumes=number;
		    }else if(entity=="faces"){
			vag_grid.number_of_faces=number;
		    }else if(entity=="edges"){
			vag_grid.number_of_edges=number;
		    }	
		    cout << "Found Number of: " << entity <<" " << number << endl;		
		} else {
		    cerr << "Wrong format: Not of after Number" << endl;
		    return;
		}
	    }else{
		// read geometry defined by vertices
		if(keyword=="Vertices"){
		    int number;
		    is >> number;	   	    
		    vag_grid.vertices.resize(3*number);// assume 3d data
		    readVector(is,vag_grid.vertices);
		}
		// here starts the reding of all pos structures
		else if(keyword=="Volumes->Faces" || keyword=="Volumes->faces"){
		    //vag_grid.volumes_to_faces=
		    int number;
		    is >> number;	
		    readPosStruct(is,number,vag_grid.volumes_to_faces);
		    cout << "Volumes->Faces: Number of " << number << endl;
		}else if(keyword=="Faces->edges" || keyword=="Faces->Edges" ||  keyword=="Faces->Edgess"){
		    int number;
		    is >> number;	   
		    //vag_grid.volumes_to_faces=
		    readPosStruct(is,number,vag_grid.faces_to_edges);
		    cout << "Faces->edges: Number of " << number << endl;
		}else if(keyword=="Faces->Vertices" || keyword=="Faces->vertices"){
		    int number;
		    is >> number;	   
		    //vag_grid.volumes_to_faces=
		    readPosStruct(is,number,vag_grid.faces_to_vertices);
		    cout << "Faces->Vertices: Number of " << number << endl;   
		}else if(keyword=="Volumes->Vertices" || keyword=="Volumes->Verticess"){
		    int number;
		    is >> number;	   
		    //vag_grid.volumes_to_faces=
		    readPosStruct(is,number,vag_grid.volumes_to_vertices);
		    cout << "Volumes->Vertices: Number of " << number << endl;   
		}
	    
		//  read simple mappings
		else if(keyword=="Edge" || keyword=="Edges"){
		    int number;
		    is >> number;	   	    
		    vag_grid.edges.resize(2*number);
		    readVector(is,vag_grid.edges);
		    cout << "Edges: Number of " << number << endl;
		}else if(keyword=="Faces->Volumes" || keyword=="Faces->Control"){
		    int number;
		    if(keyword=="Faces->Control"){
			string vol;
			is >> vol;
		    }
		    is >> number;
		    vag_grid.faces_to_volumes.resize(2*number);
		    readVector(is,vag_grid.faces_to_volumes);
		    cout << "Faces->Volumes: Number of " << number << endl;   
		}
		// read material
		else if(keyword=="Material"){
		    string snum;
		    is >> snum;
		    int number;
		    is >> number;
		    cout << "Material number  " << number << endl;
		    // we read all the rest into doubles		
		    while(!is.eof()){
			double value;
			is >> value;
			//cout << value << endl;
			vag_grid.material.push_back(value);
		    }
		}else{
		    //cout << "keyword;
		}
		//cout << "Found" << keyword << "Number of " << number << endl;
	    }
	}
    };
}
