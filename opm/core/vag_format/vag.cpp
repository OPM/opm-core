/*===========================================================================
//
// File: vag.cpp
//
// Created: 2012-06-08 15:45:53+0200
//
// Authors: Knut-Andreas Lie      <Knut-Andreas.Lie@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Xavier Raynaud        <Xavier.Raynaud@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <opm/core/vag_format/vag.hpp>
#include <opm/core/grid/cornerpoint_grid.h>	 	      
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <cassert>
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
	if(!(int(pos_struct.value.size())==pos_struct.pos[n])){
	    cerr << "Failed to read pos structure" << endl;
	    cerr << "pos_struct.value.size()" << pos_struct.value.size() << endl;
	    cerr << "pos_struct.pos[n+1]" << pos_struct.pos[n] << endl;
	}
    }
    void writePosStruct(std::ostream& os,PosStruct& pos_struct){
	using namespace std;
	//PosStruct pos_struct;
        int n=pos_struct.pos.size()-1;
	pos_struct.pos.resize(n+1);
	pos_struct.pos[0]=0;
	for(int i=0;i< n;++i){
	    int number=pos_struct.pos[i+1]-pos_struct.pos[i];
	    os << number ;
	    for(int j=0;j< number;++j){
		os << pos_struct.value[pos_struct.pos[i]+j];
	    }
	    os << endl;
	}
    }
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
    }
    void vagToUnstructuredGrid(OPM::VAG& vag_grid,UnstructuredGrid& grid){
	using namespace std;
	using namespace OPM;
	cout << "Converting grid" << endl;
	grid.dimensions=3;
	grid.number_of_cells=vag_grid.number_of_volumes;
	grid.number_of_faces=vag_grid.number_of_faces;
	grid.number_of_faces=vag_grid.number_of_faces;
	// fill face_nodes
	for(int i=0;i< int(vag_grid.faces_to_vertices.pos.size());++i){
	    grid.face_nodepos[i] = vag_grid.faces_to_vertices.pos[i];
	}	    
	for(int i=0;i< int(vag_grid.faces_to_vertices.value.size());++i){
	    grid.face_nodes[i] = vag_grid.faces_to_vertices.value[i]-1;
	}
	// fill cell_face
	for(int i=0;i< int(vag_grid.volumes_to_faces.pos.size());++i){
	    grid.cell_facepos[i] = vag_grid.volumes_to_faces.pos[i];
	}	    
	for(int i=0;i< int(vag_grid.volumes_to_faces.value.size());++i){
	    grid.cell_faces[i] = vag_grid.volumes_to_faces.value[i]-1;
	}
	// fill face_cells
	for(int i=0;i< int(vag_grid.faces_to_volumes.size());++i){
	    grid.face_cells[i] = vag_grid.faces_to_volumes[i]-1;
	}

	// fill node_cordinates. This is the only geometry given in the vag
	for(int i=0;i< int(vag_grid.vertices.size());++i){
	    grid.node_coordinates[i] = vag_grid.vertices[i];
	}
	// informations in edges, faces_to_eges, faces_to_vertices, volume_to_vertices and materials
	// is not used
	cout << "Computing geometry" << endl;
	compute_geometry(&grid);
	
    }
    void writeVagFormat(std::ostream& os,OPM::VAG& vag_grid){
	using namespace std;
	os << "File in the Vag grid format";
        os << "Number of vertices" << endl;
        os << vag_grid.number_of_vertices;
        os <<"Number of control volume " << endl;
        os << vag_grid.number_of_volumes;
        os <<"Number of faces" << endl;
        os << vag_grid.number_of_faces;
        os <<"Number of edges" << endl;
        os << vag_grid.number_of_edges;
        os <<"Vertices      "   << vag_grid.vertices.size() << endl;
        writeVector(os, vag_grid.vertices,3);
        os << "Volumes->faces   " << vag_grid.volumes_to_faces.pos.size()-1 << endl;
        writePosStruct(os, vag_grid.volumes_to_faces);
        os << "Volumes->Vertices   " << vag_grid.volumes_to_vertices.pos.size()-1 << endl;
        writePosStruct(os, vag_grid.volumes_to_vertices);
        os << "Faces->edges   " << vag_grid.faces_to_edges.pos.size()-1 << endl;
        writePosStruct(os, vag_grid.faces_to_edges);
        os << "Faces->Control volumes   " << vag_grid.faces_to_volumes.size() << endl;
        writeVector(os,vag_grid.faces_to_volumes,2);
        os << "Edges   " << vag_grid.edges.size() << endl;
        writeVector(os,vag_grid.edges,2);
        /*
        assert(vag_grid.material.size()%vag_grid.number_of_volumes==0);
        int lines= floor(vag_grid.material.size()/vag_grid.number_of_volumes);  
        os << "Material number   " << 1 << endl;
        writeVector(os,vag_grid.material,lines);
        */

    }


}
