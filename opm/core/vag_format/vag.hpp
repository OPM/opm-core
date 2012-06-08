/*===========================================================================
//
// File: vag.hpp
//
// Created: 2012-06-08 15:46:23+0200
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

#ifndef OPM_VAG_HPP_HEADER
#define OPM_VAG_HPP_HEADER




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
#endif  /* OPM_VAG_HPP_HEADER */	

