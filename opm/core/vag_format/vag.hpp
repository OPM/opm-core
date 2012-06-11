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



#include <opm/core/grid.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
namespace OPM
{
    /**
       Struct to hold maping from the natural number less than pos.size()-1 to
       a set of integers. value(pos(i):pos(i+1)-1) hold the integers corresponding to i.
       pos(end)-1==value.size();
    */               
    struct PosStruct{
	std::vector<int> pos;
        std::vector<int>  value;
    };
    /**
      Structure to represent the unstructured vag grid format. The format is only for
      3D grids.
     */
    struct VAG{
	int number_of_vertices;        
	int number_of_volumes;
	int number_of_faces;
	int number_of_edges;
        /** Vertices. The coordinates of vertice i is [vetices[3*i:3*i+2]*/
	std::vector<double> vertices;
        /** Mapping from volumes to faces */
	PosStruct volumes_to_faces;
        /** Mapping from volumes to vertices */
	PosStruct volumes_to_vertices;
        /** Mapping from faces to edges */
	PosStruct faces_to_edges;
        /** Mapping from faces to vertices */
	PosStruct faces_to_vertices;
        /** The edge i is given by the nodes edges[2*i:2*i+1] */
	std::vector<int> edges;
        /**  The two neigbours of the face i is faces_to_volumes[2*i:2*i+1] */
	std::vector<int> faces_to_volumes;
        /** A vector containing information of each volume. The size is n*number_of_volumes.
            For each i this is the information:
            material[n*i] is the volume number and should be transformed to integer
            material[n*i+1] is a tag and should  be transformed to integer
            material[n*i+2:n*(i+1)-1] represent propertices.
        */
	std::vector<double> material;
    };
    /**
       Function the vag grid format and make a vag_grid struct. This structure
       is intended to be converted to a grid.
       \param[in]  is is is stream of the file.
       \param[out] vag_grid is a reference to a vag_grid struct.
    */
    void readVagGrid(std::istream& is,OPM::VAG& vag_grid);
    /**
       Function to write vag format.
       \param[out]  is is is stream of the file.
       \param[in] vag_grid is a reference to a vag_grid struct.

    */
    void writeVagFormat(std::ostream& os,OPM::VAG& vag_grid);
    /**
       Function to read a vector of some type from a stream.
       \param[in]  os is is stream of the file.
       \param[out] vag_grid is a resized and filled vector containing the quantiy read.
    */   
    template <typename T>
    void readVector(std::istream& is,std::vector<T>& vec){
	using namespace std;
	for(int i=0;i< int(vec.size());++i){
	    is >> vec[i];
	}
    }
    /**
       Function to write a vector of some type from a stream.
       \param[in]  os is is stream of the file.
       \param[out] vag_grid is a resized and filled vector containing the quantiy read.
       \param[in]  n number of doubles on each line.
       The function will only write full lines an potentially skip numbers at end of the vector.
    */  
    template <typename T>
    void writeVector(std::ostream& os,std::vector<T>& vec,int n){
	using namespace std;
        int lines = floor(vec.size()/n);
        assert(vec.size()%n==0);
	for(int j=0;j< lines;++j){
            for(int i=0;i< n;++i){
	    os << vec[j*n+i];
            }
        }
    }

    /**
       Read pos struct type mapping from a stream
       \param[in] is is stream
       \param[in] n number of lines to read
       \param[out] pos_struct reference to PosStruct
     */
    void readPosStruct(std::istream& is,int n,PosStruct& pos_struct);
    /**
       Read pos struct type mapping from a stream
       \param[in] os is stream to write to
       \param[in] pos_struct to write
     */
    void writePosStruct(std::ostream& os,PosStruct& pos_struct);
    /**
       Fill a UnstructuredGrid from a vag_grid.
       \param[out] vag_grid s is a valid vag_grid struct.
       \param[in] grid is a grid with have allocated correct size to each pointer.
     */

    void vagToUnstructuredGrid(OPM::VAG& vag_grid,UnstructuredGrid& grid);
}
#endif  /* OPM_VAG_HPP_HEADER */	

