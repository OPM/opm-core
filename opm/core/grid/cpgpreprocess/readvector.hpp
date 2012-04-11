/*===========================================================================
//
// File: readvector.hpp
//
// Created: 2011-11-30 09:35:14+0100
//
// Author: Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
*/

#ifndef READVECTOR_HPP_HEADER
#define READVECTOR_HPP_HEADER

#include <string>

void read_vector_from_file(const std::string&, std::vector<int>& v);
void read_vector_from_file(const std::string&, std::vector<double>& v);

#endif  /* READVECTOR_HPP_HEADER */
