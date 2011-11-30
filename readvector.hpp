/*===========================================================================
//
// File: readvector.hpp
//
// Created: 2011-11-30 09:35:14+0100
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
*/

#ifndef READVECTOR_HPP_HEADER
#define READVECTOR_HPP_HEADER


template <typename T>
void read_vector_from_file(const char *fn, std::vector<T>& v);

#endif  /* READVECTOR_HPP_HEADER */
