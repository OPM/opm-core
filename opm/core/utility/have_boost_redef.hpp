/*===========================================================================
//
// File: have_boost_redef.hpp
//
// Created: 2012-07-12 11:02:22+0200
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

// Note: This file should be #include<>d after <config.h> in software
// that uses the ISTL.  This is more or less a hack to work around
// Autoconfs that do not correctly process dune-istl's `ENABLE_BOOST'
// feature.  We will happily remove this file if a better solution
// presents.

#ifndef OPM_HAVE_BOOST_REDEF_HPP_HEADER
#define OPM_HAVE_BOOST_REDEF_HPP_HEADER

#if defined(HAVE_BOOST)
#undef  HAVE_BOOST
#define HAVE_BOOST OPM_HAVE_BOOST
#endif

#endif  /* OPM_HAVE_BOOST_REDEF_HPP_HEADER */
