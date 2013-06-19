//===========================================================================
//
// File: ErrorMacros.hpp
//
// Created: Tue May 14 12:22:16 2002
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

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

#ifndef OPM_ERRORMACROS_HEADER
#define OPM_ERRORMACROS_HEADER


/// Error macros. In order to use some of them, you must also
/// include <iostream> or <exception>, so they are included by
/// this file. The compile defines NDEBUG and NVERBOSE control
/// the behaviour of these macros.

#ifndef NVERBOSE
  #include <iostream>
#endif
#include <exception>

/// Usage: REPORT;
/// Usage: MESSAGE("Message string.");
#ifdef NVERBOSE // Not verbose mode
#  ifndef REPORT
#    define REPORT
#  endif
#  ifndef MESSAGE
#    define MESSAGE(x)
#  endif
#  ifndef MESSAGE_IF
#    define MESSAGE_IF(cond, m)
#  endif
#else // Verbose mode
#  ifndef REPORT
#    define REPORT std::cerr << "\nIn file " << __FILE__ << ", line " << __LINE__ << std::endl
#  endif
#  ifndef MESSAGE
#    define MESSAGE(x) std::cerr << "\nIn file " << __FILE__ << ", line " << __LINE__ << ": " << x << std::endl
#  endif
#  ifndef MESSAGE_IF
#    define MESSAGE_IF(cond, m) do {if(cond) MESSAGE(m);} while(0)
#  endif
#endif


/// Usage: THROW("Error message string.");
#ifndef THROW
#  define THROW(x) do { MESSAGE(x); throw std::exception(); } while(0)
#endif

#define ALWAYS_ERROR_IF(condition, message) do {if(condition){ THROW(message);}} while(0)

/// Usage: ASSERT(condition)
/// Usage: ASSERT2(condition, "Error message string.")
/// Usage: DEBUG_ERROR_IF(condition, "Error message string.");
#ifdef NDEBUG // Not in debug mode
#  ifndef ASSERT
#    define ASSERT(x)
#  endif
#  ifndef ASSERT2
#    define ASSERT2(cond, x)
#  endif
#  ifndef DEBUG_ERROR_IF
#    define DEBUG_ERROR_IF(cond, x)
#  endif
#else // Debug mode
#  ifndef ASSERT
#    define ASSERT(cond) if (!(cond)) THROW("Assertion \'" #cond "\' failed.")
#  endif
#  ifndef ASSERT2
#    define ASSERT2(cond, x) do { if (!(cond)) THROW(x);} while(0)
#  endif
#  ifndef DEBUG_ERROR_IF
#    define DEBUG_ERROR_IF(cond, x) do { if (cond) THROW(x); } while(0)
#  endif
#endif


#endif // OPM_ERRORMACROS_HEADER
