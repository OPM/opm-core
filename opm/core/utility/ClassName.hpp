// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
/*****************************************************************************
 *   Copyright (C) 2012-2013 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Oliver Sander                                     *
 *   Copyright (C) 2011 by Martin Nolte                                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, version 2.                                *
 *                                                                           *
 *   This file is based on works developed by the DUNE project, see          *
 *   <http://www.dune-project.org>. The following license exception          *
 *   applies to it:                                                          *
 *                                                                           *
 *   As a special exception, you may use the DUNE library without            *
 *   restriction.  Specifically, if other files instantiate templates or     *
 *   use macros or inline functions from one or more of the DUNE source      *
 *   files, or you compile one or more of the DUNE source files and link     *
 *   them with other files to produce an executable, this does not by        *
 *   itself cause the resulting executable to be covered by the GNU          *
 *   General Public License.  This exception does not however invalidate     *
 *   any other reasons why the executable file might be covered by the       *
 *   GNU General Public License.                                             *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef OPM_CLASSNAME_HH
#define OPM_CLASSNAME_HH

/** \file
 * \brief A free function to provide the demangled class name
 *        of a given object or type as a string
 */

#include <cstdlib>
#include <string>
#include <typeinfo>

#if HAVE_CXA_DEMANGLE
#include <cxxabi.h>
#endif

namespace Opm {

template <class T>
class ClassNameHelper_
{ public:
    static std::string name()
    {
        std::string className = typeid( T ).name();
#if HAVE_CXA_DEMANGLE
        int status;
        char *demangled = abi::__cxa_demangle( className.c_str(), 0, 0, &status );
        if( demangled )
        {
          className = demangled;
          std::free( demangled );
        }
#endif // #ifdef __GNUC__
        return className;
    }
};

#if HAVE_QUAD
// specialize for quad precision floating point values to avoid
// needing a type_info structure
template <>
class ClassNameHelper_<__float128>
{ public:
    static std::string name()
    { return "quad"; }
};
#endif

/** \brief Provide the demangled class name of a given object as a string */
template <class T>
std::string className()
{
    return ClassNameHelper_<T>::name();
}

/** \brief Provide the demangled class name of a given object as a string */
template <class T>
std::string className(const T &t)
{
    return ClassNameHelper_<T>::name();
}
} // namespace Opm

#endif  // OPM_CLASSNAME_HH
