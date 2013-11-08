/*
  Copyright (c) 2013 Uni Research AS

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

#ifndef OPM_SHARED_OBJECT_HPP
#define OPM_SHARED_OBJECT_HPP

#include <memory>  // shared_ptr

namespace Opm {

/**
 * @see Opm::share_obj
 */
template <typename T>
class SharedObject : public std::shared_ptr <T> {
public:
    SharedObject (T& obj) : std::shared_ptr <T> (&obj, no_delete) { }
private:
    /// Custom deleter that does nothing
    static void no_delete (void const *) { }
};

/*!
 * Share pointer of a local object.
 *
 * Use this wrapper when an interface needs a shared_ptr, but you
 * want to pass an object that has local storage (and you know
 * that the shared_ptr client doesn't need it outside of the scope).
 *
 * \example
 * \code{.cpp}
 *  Foo obj;
 *  std::shared_ptr <Foo> ptr = share_obj (obj);
 * \endcode
 */
template <typename T> std::shared_ptr <T> share_obj (T& t) {
    return SharedObject <T> (t);
}

} // namespace Opm

#endif /* OPM_SHARED_OBJECT_HPP */
