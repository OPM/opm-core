/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Statoil ASA.

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

#ifndef OPM_COMPRESSEDPROPERTYACCESS_HPP_HEADER
#define OPM_COMPRESSEDPROPERTYACCESS_HPP_HEADER

/**
 * \file
 *
 * Facility for accessing active subset of data arrays defined for all
 * global cells.  The main component is class template \code
 * GridPropertyAccess::Compressed<> \endcode which encapsulates and
 * provides read-only access to a data array and while translating
 * active cell indices to "global" cell indices.  The data array is a
 * policy parameter for which preexisting implementations "constant"
 * and "extract from ECLIPSE input" are defined in this module.  Data
 * values in the array must be defined for all global cells.
 */

#include <opm/parser/eclipse/EclipseState/Grid/GridProperty.hpp>

#include <cassert>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Opm {
    /**
     * Nested name-space that serves no other purpose than to
     * contextualise \c Compressed class name.
     */
    namespace GridPropertyAccess {
        /**
         * Glue code for interacting with ECLIPSE input decks as
         * defined by module opm-parser.
         */
        namespace Details {
            /**
             * Implementation of property query and retrieval from an
             * ECLIPSE property container.
             */
            namespace EclPropImpl {
                /**
                 * Property existence predicate.
                 *
                 * Supported for types \c int and \c double.
                 *
                 * \tparam T Property element type.
                 */
                template <typename T>
                struct HasProperty;

                /**
                 * Property value retrieval.
                 *
                 * Supported for types \c int and \c double.
                 *
                 * \tparam T Property element type.
                 */
                template <typename T>
                struct GetProperty;

                /**
                 * Specialization of property existence predicate for
                 * type \c int.
                 */
                template <>
                struct HasProperty<int> {
                    /**
                     * Existence predicate implementation.
                     *
                     * \tparam PropertyContainer Pointer type
                     * representing collection of (global) grid
                     * properties.  Must implement method \c
                     * hasIntGridProperty.
                     *
                     * \param[in] ecl Property container.
                     *
                     * \param[in] kw ECLIPSE property keyword.
                     *
                     * \return Whether or not property \c kw exists in
                     * the container \c ecl.
                     */
                    template <class PropertyContainer>
                    static bool
                    p(PropertyContainer& ecl,
                      const std::string& kw);
                };

                template <class PropertyContainer>
                bool
                HasProperty<int>::p(PropertyContainer& ecl,
                                    const std::string& kw)
                {
                    return ecl->hasIntGridProperty(kw);
                }

                /**
                 * Specialization of property existence predicate for
                 * type \c double.
                 */
                template <>
                struct HasProperty<double> {
                    /**
                     * Existence predicate implementation.
                     *
                     * \tparam PropertyContainer Pointer type
                     * representing collection of (global) grid
                     * properties.  Must implement method \c
                     * hasDoubleGridProperty.
                     *
                     * \param[in] ecl Property container.
                     *
                     * \param[in] kw ECLIPSE property keyword.
                     *
                     * \return Whether or not property \c kw exists in
                     * the container \c ecl.
                     */
                    template <class PropertyContainer>
                    static bool
                    p(PropertyContainer& ecl,
                      const std::string& kw);
                };

                template <class PropertyContainer>
                bool
                HasProperty<double>::p(PropertyContainer& ecl,
                                       const std::string& kw)
                {
                    return ecl->hasDoubleGridProperty(kw);
                }

                /**
                 * Specialization of property value retrieval for type
                 * \c int.
                 */
                template <>
                struct GetProperty<int> {
                    /**
                     * Property value retrieval implementation.
                     *
                     * \tparam PropertyContainer Pointer type
                     * representing collection of (global) grid
                     * properties.  Must implement method \c
                     * getIntGridProperty.
                     *
                     * \param[in] ecl Property container.
                     *
                     * \param[in] kw ECLIPSE property keyword.
                     *
                     * \return Data values for property \c kw.
                     */
                    template <class PropertyContainer>
                    static std::shared_ptr< GridProperty<int> >
                    value(PropertyContainer& ecl,
                          const std::string& kw);
                };

                template <class PropertyContainer>
                std::shared_ptr< GridProperty<int> >
                GetProperty<int>::value(PropertyContainer& ecl,
                                        const std::string& kw)
                {
                    assert (HasProperty<int>::p(ecl, kw));

                    return ecl->getIntGridProperty(kw);
                }

                /**
                 * Specialization of property value retrieval for type
                 * \c double.
                 */
                template <>
                struct GetProperty<double> {
                    /**
                     * Property value retrieval implementation.
                     *
                     * \tparam PropertyContainer Pointer type
                     * representing collection of (global) grid
                     * properties.  Must implement method \c
                     * getDoubleGridProperty.
                     *
                     * \param[in] ecl Property container.
                     *
                     * \param[in] kw ECLIPSE property keyword.
                     *
                     * \return Data values for property \c kw.
                     */
                    template <class PropertyContainer>
                    static std::shared_ptr< GridProperty<double> >
                    value(PropertyContainer& ecl,
                          const std::string& kw);
                };

                template <class PropertyContainer>
                std::shared_ptr< GridProperty<double> >
                GetProperty<double>::value(PropertyContainer& ecl,
                                           const std::string& kw)
                {
                    assert (HasProperty<double>::p(ecl, kw));

                    return ecl->getDoubleGridProperty(kw);
                }
            } // namespace EclPropImpl

            /**
             * Conditional retrieval of property values from an
             * ECLIPSE input deck.
             *
             * Supported for types \c int and \c double.
             *
             * \tparam T Property element type.
             */
            template <typename T>
            struct EclipsePropertyArray {
                /**
                 * Retrieve property values if present in container.
                 *
                 * \tparam PropertyContainer Pointer type representing
                 * collection of (global) grid properties.
                 *
                 * \param[in] ecl Property container.
                 *
                 * \param[in] kw ECLIPSE property keyword.
                 *
                 * \return Data values for property \c kw if present,
                 * an empty \code shared_ptr<> \endcode if not.
                 */
                template <class PropertyContainer>
                static std::shared_ptr< GridProperty<T> >
                value(PropertyContainer& ecl,
                      const std::string& kw);
            };

            template <typename T>
            template <class PropertyContainer>
            std::shared_ptr< GridProperty<T> >
            EclipsePropertyArray<T>::value(PropertyContainer& ecl,
                                           const std::string& kw)
            {
                std::shared_ptr< GridProperty<T> > x;

                if (EclPropImpl::HasProperty<T>::p(ecl, kw)) {
                    x = EclPropImpl::GetProperty<T>::value(ecl, kw);
                }

                return x;
            }
        } // namespace Details

        /**
         * Predefined data array policies for use with class template
         * \code Compressed<> \endcode.
         */
        namespace ArrayPolicy {
            /**
             * Data array policy that extracts the array values from
             * an ECLIPSE input deck or returns a user specified
             * default value if the data vector is not present in a
             * particular input deck.
             *
             * Provides read-only access to the underlying data.
             *
             * \tparam T Array element type.  Must be \c int or \c
             * double.
             */
            template <typename T>
            class ExtractFromDeck {
            public:
                /**
                 * Constructor.
                 *
                 * \tparam PropertyContainer Pointer type representing
                 * collection of (global) grid properties.  Typically
                 * \c EclipseStatePtr or \c EclipseStateConstPtr.
                 * Must implement methods \c hasIntGridProperty and \c
                 * getIntGridProperty if \c T is \c int, or \c
                 * hasDoubleGridProperty and \c getDoubleGridProperty
                 * if \c T is \c double.
                 *
                 * \param[in] ecl Property container.
                 *
                 * \param[in] kw ECLIPSE keyword from which to extract
                 * data array.
                 *
                 * \param[in] dlft Default/fall-back data array value
                 * if \c kw is not defined.
                 */
                template <class PropertyContainer>
                ExtractFromDeck(PropertyContainer& ecl,
                                const std::string& kw,
                                const T            dflt)
                    : x_   (Details::EclipsePropertyArray<T>::value(ecl, kw))
                    , dflt_(dflt)
                {}

                /**
                 * Publicly accessible data array element type.
                 */
                typedef T value_type;

                /**
                 * Index type for accessing data array.
                 */
                typedef std::size_t size_type;

                /**
                 * Read-only data array access.
                 *
                 * \param[in] i Array index.  Assumed to identify a
                 * global (uncompressed) cell.
                 *
                 * \return Data array element at global index \c i if
                 * present in input or user specified fall-back value
                 * if not.
                 */
                value_type
                operator[](const size_type i) const
                {
                    if (x_) {
                        return x_->iget(i);
                    }
                    else {
                        return dflt_;
                    }
                }

            private:
                /**
                 * Grid property handle.
                 *
                 * Null if data not defined.
                 */
                std::shared_ptr< GridProperty<T> > x_;

                /**
                 * Fall-back data element value if data not defined.
                 */
                T dflt_;
            };

            /**
             * Data array policy that returns a single, constant user
             * specified value for every global cell.
             *
             * Provides read-only access to the underlying data.
             *
             * \tparam T Array element type.
             */
            template <typename T>
            class Constant {
            public:
                /**
                 * Constructor
                 *
                 * \param[in] c Constant property value used for all
                 * global cells.
                 */
                Constant(const T c)
                    : c_(c)
                {}

                /**
                 * Publicly accessible data array element type.
                 */
                typedef T value_type;

                /**
                 * Index type for accessing data array.
                 */
                typedef std::size_t size_type;

                /**
                 * Read-only data array access.
                 *
                 * \param[in] i Array index.  Assumed to identify a
                 * global (uncompressed) cell.  Unused.
                 *
                 * \return User specified constant value for every
                 * (global) cell.
                 */
                value_type
                operator[](const size_type i) const
                {
                    static_cast<void>(i); // Suppress "unused parameter"

                    return c_;
                }

            private:
                /**
                 * Constant, user specified property value.
                 */
                T c_;
            };
        } // namespace ArrayPolicy

        /**
         * Collection of tags to help enforce semantic type checks
         * when using class \code Compressed<> \endcode.
         */
        namespace Tag {
            /**
             * Default tag that implies no restriction.
             */
            struct Any {};

            /**
             * Tag that restricts usage to NTG (net-to-gross)
             * contexts.
             */
            struct NTG : public Any {};
        } // namespace Tag

        /**
         * Provide compressed (active cell) read-only access to
         * globally defined data array.
         *
         * \tparam DataArray Type representing an array of data
         * values, one value for each global (uncompressed) cell in a
         * model.  Must implement value semantics.  Typically one of
         * the array policies of name space \c ArrayPolicy.  Must
         * provide public type \c value_type to infer the data element
         * type and \code operator[](i) \endcode to access the
         * property value of the \c i'th global cell.
         *
         * \tparam PropertyTag Type tag that can be used to restrict
         * applicability of the resulting \c Compressed array, e.g.,
         * to enforce net-to-gross ratios only.  Default: No
         * restriction.
         */
        template <class DataArray, class PropertyTag = Tag::Any>
        class Compressed {
        public:
            /**
             * Constructor
             *
             * \param[in] x Preconfigured global property value array.
             * The \c Compressed array creates a private copy of this
             * object.
             *
             * \param[in] gc Compressed-to-global cell map.  Typically
             * the \c global_cell field of an \c UnstructuredGrid or
             * something very similar.  If null, interpreted as
             * identity mapping, i.e., as if all cells are active.
             */
            Compressed(const DataArray& x,
                       const int*       gc)
                : x_ (x)
                , gc_(gc)
            {}

            /**
             * Property value type.
             */
            typedef typename DataArray::value_type value_type;

            /**
             * Read-only data array access.
             *
             * \param[in] c Active cell index.
             *
             * \return Property value in active cell \c c.
             */
            value_type
            operator[](const int c) const
            {
                return x_[ (gc_ == 0) ? c : gc_[c] ];
            }

        private:
            /**
             * Global property value array.
             *
             * Value semantics to support putting \c Compressed arrays
             * into standard containers.
             */
            DataArray x_;

            /**
             * Compressed-to-global cell index map.  \c Null if all
             * cells active.
             */
            const int* gc_;
        };
    } // namespace GridPropertyAccess
} // namespace Opm

#endif  /* OPM_COMPRESSEDPROPERTYACCESS_HPP_HEADER */
