//===========================================================================
//
// File: Factory.hpp
//
// Created: Mon Nov  6 10:00:43 2000
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

#ifndef OPM_FACTORY_HEADER
#define OPM_FACTORY_HEADER

#include <map>
#include <memory>

namespace Opm
{

    /** This is an object factory for creating objects of some type
     *  requested by the user, with a shared base class.  The user
     *  need only interact with the factory through the static
     *  template member addCreator() and the static member function
     *  createObject().
     */
    template<class Base>
    class Factory
    {
    public:
        /// The type of pointer returned by createObject().
        typedef std::shared_ptr<Base> ProductPtr;

        /// Creates a new object of the class associated with the given type string,
        /// and returns a pointer to it.
        /// \param type the type string of the class that the user wants to have
        ///             constructed.
        /// \return (smart) pointer to the created object.
        static ProductPtr createObject(const std::string& type)
        {
            return instance().doCreateObject(type);
        }

        /// Clones an new object of the class associated with the given type string,
        /// and returns a pointer to it.
        /// \param type the type string of the class that the user wants to have
        ///             constructed.
        /// \param original (smart) pointer to object to be cloned.
        /// \return (smart) pointer to the created object.
        static ProductPtr cloneObject(const std::string& type,
                                      const ProductPtr original)
        {
            return instance().doCloneObject(type, original);
        }

        /// Add a creator to the Factory.
        /// After the call, the user may obtain new objects of the Derived type by
        /// calling createObject() with the given type string as an argument.
        /// \tparam Derived the class we want to add a creator for, must inherit
        ///                 the class template parameter Base.
        /// \param type the type string with which we want the Factory to associate
        ///             the class Derived.
        template<class Derived>
        static void addCreator(const std::string& type)
        {
            instance().template doAddCreator<Derived>(type);
        }

    private:
        // The method that implements the singleton pattern,
        // using the Meyers singleton technique.
        static Factory& instance()
        {
            static Factory singleton;
            return singleton;
        }

        // Private constructor, to keep users from creating a Factory.
        Factory()
        {
        }

        // Abstract base class for Creators.
        class Creator
        {
        public:
            virtual ProductPtr create() = 0;
            virtual ProductPtr clone(ProductPtr original) = 0;
            virtual ~Creator() {}
        };

        /// This is the concrete Creator subclass for generating Derived objects.
        template<class Derived>
        class ConcreteCreator: public Creator
        {
        public:
            virtual ProductPtr create()
            {
                return ProductPtr(new Derived);
            }
            virtual ProductPtr clone(const ProductPtr original)
            {
                const Derived& orig = dynamic_cast<const Derived&>(*original);
                return ProductPtr(new Derived(orig));
            }
        };

        typedef std::shared_ptr<Creator> CreatorPtr;
        typedef std::map<std::string, CreatorPtr> CreatorMap;
        // This map contains the whole factory, i.e. all the Creators.
        CreatorMap string_to_creator_;

        // Actually creates the product object.
        ProductPtr doCreateObject(const std::string& type)
        {
            typename CreatorMap::iterator it;
            it = string_to_creator_.find(type);
            if (it == string_to_creator_.end()) {
                THROW("Creator type " << type
                      << " is not registered in the factory.");
            }
            return it->second->create();
        }

        // Actually clones the product object.
        ProductPtr doCloneObject(const std::string& type,
                                 const ProductPtr original)
        {
            typename CreatorMap::iterator it;
            it = string_to_creator_.find(type);
            if (it == string_to_creator_.end()) {
                THROW("Creator type " << type
                      << " is not registered in the factory.");
            }
            return it->second->clone(original);
        }

        // Actually adds the creator.
        template<class Derived>
        void doAddCreator(const std::string& type)
        {
            CreatorPtr c(new ConcreteCreator<Derived>);
            string_to_creator_[type] = c;
        }
    };

} // namespace Opm

#endif // OPM_FACTORY_HEADER
