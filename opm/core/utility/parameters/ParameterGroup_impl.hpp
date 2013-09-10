//===========================================================================
//
// File: ParameterGroup_impl.hpp
//
// Created: Tue Jun  2 19:06:46 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPM_PARAMETERGROUP_IMPL_HEADER
#define OPM_PARAMETERGROUP_IMPL_HEADER

#include <iostream>
#include <string>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/parameters/ParameterStrings.hpp>
#include <opm/core/utility/parameters/ParameterTools.hpp>
#include <opm/core/utility/parameters/Parameter.hpp>

namespace Opm {
    namespace parameter {

	template<>
	struct ParameterMapItemTrait<ParameterGroup> {
	    static ParameterGroup
            convert(const ParameterMapItem& item,
                    std::string& conversion_error)
            {
		std::string tag = item.getTag();
		if (tag != ID_xmltag__param_grp) {
		    conversion_error = "The XML tag was '" + tag +
                                       "' but should be '" +
                                       ID_xmltag__param_grp + "'.\n";
		    return ParameterGroup("", 0);
		}
		conversion_error = "";
		const ParameterGroup& pg = dynamic_cast<const ParameterGroup&>(item);
		return pg;
	    }
	    static std::string type() {return "ParameterGroup";}
	};

	namespace {
	    template <typename T>
	    inline std::string
	    to_string(const T& val)
	    {
		std::ostringstream os;
		os << val;
		return os.str();
	    }

	    inline std::string
	    to_string(const bool b) {
		if (b) {
		    return ID_true;
		} else {
		    return ID_false;
		}
	    }

	    inline std::string
	    to_string(const ParameterGroup&)
	    {
		return std::string("<parameter group>");
	    }

	    std::pair<std::string, std::string>
            filename_split(const std::string& filename)
            {
		int fpos = filename.rfind('.');
		std::string name = filename.substr(0, fpos);
		std::string type = filename.substr(fpos+1);
		return std::make_pair(name, type);
	    }

	}



        template <typename StringArray>
	ParameterGroup::ParameterGroup(int argc, StringArray argv, bool verify_syntax)
            : path_(ID_path_root), parent_(0), output_is_enabled_(true)
	{
	    if (verify_syntax && (argc < 2)) {
		std::cerr << "Usage: " << argv[0] << " "
			  << "[paramfilename1.{xml,param}] "
			  << "[paramfilename2.{xml,param}] "
			  << "[overridden_arg1=value1] "
			  << "[overridden_arg2=value2] "
			  << "[...]" << std::endl;
		exit(EXIT_FAILURE);
	    }
	    this->parseCommandLineArguments(argc, argv);
	}

        template <typename StringArray>
	void ParameterGroup::parseCommandLineArguments(int argc, StringArray argv)
        {
	    std::vector<std::string> files;
	    std::vector<std::pair<std::string, std::string> > assignments;
	    for (int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		int fpos = arg.find(ID_delimiter_assignment);
		if (fpos == int(std::string::npos)) {
		    std::string filename = arg.substr(0, fpos);
		    files.push_back(filename);
		    continue;
		}
		int pos = fpos + ID_delimiter_assignment.size();
		int spos = arg.find(ID_delimiter_assignment, pos);
		if (spos == int(std::string::npos)) {
		    std::string name = arg.substr(0, fpos);
		    std::string value = arg.substr(pos, spos);
		    assignments.push_back(std::make_pair(name, value));
		    continue;
		}
		std::cout << "WARNING: Too many assignements  (' "
                          << ID_delimiter_assignment
                          << "') detected in argument " << i << ".\n";
	    }
	    for (int i = 0; i < int(files.size()); ++i) {
		std::pair<std::string, std::string> file_type = filename_split(files[i]);
		if (file_type.second == "xml") {
		    this->readXML(files[i]);
		} else if (file_type.second == "param") {
		    this->readParam(files[i]);
		} else {
		    std::cout << "WARNING: Ignoring file '"
                              << files[i] << "' with unknown extension.\n"
			      << "Valid filename extensions are 'xml' and 'param'.\n";
		}
	    }
	    for (int i = 0; i < int(assignments.size()); ++i) {
		this->insertParameter(assignments[i].first, assignments[i].second);
	    }
	}


	template<typename T>
	inline T ParameterGroup::get(const std::string& name) const
        {
	    return this->get<T>(name, ParameterRequirementNone());
	}

	template<typename T, class Requirement>
	inline T ParameterGroup::get(const std::string& name,
				     const Requirement& r) const
        {
	    setUsed();
	    std::pair<std::string, std::string> name_path = split(name);
	    map_type::const_iterator it = map_.find(name_path.first);
	    if (it == map_.end()) {
		if (parent_ != 0) {
		    // If we have a parent, ask it instead.
		    if (output_is_enabled_) {
			//TermColors::Red();
			std::cout << name;
			//TermColors::Normal();
			std::cout << " not found at "
                                  << (path() + ID_delimiter_path)
                                  << ", asking parent." << std::endl;
		    }
		    return parent_->get<T>(name, r);
		} else {
		    // We are at the top, name has not been found.
		    std::cerr << "ERROR: The group '"
			      << this->path()
			      << "' does not contain an element named '"
			      << name
			      << "'.\n";
		    throw NotFoundException();
		}
	    }
	    if (name_path.second == "") {
		T val = this->translate<T>(*it, r);
		it->second->setUsed();
		if (output_is_enabled_) {
		    //TermColors::Green();
		    std::cout << name;
		    //TermColors::Normal();
		    std::cout << " found at " << (path() + ID_delimiter_path)
			      << ", value is " << to_string(val) << std::endl;
		}
		return val;
	    } else {
		ParameterGroup& pg = dynamic_cast<ParameterGroup&>(*(*it).second);
		pg.setUsed();
		return pg.get<T>(name_path.second, r);
	    }
	}

	template<typename T>
	inline T ParameterGroup::getDefault(const std::string& name,
					    const T& default_value) const
        {
	    return this->getDefault<T>(name, default_value, ParameterRequirementNone());
	}

	template<typename T, class Requirement>
	inline T ParameterGroup::getDefault(const std::string& name,
					    const T& default_value,
					    const Requirement& r) const
        {
	    setUsed();
	    std::pair<std::string, std::string> name_path = split(name);
	    map_type::const_iterator it = map_.find(name_path.first);
	    if (it == map_.end()) {
		if (parent_ != 0) {
		    // If we have a parent, ask it instead.
		    if (output_is_enabled_) {
			//TermColors::Red();
			std::cout << name;
			//TermColors::Normal();
			std::cout <<  " not found at " << (path() + ID_delimiter_path)
				  << ", asking parent." << std::endl;
		    }
		    return parent_->getDefault<T>(name, default_value, r);
		} else {
		    // We check the requirement for the default value
		    std::string requirement_result = r(default_value);
		    if (requirement_result != "") {
			std::cerr << "ERROR: The default value for the "
				  << " element named '"
				  << name
				  << "' in the group '"
				  << this->path()
				  << "' failed to meet a requirenemt.\n";
			std::cerr << "The requirement enforcer returned the following message:\n"
				  << requirement_result
				  << "\n";
			throw RequirementFailedException<Requirement>();
		    }
		}
		if (output_is_enabled_) {
		    //TermColors::Blue();
		    std::cout << name;
		    //TermColors::Normal();
		    std::cout << " not found. Using default value '"
                              << to_string(default_value) << "'." << std::endl;
		}
		return default_value;
	    }
	    if (name_path.second == "") {
		T val = this->translate<T>(*it, r);
		it->second->setUsed();
		if (output_is_enabled_) {
		    //TermColors::Green();
		    std::cout << name;
		    //TermColors::Normal();
		    std::cout << " found at " << (path() + ID_delimiter_path)
			      << ", value is '" << to_string(val) << "'." << std::endl;
		}
		return val;
	    } else {
		ParameterGroup& pg = dynamic_cast<ParameterGroup&>(*(*it).second);
		pg.setUsed();
		return pg.getDefault<T>(name_path.second, default_value, r);
	    }
	}

	template<typename T, class Requirement>
	inline T ParameterGroup::translate(const pair_type& named_data,
					   const Requirement& chk) const
        {
	    const std::string& name = named_data.first;
	    const data_type data = named_data.second;
	    std::string conversion_error;
	    T value = ParameterMapItemTrait<T>::convert(*data, conversion_error);
	    if (conversion_error != "") {
		std::cerr << "ERROR: Failed to convert the element named '"
			  << name
			  << "' in the group '"
			  << this->path()
			  << "' to the type '"
			  << ParameterMapItemTrait<T>::type()
			  << "'.\n";
		std::cerr << "The conversion routine returned the following message:\n"
			  << conversion_error
			  << "\n";
		throw WrongTypeException();
	    }
	    std::string requirement_result = chk(value);
	    if (requirement_result != "") {
		std::cerr << "ERROR: The element named '"
			  << name
			  << "' in the group '"
			  << this->path()
			  << "' of type '"
			  << ParameterMapItemTrait<T>::type()
			  << "' failed to meet a requirenemt.\n";
		std::cerr << "The requirement enforcer returned the following message:\n"
			  << requirement_result
			  << "\n";
		throw RequirementFailedException<Requirement>();
	    }
	    return value;
	}
    } // namespace parameter
} // namespace Opm

#endif // OPM_PARAMETERGROUP_IMPL_HEADER
