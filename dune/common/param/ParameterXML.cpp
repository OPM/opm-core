//===========================================================================
//
// File: ParameterXML.cpp
//
// Created: Tue Jun  2 18:57:25 2009
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

This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/param/ParameterXML.hpp>
#include <dune/common/param/Parameter.hpp>
#include <dune/common/param/ParameterStrings.hpp>

#include <exception>
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>

#include <libxml/parser.h> // libxml2
#include <libxml/tree.h>   // libxml2


namespace Dune {
    namespace parameter {

	namespace libxml2 {
	    void read_xml(ParameterGroup& pg, const std::string filename);
	    void fill_tree(ParameterGroup& pg, xmlNodePtr root,
                           const std::string& xml_file_path);
	}

	void fill_xml(ParameterGroup& pg, const std::string filename) {
	    libxml2::read_xml(pg, filename);
	}

	namespace libxml2 {
	    std::string getProperty(const std::string& property,
                                    const xmlNodePtr& node_ptr)
            {
		const char* prop_value_ptr = (const char*) xmlGetProp(node_ptr, (const xmlChar*)property.c_str());
		std::string property_value(prop_value_ptr ? prop_value_ptr : "");
		xmlFree((xmlChar*) prop_value_ptr);
		return property_value;
	    }

	    void read_xml(ParameterGroup& pg, const std::string filename)
            {
		// Macro to check that the libxml version in use is compatible
		// with the version the software has been compiled against
		LIBXML_TEST_VERSION; // defined by libxml2

		xmlDocPtr doc = xmlParseFile(filename.c_str());
		if (doc == NULL) {
		    xmlFreeDoc(doc);
		    xmlCleanupParser();
		    std::cerr << "ERROR: Failed to open XML file '" << filename << "'\n";
		    throw std::exception();
		}
		xmlNodePtr root = xmlDocGetRootElement(doc);
		std::string xml_file_path = boost::filesystem::path(filename).branch_path().string();
		fill_tree(pg, root, xml_file_path);
		xmlFreeDoc(doc);
		xmlCleanupParser();
	    }

	    void fill_tree(ParameterGroup& pg, xmlNodePtr root,
                           const std::string& xml_file_path)
            {
		//std::cout << "GROUP '" << value << "' BEGIN\n";
		for (xmlNodePtr cur = root->children; cur; cur = cur->next) {
		    if (cur->type == XML_ELEMENT_NODE) {
			const char* tag_ptr = (const char*) cur->name;
			std::string tag_name(tag_ptr ? tag_ptr : "");
			if (tag_name == ID_xmltag__file_params) {
			    std::string relative_filename = getProperty(ID_xmlatt__value, cur);
			    std::string filename = (boost::filesystem::path(xml_file_path) / relative_filename).string();
			    fill_xml(pg, filename);
			    continue;
			}
			std::string name = getProperty(ID_xmlatt__name, cur);
			std::tr1::shared_ptr<ParameterMapItem> data;
			if (tag_name == ID_xmltag__param) {
			    std::string value = getProperty(ID_xmlatt__value, cur);
			    std::string type = getProperty(ID_xmlatt__type, cur);
			    if (type == ID_param_type__file) {
				value = (boost::filesystem::path(xml_file_path) / value).string();
				type  = ID_param_type__string;
			    }
			    data.reset(new Parameter(value, type));
			} else if (tag_name == ID_xmltag__param_grp) {
			    std::string child_path  = pg.path() + ID_delimiter_path + name;
			    data.reset(new ParameterGroup(child_path, &pg));
			    fill_tree(dynamic_cast<ParameterGroup&>(*data), cur, xml_file_path);
			} else if (tag_name == ID_xmltag__file_param_grp) {
			    std::string child_path  = pg.path() + ID_delimiter_path + name;
			    data.reset(new ParameterGroup(child_path, &pg));
			    std::string relative_filename = getProperty(ID_xmlatt__value, cur);
			    std::string filename = (boost::filesystem::path(xml_file_path) / relative_filename).string();
			    fill_xml(dynamic_cast<ParameterGroup&>(*data), filename);
			} else {
			    std::cerr << "ERROR: '" << tag_name << "' is an unknown xml tag.\n";
			    throw std::exception();
			}
			pg.insert(name, data);
		    }
		}
		//std::cout << "GROUP '" << id << "' END\n";
	    }
	}
    } // namespace parameter
} // namespace Dune
