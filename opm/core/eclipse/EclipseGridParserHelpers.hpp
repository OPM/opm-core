//===========================================================================
//
// File: EclipseGridParserHelpers.hpp
//
// Created: Tue Dec 22 11:35:32 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_ECLIPSEGRIDPARSERHELPERS_HEADER
#define OPENRS_ECLIPSEGRIDPARSERHELPERS_HEADER

#include <limits>
#include <string>
#include <istream>
#include <vector>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linInt.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

namespace Opm
{

namespace
{

    inline std::istream& ignoreLine(std::istream& is)
    {
	is.ignore(std::numeric_limits<int>::max(), '\n');
	return is;
    }

    inline std::istream& ignoreSlashLine(std::istream& is)
    {
	is.ignore(std::numeric_limits<int>::max(), '/');
	is.ignore(std::numeric_limits<int>::max(), '\n');
	return is;
    }

    inline std::istream& ignoreWhitespace(std::istream& is)
    {
	// Getting the character type facet for is()
	// We use the classic (i.e. C) locale.
	const std::ctype<char>& ct = std::use_facet< std::ctype<char> >(std::locale::classic());
	char c;
	while (is.get(c)) {
	    if (!ct.is(std::ctype_base::space, c)) {
		is.putback(c);
		break;
	    }
	}
	return is;
    }

    inline std::string readString(std::istream& is)
    {
	std::string string_candidate;
	is >> string_candidate;
        const char quote('\'');
        int beg = string_candidate[0] == quote ? 1 : 0;
        int len = string_candidate[0] == quote ? string_candidate.size() - 2 : string_candidate.size();
        return string_candidate.substr(beg, len);
    }

    // Reads data until '/' or an error is encountered.
    template<typename T>
    inline void readVectorData(std::istream& is, std::vector<T>& data)
    {
	data.clear();
	while (is) {
	    T candidate;
	    is >> candidate;
	    if (is.rdstate() & std::ios::failbit) {
		is.clear(is.rdstate() & ~std::ios::failbit);
                is >> ignoreWhitespace;
		char dummy;
		is >> dummy;
		if (dummy == '/') {
                    is >> ignoreLine;
		    break;
		} else if (dummy == '-') {  // "comment test"
		    is >> ignoreLine; // This line is a comment
		} else {
                    char buffer[1000];
                    is.getline(buffer, sizeof(buffer));
                    std::cout << buffer<<std::endl;
                    THROW("Encountered format error while reading data values. Value = " << dummy);
		}
	    } else {
		if (is.peek() == int('*')) {
		    is.ignore(); // ignore the '*'
		    int multiplier = int(candidate);
		    is >> candidate;
		    data.insert(data.end(), multiplier, candidate);
		} else {
		    data.push_back(candidate);
		}
	    }
	}
	if (!is) {
	    THROW("Encountered error while reading data values.");
	}
    }


    // Reads data items of type T. Not more than 'max_values' items.
    // Asterisks may be used to signify 'repeat counts'. 5*3.14 will
    // insert 3.14 five times. Asterisk followed by a space is used to
    // signify default values.  n* will default n consecutive quantities.
    template<class Vec>
    inline int readDefaultedVectorData(std::istream& is, Vec& data, int max_values)
    {
        ASSERT(int(data.size()) >= max_values);
	const std::ctype<char>& ct = std::use_facet< std::ctype<char> >(std::locale::classic());
        int num_values = 0;
        while (is) {
            typename Vec::value_type candidate;
            is >> candidate;
            if (is.rdstate() & std::ios::failbit) {
                is.clear(is.rdstate() & ~std::ios::failbit);
                std::string dummy;
                is >> dummy;
                if (dummy == "/") {
                    is >> ignoreLine;	
                    break;
		} else if (dummy[0] == '-') {  // "comment test"
                    is >> ignoreLine;   // This line is a comment
                } else {
                    THROW("Encountered format error while reading data values. Value = " << dummy);
                }
            } else {
                if (is.peek() == int('*')) {
                    is.ignore(); // ignore the '*'
                    int multiplier = (int)candidate;
                    if (ct.is(std::ctype_base::space, is.peek())) {
                        num_values += multiplier;  // Use default value(s)
                    } else {
                        is >> candidate;         // Use candidate 'multipler' times
                        for (int i=0; i<multiplier; ++i, ++num_values) {
                            data[num_values] = candidate;
                        }
                    }
                } else {
                    data[num_values] = candidate;
                    ++num_values;
                }
            }
            if (num_values >= max_values) {
		//is >> ignoreLine;
                break;
            }
        }
        if (!is) {
            THROW("Encountered error while reading data values.");
        }
        return num_values;
    }


    // Keywords SGOF and SWOF. Reads data until '/' or an error is encountered.
    // Default values represented by 1* is replaced by -1. Use linear interpolation
    // outside this function to replace -1.
    template<typename T>
    inline void readRelPermTable(std::istream& is, std::vector<T>& data)
    {
	data.clear();
	while (is) {
	    T candidate;
	    is >> candidate;
	    if (is.rdstate() & std::ios::failbit) {
		is.clear(is.rdstate() & ~std::ios::failbit);
		std::string dummy;
		is >> dummy;
		if (dummy == "/") {
                    is >> ignoreLine;
		    break;
		} else if (dummy[0] == '-') {  // "comment test"
		    is >> ignoreLine; // This line is a comment
		} else {
                    THROW("Encountered format error while reading data values. Value = " << dummy);
		}
	    } else {
		if (is.peek() == int('*')) {
		    is.ignore(); // ignore the '*'
		    ASSERT(int(candidate) == 1);
		    data.push_back(-1); // Set new flag for interpolation.
		} else {
		    data.push_back(candidate);
		}
	    }
	}
	if (!is) {
	    THROW("Encountered error while reading data values.");
	}
    }


    // Returns month number 1-12. Returns 0 if illegal month name. 
    inline int getMonthNumber(const std::string& month_name)
    {
        const int num_months = 12;
        std::string months[num_months] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                                          "JLY", "AUG", "SEP", "OCT", "NOV", "DEC"};
        if (month_name == "JUL") {
            return 7;    // "JUL" is an acceptable alternative to 'JLY'
        }
        int m = 0;
        for (int i=0; i<num_months; ++i) {
            if (month_name == months[i]) {
                m = i+1;
                break;
            }
        }
        return m;
    }


    inline boost::gregorian::date readDate(std::istream& is)
    {
        while (is.peek() == int('-')) {
            is >> ignoreLine;   // This line is a comment
        }
        int day, year;
        std::string month_name;
        is >> day;
        month_name = readString(is);
        is >> year;
        ignoreSlashLine(is);
        int month = getMonthNumber(month_name);

        return boost::gregorian::date(boost::gregorian::greg_year(year),
                                      boost::gregorian::greg_month(month),
                                      boost::gregorian::greg_day(day));
    }

    // Get first character after whitespace and comments, and decide
    // next action.
    // NB! Will fail for negative number in column one.
    inline int next_action(std::istream& is)
    {
	int task(0);  // 0:continue  1:return  2:throw
	const std::ctype<char>& ct =
	    std::use_facet< std::ctype<char> >(std::locale::classic());

	while (!is.eof()) {
	    is >> ignoreWhitespace;
	    char c;
	    is.get(c);
	    if (is.eof()) {
		return task;
	    }
	    is.putback(c);
	    if (ct.is(std::ctype_base::digit, c) || c== '.') {
		task = 0;   // Decimal digit. Read more data.
		break;
	    }
	    if (ct.is(std::ctype_base::alpha, c)) {
		task = 1;   // Alphabetic char. Read next keyword.
		break;
	    }
	    if (c == '-') {
		is >> ignoreLine;  // This line is a comment
	    } else {
		task = 2;  // Read error. Unexpected character.
		break;
	    }
	}
	return task;
    }

    // Reads keywords PVTG, PVTO
    typedef std::vector<std::vector<std::vector<double> > > table_t;
    inline void readPvtTable(std::istream& is, table_t& pvt_table,
			     const std::string& field_name)
    {
	const std::ctype<char>& ct =
	    std::use_facet< std::ctype<char> >(std::locale::classic());
	std::vector<double> record;
	std::vector<std::vector<double> > table;	
	while (!is.eof()) {
	    record.clear();
	    readVectorData(is, record);
	    table.push_back(record);
	    while (!is.eof()) {
		is >> ignoreWhitespace;
		char c;
		is.get(c);
                if (is.eof()) {
                    // Reached end of file, we should have pushed
                    // the last table, and emptied it. If not,
                    // we have an error.
                    if (!table.empty()) {
                        THROW("Reached EOF while still building PVT table. Missing end-of-table (slash)?");
                    }
                    return;
                }
		is.putback(c);
		if (ct.is(std::ctype_base::digit, c) || c== '.') {
		    break;   // Decimal digit. Read more records.
		}
		if (ct.is(std::ctype_base::alpha, c)) {
		    return;   // Alphabetic char. Read next keyword.
		}
		if (c == '-') {
		    is >> ignoreLine;  // This line is a comment
		    continue;
		}
		if (c == '/') {
		    is >> ignoreLine;
		    pvt_table.push_back(table);
		    table.clear();
		} else {
		    std::ostringstream oss;
		    oss << "Error reading " << field_name
		       << ". Next character is " <<  (char)is.peek();
		    THROW(oss.str());
		}
	    }
	}
    }

    // Reads keywords PVDG, PVDO, ROCKTAB
    inline void readPvdTable(std::istream& is, table_t& pvd_table,
			     const std::string& field_name, int ncol)
    {
	std::vector<double> record;
	std::vector<std::vector<double> > table(ncol);	
	while (!is.eof()) {
	    record.clear();
	    readVectorData(is, record);
	    const int rec_size = record.size()/ncol;
	    for (int k=0; k<ncol; ++k) {
		table[k].resize(rec_size);
	    }
	    for (int i=0, n=-1; i<rec_size; ++i) {
		for (int k=0; k<ncol; ++k) {
		    table[k][i] = record[++n];
		}
	    }
	    pvd_table.push_back(table);

	    int action = next_action(is); // 0:continue  1:return  2:throw
	    if (action == 1) {
		return;     // Alphabetic char. Read next keyword.
	    } else if (action == 2) {
		std::ostringstream oss;
		oss << "Error reading " << field_name
		    << ". Next character is " <<  (char)is.peek();
		THROW(oss.str());
	    }
	}
    }

    // Replace default values -1 by linear interpolation
    inline void insertDefaultValues(std::vector<std::vector<double> >& table, int ncol)
    {
	const int sz = table[0].size();
	for (int k=1; k<ncol; ++k) {
	    std::vector<int> indx;
	    std::vector<double> x;
	    for (int i=0; i<sz; ++i) {
		if (table[k][i] == -1) {
		    indx.push_back(i);
		    x.push_back(table[0][i]);    
		}
	    }
	    if (!indx.empty()) {
		std::vector<double> xv, yv;
		for (int i=0; i<sz; ++i) {
		    if (table[k][i] != -1) {
			xv.push_back(table[0][i]);		    
			yv.push_back(table[k][i]);		    
		    }
		}
		// Interpolate
		for (int i=0; i<int(indx.size()); ++i) {
		    table[k][indx[i]] = linearInterpolationExtrap(xv, yv, x[i]);
		}
	    }
	}
    }

    // Reads keywords SGOF and SWOF
    inline void readSGWOF(std::istream& is, table_t& relperm_table,
			  const std::string& field_name, int ncol)
    {
	std::vector<double> record;
	std::vector<std::vector<double> > table(ncol);	
	while (!is.eof()) {
	    record.clear();
	    readRelPermTable(is, record);
	    const int rec_size = record.size()/ncol;
	    for (int k=0; k<ncol; ++k) {
		table[k].resize(rec_size);
	    }
	    for (int i=0, n=-1; i<rec_size; ++i) {
		for (int k=0; k<ncol; ++k) {
		    table[k][i] = record[++n];
		}
	    }
	    insertDefaultValues(table, ncol);
	    relperm_table.push_back(table);

	    int action = next_action(is); // 0:continue  1:return  2:throw
	    if (action == 1) {
		return;     // Alphabetic char. Read next keyword.
	    } else if (action == 2) {
		std::ostringstream oss;
		oss << "Error reading " << field_name
		    << ". Next character is " <<  (char)is.peek();
		THROW(oss.str());
	    }
	}
    }

} // anon namespace

} // namespace Opm


#endif // OPENRS_ECLIPSEGRIDPARSERHELPERS_HEADER
