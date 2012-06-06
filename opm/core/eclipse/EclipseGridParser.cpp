//===========================================================================
//
// File: EclipseGridParser.C
//
// Created: Thu Dec  6 08:46:05 2007
//
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
// Revision: $Id: EclipseGridParser.C,v 1.4 2008/08/18 14:16:14 atgeirr Exp $
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
#include <iostream>
#include <fstream>
#include <exception>
#include <algorithm>
#include <limits>
#include <numeric>
#include <cfloat>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/eclipse/EclipseGridParserHelpers.hpp>
#include <opm/core/eclipse/SpecialEclipseFields.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <boost/filesystem.hpp>
#include <opm/core/utility/Units.hpp>

using namespace std;

//#define VERBOSE

namespace Opm
{

// ---------- List of supported keywords ----------

namespace EclipseKeywords
{
    string integer_fields[] =
        { string("ACTNUM"),
          string("SATNUM"),
          string("EQLNUM"),
          string("REGNUM"),
          string("ROCKTYPE"),
          string("DIMENS"),
          string("REGDIMS"),
          string("WELLDIMS"),
          string("TABDIMS"),
          string("FIPNUM"),
          string("GRIDFILE")
        };
    const int num_integer_fields = sizeof(integer_fields) / sizeof(integer_fields[0]);

    string floating_fields[] =
        { string("COORD"),    string("ZCORN"),      string("PERMX"),
          string("PERMY"),    string("PERMZ"),      string("PERMXX"),
          string("PERMYY"),   string("PERMZZ"),     string("PERMXY"),
          string("PERMYZ"),   string("PERMZX"),     string("PORO"),
          string("BULKMOD"),  string("YOUNGMOD"),   string("LAMEMOD"),
          string("SHEARMOD"), string("POISSONMOD"), string("PWAVEMOD"),
          string("MULTPV"),   string("PRESSURE"),   string("SGAS"),
          string("SWAT"),     string("SOIL"),       string("RS"),
          string("DXV"),      string("DYV"),        string("DZV"),
          string("DEPTHZ")
        };
    const int num_floating_fields = sizeof(floating_fields) / sizeof(floating_fields[0]);

    string special_fields[] =
        { string("SPECGRID"), string("FAULTS"), string("MULTFLT"),
          string("TITLE"),    string("START"),  string("DATES"),
          string("DENSITY"),  string("PVDG"),   string("PVDO"),
          string("PVTG"),     string("PVTO"),   string("PVTW"),
          string("SGOF"),     string("SWOF"),   string("ROCK"),
          string("ROCKTAB"),  string("WELSPECS"), string("COMPDAT"),
          string("WCONINJE"), string("WCONPROD"), string("WELTARG"),
          string("WELOPEN"),
          string("EQUIL"),    string("PVCDO"),    string("TSTEP"),
          string("PLYVISC"),  string("PLYROCK"),  string("PLYADS"),
          string("PLYMAX"),   string("TLMIXPAR"), string("WPOLYMER"),
          string("GRUPTREE"), string("GCONINJE"), string("GCONPROD"),
          string("WGRUPCON"),
          // The following fields only have a dummy implementation
          // that allows us to ignore them.
          string("SWFN"),
          string("SOF2"),
          string("TUNING")
        };

    const int num_special_fields = sizeof(special_fields) / sizeof(special_fields[0]);

    string ignore_with_data[] =
        { string("MAPUNITS"), string("MAPAXES"),  string("GRIDUNIT"),
          string("NTG"),      string("REGDIMS"),  string("WELLDIMS"),
          string("NSTACK"),   string("SATNUM"),
          string("RPTRST"),   string("ROIP"),     string("RWIP"),
          string("RWSAT"),    string("RPR"),      string("WBHP"),
          string("WOIR"),     string("BOX"),
          string("COORDSYS"), string("PBVD")
        };
    const int num_ignore_with_data = sizeof(ignore_with_data) / sizeof(ignore_with_data[0]);

    string ignore_no_data[] =
        { string("RUNSPEC"), string("WATER"),    string("OIL"),
          string("METRIC"),  string("FMTIN"),    string("FMTOUT"),
          string("GRID"),    string("INIT"),     string("NOECHO"),
          string("ECHO"),    string("EDIT"),     string("PROPS"),
          string("REGIONS"), string("SOLUTION"), string("SUMMARY"),
          string("FPR"),     string("FOIP"),     string("FWIP"),
          string("RUNSUM"),  string("EXCEL"),    string("SCHEDULE"),
          string("END"),     string("ENDBOX"),   string("CONTINUE"),
          string("NONNC"),   string("GAS"),      string("DISGAS"),
          string("FIELD"),   string("POLYMER")
        };
    const int num_ignore_no_data = sizeof(ignore_no_data) / sizeof(ignore_no_data[0]);

    string include_keywords[] = { string("INCLUDE") };
    const int num_include_keywords = sizeof(include_keywords) / sizeof(include_keywords[0]);


} // namespace EclipseKeywords

namespace {

    enum FieldType {
        Integer,
        FloatingPoint,
        Timestepping,
        SpecialField,
        IgnoreWithData,
        IgnoreNoData,
        Include,
        Unknown
    };

    inline FieldType classifyKeyword(const string& keyword)
    {
        using namespace EclipseKeywords;
        if (count(integer_fields, integer_fields + num_integer_fields, keyword)) {
            return Integer;
        } else if (count(floating_fields, floating_fields + num_floating_fields, keyword)) {
            return FloatingPoint;
        } else if (keyword == "TSTEP" || keyword == "DATES") {
            return Timestepping;
        } else if (count(special_fields, special_fields + num_special_fields, keyword)) {
            return SpecialField;
        } else if (count(ignore_with_data, ignore_with_data + num_ignore_with_data, keyword)) {
            return IgnoreWithData;
        } else if (count(ignore_no_data, ignore_no_data + num_ignore_no_data, keyword)) {
            return IgnoreNoData;
        } else if (count(include_keywords, include_keywords + num_include_keywords, keyword)) {
            return Include;
        } else {
            return Unknown;
        }
    }

    inline std::string upcase(const std::string& s)
    {
        std::string us(s);
        // Getting the character type facet for toupper().
        // We use the classic (i.e. C) locale.
        const std::ctype<char>& ct = std::use_facet< std::ctype<char> >(std::locale::classic());
        for (int i = 0; i < int(s.size()); ++i) {
            us[i] = ct.toupper(s[i]);
        }
        return us;
    }

    inline std::string readKeyword(std::istream& is)
    {
        std::string keyword_candidate;
        while (!is.eof()) {
            is >> keyword_candidate;
            if(keyword_candidate.find("--") == 0) {
                is >> ignoreLine;  // This line is a comment
            } else {
                return upcase(keyword_candidate);
            }
        }
        return "CONTINUE";  // Last line in included file is a comment
    }

    inline bool readKeywordNew(std::istream& is, std::string& keyword)
    {
        char buf[9];
        int i, j;
        char c;
        /* Clear buf */
        for (i=0; i<9; ++i) {
            buf[i] = '\0';
        }

        /* Read first character and check if it is uppercase*/
        //buf[0] = fgetc(fp);
        is.get(buf[0]);
        if ( !isupper( buf[0] ) ) {
            is.unget();
            return false;          /* NOT VALID CHARACTER */
        }

        /* Scan as much as possible possible keyword, 8 characters long */
        i = 1;
        is.get(c);
        while ( (is.good()) &&
                (c != EOF     ) &&
                (!isblank(c)   ) &&
                (isupper(c) || isdigit(c)) &&
                (c != '\n'    ) &&
                (c != '/'     ) &&
                (i < 8        )) {
            buf[i++] = c;
            is.get(c);
        }

        /* Skip rest of line */
        if (c != '\n'){
            is.get(c);
            while ( (is.good()) &&
                    (c != EOF     ) &&
                    (c != '\n'    )) {
                is.get(c);
            }
        }
        if(c == '\n') {
            is.unget();
        }

        /* Find first non-uppercase or non-digit character */
        for (i=0; i<8; ++i) {
            if ( !(isupper(buf[i]) || isdigit(buf[i])) ) {
                break;
            }
        }

        /* Check if remaining characters are blank */
        for (j = i; j<8; ++j) {
            if(!isspace(buf[j]) && buf[j] != '\0') {
                return false; /* CHARACTER AFTER SPACE OR INVALID CHARACTER */
            }
            buf[j] = '\0';
        }
        keyword = std::string(buf);
        std::string::size_type end = keyword.find_last_of('\0');
        if(end != keyword.npos)
        keyword = keyword.substr(0, end+1);
        return true;
    }

} // anon namespace



// ---------- Member functions ----------

/// Default constructor.
//---------------------------------------------------------------------------
EclipseGridParser::EclipseGridParser()
//---------------------------------------------------------------------------
    : current_reading_mode_(Regular),
      start_date_(boost::date_time::not_a_date_time),
      current_time_days_(0.0),
      current_epoch_(0)
{
}


/// Constructor taking an eclipse filename.
//---------------------------------------------------------------------------
EclipseGridParser::EclipseGridParser(const string& filename, bool convert_to_SI)
//---------------------------------------------------------------------------
    : current_reading_mode_(Regular),
      start_date_(boost::date_time::not_a_date_time),
      current_time_days_(0.0),
      current_epoch_(0)
{
    // Store directory of filename
    boost::filesystem::path p(filename);
    directory_ = p.parent_path().string();
    ifstream is(filename.c_str());
    if (!is) {
        cerr << "Unable to open file " << filename << endl;
        throw exception();
    }
    read(is, convert_to_SI);
}


/// Read the given stream, overwriting any previous data.
//---------------------------------------------------------------------------
void EclipseGridParser::read(istream& is, bool convert_to_SI)
//---------------------------------------------------------------------------
{
    integer_field_map_.clear();
    floating_field_map_.clear();
    special_field_by_epoch_.clear();
    special_field_by_epoch_.push_back(SpecialMap());

    readImpl(is);

    current_epoch_ = 0;

    computeUnits();
    if (convert_to_SI) {
        convertToSI();
    }

#define VERBOSE_LIST_FIELDS 1
#if VERBOSE_LIST_FIELDS
    std::cout << "\nInteger fields:\n";
    for (std::map<string, std::vector<int> >::iterator
             i = integer_field_map_.begin(); i != integer_field_map_.end(); ++i)
        std::cout << '\t' << i->first << '\n';

    std::cout << "\nFloat fields:\n";
    for (std::map<string, std::vector<double> >::iterator
             i = floating_field_map_.begin(); i != floating_field_map_.end(); ++i)
        std::cout << '\t' << i->first << '\n';

    std::cout << "\nSpecial fields:\n";
    for (int epoch = 0; epoch < numberOfEpochs(); ++epoch) {
        std::cout << "Epoch " << epoch << '\n';
        const SpecialMap& sm = special_field_by_epoch_[epoch];
        for (SpecialMap::const_iterator i = sm.begin(); i != sm.end(); ++i) {
            std::cout << '\t' << i->first << '\n';
        }
    }
#endif
}

//---------------------------------------------------------------------------
void EclipseGridParser::readImpl(istream& is)
//---------------------------------------------------------------------------
{
    if (!is) {
        cerr << "Could not read given input stream." << endl;
        throw exception();
    }

    // Make temporary maps that will at the end be swapped with the
    // member maps
    // NOTE: Above is no longer true, for easier implementation of
    //       the INCLUDE keyword. We lose the strong exception guarantee,
    //       though (of course retaining the basic guarantee).
    map<string, vector<int> >& intmap = integer_field_map_;
    map<string, vector<double> >& floatmap = floating_field_map_;

    // Actually read the data
    std::string keyword;
    while (is.good()) {
        is >> ignoreWhitespace;
        bool ok = readKeywordNew(is, keyword);
        if (ok) {
            //#ifdef VERBOSE
            cout << "Keyword found: " << keyword << endl;
            //#endif
            FieldType type = classifyKeyword(keyword);
            // std::cout << "Classification: " << type << std::endl;
            switch (type) {
            case Integer: {
                readVectorData(is, intmap[keyword]);
                break;
            }
            case FloatingPoint: {
                readVectorData(is, floatmap[keyword]);
                break;
            }
            case Timestepping: {
                SpecialMap& sm = special_field_by_epoch_[current_epoch_];
                if (start_date_.is_not_a_date()) {
                    // Set it to START date, or default if no START.
                    // This will only ever happen in the first epoch,
                    // upon first encountering a timestepping keyword.
                    SpecialMap::const_iterator it = sm.find("START");
                    if (hasField("START")) {
                        start_date_ = getSTART().date;
                    } else {
                        start_date_ = boost::gregorian::date(1983, 1, 1);
                    }
                }
                if (current_reading_mode_ == Regular) {
                    current_reading_mode_ = Timesteps;
                }
                // Get current epoch's TSTEP, if it exists, create new if not.
                SpecialMap::iterator it = sm.find("TSTEP");
                TSTEP* tstep = 0;
                if (it != sm.end()) {
                    tstep = dynamic_cast<TSTEP*>(it->second.get());
                } else {
                    SpecialFieldPtr sb_ptr(new TSTEP());
                    tstep = dynamic_cast<TSTEP*>(sb_ptr.get());
                    sm["TSTEP"] = sb_ptr;
                }
                ASSERT(tstep != 0);
                // Append new steps to current TSTEP object
                if (keyword == "TSTEP") {
                    const int num_steps_old = tstep->tstep_.size();
                    tstep->read(is); // This will append to the TSTEP object.
                    const double added_days
                        = std::accumulate(tstep->tstep_.begin() + num_steps_old, tstep->tstep_.end(), 0.0);
                    current_time_days_ += added_days;
                } else if (keyword == "DATES") {
                    DATES dates;
                    dates.read(is);
                    for (std::size_t dix = 0; dix < dates.dates.size(); ++dix) {
                        boost::gregorian::date_duration since_start = dates.dates[dix] - start_date_;
                        double step = double(since_start.days()) - current_time_days_;
                        tstep->tstep_.push_back(step);
                        current_time_days_ = double(since_start.days());
                    }
                } else {
                    THROW("Keyword " << keyword << " cannot be handled here.");
                }
                break;
            }
            case SpecialField: {
                if (current_reading_mode_ == Timesteps) {
                    // We have been reading timesteps, but have
                    // now encountered something else.
                    // That means we are in a new epoch.
                    current_reading_mode_ = Regular;
                    special_field_by_epoch_.push_back(SpecialMap());
                    ++current_epoch_;
                    ASSERT(int(special_field_by_epoch_.size()) == current_epoch_ + 1);
                }
                SpecialFieldPtr sb_ptr = createSpecialField(is, keyword);
                if (sb_ptr) {
                    special_field_by_epoch_[current_epoch_][keyword] = sb_ptr;
                } else {
                    THROW("Could not create field " << keyword);
                }
                break;
            }
            case IgnoreWithData: {
                ignored_fields_.insert(keyword);
                //is >> ignoreSlashLine;
                //#ifdef VERBOSE
                // cout << "(ignored)" << endl;
                //#endif
                break;
            }
            case IgnoreNoData: {
                ignored_fields_.insert(keyword);
                //is >> ignoreLine;
                //#ifdef VERBOSE
                // cout << "(ignored)" << endl;
                //#endif
                break;
            }
            case Include: {
                string include_filename = readString(is);
                if (!directory_.empty()) {
                    include_filename = directory_ + '/' + include_filename;
                }
                ifstream include_is(include_filename.c_str());
                if (!include_is) {
                    THROW("Unable to open INCLUDEd file " << include_filename);
                }
                readImpl(include_is);
                //              is >> ignoreSlashLine;
                break;
            }
            case Unknown:
            default:
                ignored_fields_.insert(keyword);
                cout << "*** Warning: keyword " << keyword << " is unknown." << endl;
                //is >> ignoreSlashLine;
                //throw exception();
            }
        } else {
            // if (!ok)
            is >> ignoreLine;
        }
    }
}



//---------------------------------------------------------------------------
void EclipseGridParser::convertToSI()
//---------------------------------------------------------------------------
{
    // Convert all special fields.
    typedef SpecialMap::iterator SpecialIt;
    for (int epoch = 0; epoch < numberOfEpochs(); ++epoch) {
        SpecialMap& sm = special_field_by_epoch_[epoch];
        for (SpecialIt i = sm.begin(); i != sm.end(); ++i) {
            i->second->convertToSI(units_);
        }
    }

    // Convert all floating point fields.
    typedef std::map<string, std::vector<double> >::iterator FloatIt;
    for (FloatIt i = floating_field_map_.begin(); i != floating_field_map_.end(); ++i) {
        const std::string& key = i->first;
        std::vector<double>& field = i->second;
        // Find the right unit.
        double unit = 1e100;
        bool do_convert = true;
        if (key == "COORD" || key == "ZCORN" ||
            key == "DXV"   || key == "DYV"   || key == "DZV" ||
            key == "DEPTHZ") {
            unit = units_.length;
        } else if (key == "PERMX"  || key == "PERMY"  || key == "PERMZ"  ||
                   key == "PERMXX" || key == "PERMYY" || key == "PERMZZ" ||
                   key == "PERMXY" || key == "PERMYZ" || key == "PERMZX") {
            unit = units_.permeability;
        } else if (key == "PORO"     || key == "BULKMOD"  || key == "YOUNGMOD" ||
                   key == "LAMEMOD"  || key == "SHEARMOD" || key == "POISSONMOD" ||
                   key == "PWAVEMOD" || key == "MULTPV"   || key == "PWAVEMOD" ||
                   key == "SGAS"     || key == "SWAT"     || key == "SOIL"     ||
                   key == "RS") {
            unit = 1.0;
            do_convert = false; // Dimensionless keywords...
        } else if (key == "PRESSURE") {
            unit = units_.pressure;
        } else {
            THROW("Units for field " << key << " not specified. Cannon convert to SI.");
        }

        if (do_convert) {
            for (std::vector<double>::size_type j = 0; j < field.size(); ++j) {
                field[j] = unit::convert::from(field[j], unit);
            }
        }
    }

    // Set all units to one.
    units_.setToOne();
}


/// Returns true if the given keyword corresponds to a field that
/// was found in the file.
//---------------------------------------------------------------------------
bool EclipseGridParser::hasField(const string& keyword) const
//---------------------------------------------------------------------------
{
    string ukey = upcase(keyword);
    return integer_field_map_.count(ukey) || floating_field_map_.count(ukey) ||
        special_field_by_epoch_[current_epoch_].count(ukey) || ignored_fields_.count(ukey);
}


/// Returns true if all the given keywords correspond to fields
/// that were found in the file.
//---------------------------------------------------------------------------
bool EclipseGridParser::hasFields(const vector<string>& keywords) const
//---------------------------------------------------------------------------
{
    int num_keywords = keywords.size();
    for (int i = 0; i < num_keywords; ++i) {
        if (!hasField(keywords[i])) {
            return false;
        }
    }
    return true;
}

//---------------------------------------------------------------------------
vector<string> EclipseGridParser::fieldNames() const
//---------------------------------------------------------------------------
{
    vector<string> names;
    names.reserve(integer_field_map_.size() +
                  floating_field_map_.size() +
                  special_field_by_epoch_[current_epoch_].size() +
                  ignored_fields_.size());
    {
        map<string, vector<int> >::const_iterator it = integer_field_map_.begin();
        for (; it != integer_field_map_.end(); ++it) {
            names.push_back(it->first);
        }
    }
    {
        map<string, vector<double> >::const_iterator it = floating_field_map_.begin();
        for (; it != floating_field_map_.end(); ++it) {
            names.push_back(it->first);
        }
    }
    {
        SpecialMap::const_iterator it = special_field_by_epoch_[current_epoch_].begin();
        for (; it != special_field_by_epoch_[current_epoch_].end(); ++it) {
            names.push_back(it->first);
        }
    }
    {
        set<string>::const_iterator it = ignored_fields_.begin();
        for (; it != ignored_fields_.end(); ++it) {
            names.push_back(*it);
        }
    }
    return names;
}


//---------------------------------------------------------------------------
int EclipseGridParser::numberOfEpochs() const
//---------------------------------------------------------------------------
{
    return special_field_by_epoch_.size();
}


//---------------------------------------------------------------------------
void EclipseGridParser::setCurrentEpoch(int epoch)
//---------------------------------------------------------------------------
{
    ASSERT(epoch >= 0 && epoch < numberOfEpochs());
    current_epoch_ = epoch;
}


//---------------------------------------------------------------------------
const std::vector<int>& EclipseGridParser::getIntegerValue(const std::string& keyword) const
//---------------------------------------------------------------------------
{
    if (keyword == "SPECGRID") {
        cerr << "\nERROR. Interface has changed!\n"
             << "const vector<int>& dim = parser.getIntegerValue(""SPECGRID"") is deprecated.\n"
             << "Use:\n"
             << "const SPECGRID& specgrid = parser.getSPECGRID();\n"
             << "const vector<int>& dim = specgrid.dimensions;\n\n";
        throw exception();
    }

    map<string, vector<int> >::const_iterator it
        = integer_field_map_.find(keyword);
    if (it == integer_field_map_.end()) {
        THROW("No such field: " << keyword);
    } else {
        return it->second;
    }
}

//---------------------------------------------------------------------------
const std::vector<double>& EclipseGridParser::getFloatingPointValue(const std::string& keyword) const
//---------------------------------------------------------------------------
{
    map<string, vector<double> >::const_iterator it
        = floating_field_map_.find(keyword);
    if (it == floating_field_map_.end()) {
        THROW("No such field: " << keyword);
    } else {
        return it->second;
    }
}


//---------------------------------------------------------------------------
const std::tr1::shared_ptr<SpecialBase> EclipseGridParser::getSpecialValue(const std::string& keyword) const
//---------------------------------------------------------------------------
{
    SpecialMap::const_iterator it = special_field_by_epoch_[current_epoch_].find(keyword);
    if (it == special_field_by_epoch_[current_epoch_].end()) {
        THROW("No such field: " << keyword);
    } else {
        return it->second;
    }
}

//---------------------------------------------------------------------------
std::tr1::shared_ptr<SpecialBase>
EclipseGridParser::createSpecialField(std::istream& is,
                                      const std::string& fieldname)
//---------------------------------------------------------------------------
{
    string ukey = upcase(fieldname);
    std::tr1::shared_ptr<SpecialBase> spec_ptr
        = Factory<SpecialBase>::createObject(fieldname);
    is >> ignoreWhitespace;
    spec_ptr->read(is);
    return spec_ptr;
}

//---------------------------------------------------------------------------
void EclipseGridParser::setIntegerField(const std::string& keyword,
                                        const std::vector<int>& field)
//---------------------------------------------------------------------------
{
    integer_field_map_[keyword] = field;
}

//---------------------------------------------------------------------------
void EclipseGridParser::setFloatingPointField(const std::string& keyword,
                                              const std::vector<double>& field)
//---------------------------------------------------------------------------
{
    floating_field_map_[keyword] = field;
}

//---------------------------------------------------------------------------
void EclipseGridParser::setSpecialField(const std::string& keyword,
                                        std::tr1::shared_ptr<SpecialBase> field)
//---------------------------------------------------------------------------
{
    special_field_by_epoch_[current_epoch_][keyword] = field;
}

//---------------------------------------------------------------------------
const EclipseUnits& EclipseGridParser::units() const
//---------------------------------------------------------------------------
{
    return units_;
}

//---------------------------------------------------------------------------
void EclipseGridParser::computeUnits()
//---------------------------------------------------------------------------
{
    // Decide unit family.
    enum EclipseUnitFamily { Metric = 0, Field = 1, Lab = 2, Pvtm = 3 };
    EclipseUnitFamily unit_family = Metric; // The default.
    if (hasField("FIELD")) unit_family = Field;
    if (hasField("LAB")) unit_family = Lab;
    if (hasField("PVT-M")) unit_family = Pvtm;

    // Set units.
    using namespace prefix;
    using namespace unit;
    switch (unit_family) {
    case Metric:
        units_.length = meter;
        units_.time = day;
        units_.density = kilogram/cubic(meter);
        units_.polymer_density = kilogram/cubic(meter);
        units_.pressure = barsa;
        units_.compressibility = 1.0/barsa;
        units_.viscosity = centi*Poise;
        units_.permeability = milli*darcy;
        units_.liqvol_s = cubic(meter);
        units_.liqvol_r = cubic(meter);
        units_.gasvol_s = cubic(meter);
        units_.gasvol_r = cubic(meter);
        units_.transmissibility = centi*Poise * cubic(meter) / (day * barsa);
        break;
    case Field:
        units_.length = feet;
        units_.time = day;
        units_.density = pound/cubic(feet);
        units_.polymer_density = pound/stb;
        units_.pressure = psia;
        units_.compressibility = 1.0/psia;
        units_.viscosity = centi*Poise;
        units_.permeability = milli*darcy;
        units_.liqvol_s = stb;
        units_.liqvol_r = stb;
        units_.gasvol_s = 1000*cubic(feet);  // Prefix 'M' is 1000
        units_.gasvol_r = stb;
        units_.transmissibility = centi*Poise * stb / (day * psia);
        break;
    case Lab:
        THROW("Unhandled unit family " << unit_family);
        break;
    case Pvtm:
        THROW("Unhandled unit family " << unit_family);
        break;
    default:
        THROW("Unknown unit family " << unit_family);
    }
}

} // namespace Opm
