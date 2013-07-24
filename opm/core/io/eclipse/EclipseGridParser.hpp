//===========================================================================
//                                                                           
// File: EclipseGridParser.h                                                 
//                                                                           
// Created: Wed Dec  5 17:05:13 2007                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//                                                                           
// Revision: $Id: EclipseGridParser.h,v 1.3 2008/08/18 14:16:13 atgeirr Exp $
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

#ifndef OPM_ECLIPSEGRIDPARSER_HEADER
#define OPM_ECLIPSEGRIDPARSER_HEADER

#include <iosfwd>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <boost/shared_ptr.hpp>
#include <opm/core/io/eclipse/SpecialEclipseFields.hpp>
#include <opm/core/io/eclipse/EclipseUnits.hpp>
#include <opm/core/utility/Factory.hpp>

#include <opm/core/grid/cornerpoint_grid.h>
#ifdef HAVE_ERT
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_grid.h>
#endif


namespace Opm
{

/**
   @brief A class for reading and parsing all fields of an eclipse file.
   
   This object is constructed using an Eclipse .grdecl-file. All data
   fields are extracted upon construction and written to vector data
   structures, which can then be read out afterwards via
   convenience functions.

   There is also a convenience function to easily check which fields
   were successfully parsed.

   @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
   @date 2007/12/06 13:00:03

*/

  enum FieldType {
    Integer,
    FloatingPoint,
    Timestepping,
    SpecialField,
    IgnoreWithData,
    IgnoreNoData,
    Include,
    Import,
    Unknown
  };



  class EclipseGridParser
  {
  public:
    /// Default constructor.
    EclipseGridParser();
    /// Constructor taking an eclipse filename. Unless the second
    /// argument 'convert_to_SI' is false, all fields will be
    /// converted to SI units.
    explicit EclipseGridParser(const std::string& filename, bool convert_to_SI = true);

    static FieldType classifyKeyword(const std::string& keyword);
    static bool readKeyword(std::istream& is, std::string& keyword);


    /// Read the given stream, overwriting any previous data.  Unless
    /// the second argument 'convert_to_SI' is false, all fields will
    /// be converted to SI units.
    void read(std::istream& is, bool convert_to_SI = true);

    /// Convert all data to SI units, according to unit category as
    /// specified in the eclipse file (METRIC, FIELD etc.).
    /// It is unnecessary to call this if constructed with
    /// 'convert_to_SI' equal to true, but it is not an error.
    void convertToSI();

    /// Returns true if the given keyword corresponds to a field that
    /// was found in the file.
    bool hasField(const std::string& keyword) const;
    /// Returns true if all the given keywords correspond to fields
    /// that were found in the file.
    bool hasFields(const std::vector<std::string>& keywords) const;
    /// The keywords/fields found in the file.
    std::vector<std::string> fieldNames() const;

    /// Returns the number of distinct epochs found in the deck SCHEDULE.
    int numberOfEpochs() const;

    /// Sets the current epoch.
    /// Valid arguments are in [0, ..., numberOfEpochs() - 1].
    /// After reading, current epoch always starts at 0.
    void setCurrentEpoch(int epoch);

    /// Returns the start_date_
    boost::gregorian::date getStartDate() const;

    /// Returns a reference to a vector containing the values
    /// corresponding to the given integer keyword.
    const std::vector<int>& getIntegerValue(const std::string& keyword) const;

    /// Returns a reference to a vector containing the values
    /// corresponding to the given floating-point keyword.
    const std::vector<double>& getFloatingPointValue(const std::string& keyword) const;

    typedef boost::shared_ptr<SpecialBase> SpecialFieldPtr;

    /// Returns a reference to a vector containing pointers to the values 
    /// corresponding to the given keyword when the values are not only integers
    /// or floats.
    const SpecialFieldPtr getSpecialValue(const std::string& keyword) const;

    // This macro implements support for a special field keyword. It requires that a subclass
    // of SpecialBase exists, that has the same name as the keyword.
    // After using SPECIAL_FIELD(KEYWORD), the public member getKEYWORD will be available.
#define SPECIAL_FIELD(keyword)                                                                   \
private:                                                                                         \
    struct X##keyword { X##keyword() { Factory<SpecialBase>::addCreator<keyword>(#keyword); } }; \
    X##keyword x##keyword;                                                                       \
public:                                                                                          \
    const keyword& get##keyword() const                                                          \
    { return dynamic_cast<const keyword&>(*getSpecialValue(#keyword)); }

    // Support for special fields.
    SPECIAL_FIELD(SPECGRID)
    SPECIAL_FIELD(FAULTS)
    SPECIAL_FIELD(MULTFLT)
    SPECIAL_FIELD(TITLE)
    SPECIAL_FIELD(START)
    SPECIAL_FIELD(DATES)
    SPECIAL_FIELD(DENSITY)
    SPECIAL_FIELD(PVDG)
    SPECIAL_FIELD(PVDO)
    SPECIAL_FIELD(PVTG)
    SPECIAL_FIELD(PVTO)
    SPECIAL_FIELD(PVTW)
    SPECIAL_FIELD(SGOF)
    SPECIAL_FIELD(SWOF)
    SPECIAL_FIELD(ROCK)
    SPECIAL_FIELD(ROCKTAB)
    SPECIAL_FIELD(WELSPECS)
    SPECIAL_FIELD(COMPDAT)
    SPECIAL_FIELD(WCONINJE)
    SPECIAL_FIELD(WCONPROD)
    SPECIAL_FIELD(WELTARG)
    SPECIAL_FIELD(WELOPEN)
    SPECIAL_FIELD(EQUIL)
    SPECIAL_FIELD(PVCDO)
    SPECIAL_FIELD(TSTEP)
    SPECIAL_FIELD(PLYVISC)
    SPECIAL_FIELD(PLYROCK)
    SPECIAL_FIELD(PLYADS)
    SPECIAL_FIELD(PLYMAX)
    SPECIAL_FIELD(TLMIXPAR)
    SPECIAL_FIELD(WPOLYMER)
    SPECIAL_FIELD(GRUPTREE)
    SPECIAL_FIELD(GCONINJE)
    SPECIAL_FIELD(GCONPROD)
    SPECIAL_FIELD(WGRUPCON)
    SPECIAL_FIELD(ENDSCALE)
    SPECIAL_FIELD(SCALECRS)
    SPECIAL_FIELD(ENPTVD)
    SPECIAL_FIELD(ENKRVD)

    // The following fields only have a dummy implementation
    // that allows us to ignore them.
    SPECIAL_FIELD(SWFN)
    SPECIAL_FIELD(SOF2)
    SPECIAL_FIELD(TUNING)

#undef SPECIAL_FIELD


    /// Sets an integer field to have a particular value.
    void setIntegerField(const std::string& keyword, const std::vector<int>& field);

    /// Sets a floating point field to have a particular value.
    void setFloatingPointField(const std::string& keyword, const std::vector<double>& field);

    /// Sets a special field to have a particular value.
    void setSpecialField(const std::string& keyword, SpecialFieldPtr field);

    /// Compute the units used by the deck, depending on the presence
    /// of keywords such as METRIC, FIELD etc.  It is an error to call
    /// this after conversion to SI has taken place.
    void computeUnits();

    /// The units specified by the eclipse file read.
    const EclipseUnits& units() const;

    struct grdecl get_grdecl() const;

    /// Save grid parts of deck in EGRID format. 
    void saveEGRID(const std::string & filename, std::vector<int>& actnum) const;

#ifdef HAVE_ERT
  void saveEGRID_INIT( const std::string& output_dir , const std::string& basename, bool fmt_file = false);
  void saveINIT( const std::string & filename , const ecl_grid_type * ecl_grid);
  ecl_grid_type * newGrid( );
#endif


private:

#ifdef HAVE_ERT
  ecl_kw_type * newEclKW(const std::string &keyword , ecl_type_enum ecl_type) const;
  void          save_kw( fortio_type * fortio , const std::string & kw , ecl_type_enum ecl_type);
#endif

    SpecialFieldPtr createSpecialField(std::istream& is, const std::string& fieldname);
    SpecialFieldPtr cloneSpecialField(const std::string& fieldname,
                                      const boost::shared_ptr<SpecialBase> original);
    void readImpl(std::istream& is);
    void getNumericErtFields(const std::string& filename);


    std::string directory_;
    std::map<std::string, std::vector<int> > integer_field_map_;
    std::map<std::string, std::vector<double> > floating_field_map_;
    // std::map<std::string, SpecialFieldPtr> special_field_map_;
    std::set<std::string> ignored_fields_;
    EclipseUnits units_;

    // For SCHEDULE handling.
    enum ReadingMode { Regular, Timesteps };
    ReadingMode current_reading_mode_;
    boost::gregorian::date start_date_;
    double current_time_days_;
    int current_epoch_;
    typedef std::map<std::string, SpecialFieldPtr> SpecialMap;
    std::vector<SpecialMap> special_field_by_epoch_;
};



} // namespace Opm

#endif // OPM_ECLIPSEGRIDPARSER_HEADER
