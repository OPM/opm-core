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

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE

#define BOOST_TEST_MODULE CompressedPropertyTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <opm/core/utility/CompressedPropertyAccess.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid.h>

struct SetupSimple {
    SetupSimple()
    {
        Opm::ParserPtr parser(new Opm::Parser());
        deck = parser->parseFile("compressed_gridproperty.data");
        ecl.reset(new Opm::EclipseState(deck));
    }

    Opm::DeckConstPtr         deck;
    Opm::EclipseStateConstPtr ecl;
};


template <class Setup>
struct TestFixture : public Setup
{
    TestFixture()
        : Setup ()
        , grid  (deck)
        , reltol(1.0e-10)
    {
    }

    using Setup::deck;
    using Setup::ecl;

    Opm::GridManager grid;
    const double     reltol;
};


BOOST_AUTO_TEST_SUITE(CompressedPropertyHandling)

// Construct global array by extracting an "undefined" (unspecified)
// double vector from the input deck.
BOOST_FIXTURE_TEST_CASE(APExtractDoubleUndefined,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<double> ECLGlobalDoubleArray;

    ECLGlobalDoubleArray multpv(ecl, "MULTPV", 1.0);

    BOOST_REQUIRE_CLOSE(multpv[0], 1.0, reltol);
    BOOST_REQUIRE_CLOSE(multpv[1], 1.0, reltol);
    BOOST_REQUIRE_CLOSE(multpv[2], 1.0, reltol);
    BOOST_REQUIRE_CLOSE(multpv[3], 1.0, reltol);
}


// Construct global array by extracting a fully defined (specified)
// double vector from the input deck.
BOOST_FIXTURE_TEST_CASE(APExtractDoubleDefined,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<double> ECLGlobalDoubleArray;

    ECLGlobalDoubleArray ntg(ecl, "NTG", 1.0);

    BOOST_REQUIRE_CLOSE(ntg[0], 0.1, reltol);
    BOOST_REQUIRE_CLOSE(ntg[1], 0.2, reltol);
    BOOST_REQUIRE_CLOSE(ntg[2], 0.3, reltol);
    BOOST_REQUIRE_CLOSE(ntg[3], 0.4, reltol);
}


// Construct global array by extracting an undefined (unspecified)
// integer (int) vector from the input deck.
BOOST_FIXTURE_TEST_CASE(APExtractIntUndefined,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<int> ECLGlobalIntArray;

    ECLGlobalIntArray imbnum(ecl, "IMBNUM", 1);

    BOOST_REQUIRE_EQUAL(imbnum[0], 1);
    BOOST_REQUIRE_EQUAL(imbnum[1], 1);
    BOOST_REQUIRE_EQUAL(imbnum[2], 1);
    BOOST_REQUIRE_EQUAL(imbnum[3], 1);
}


// Construct global array by extracting a fully defined (specified)
// integer (int) vector from the input deck.
BOOST_FIXTURE_TEST_CASE(APExtractIntDefined,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<int> ECLGlobalIntArray;

    ECLGlobalIntArray satnum(ecl, "SATNUM", 1);

    BOOST_REQUIRE_EQUAL(satnum[0], 4);
    BOOST_REQUIRE_EQUAL(satnum[1], 3);
    BOOST_REQUIRE_EQUAL(satnum[2], 2);
    BOOST_REQUIRE_EQUAL(satnum[3], 1);
}


// Construct global, infinitely sized, double array by specifying a
// single constant for all cells.
BOOST_FIXTURE_TEST_CASE(APConstantDouble,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::Constant<double> ConstantDoubleArray;

    const double c = 1.234e5;

    ConstantDoubleArray x(c);

    BOOST_REQUIRE_CLOSE(x[    0], c, reltol);
    BOOST_REQUIRE_CLOSE(x[    1], c, reltol);
    BOOST_REQUIRE_CLOSE(x[    2], c, reltol);
    BOOST_REQUIRE_CLOSE(x[10000], c, reltol); // Infinite size array
}


// Construct global, infinitely sized, integer (int) array by
// specifying a single constant for all cells.
BOOST_FIXTURE_TEST_CASE(APConstantInt,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::Constant<int> ConstantIntArray;

    const int i = 12345;

    ConstantIntArray a(i);

    BOOST_REQUIRE_EQUAL(a[    0], i);
    BOOST_REQUIRE_EQUAL(a[    1], i);
    BOOST_REQUIRE_EQUAL(a[    2], i);
    BOOST_REQUIRE_EQUAL(a[10001], i); // Inifinite size array
}


// Construct compressed double array based on global, undefined array
// extracted from input deck.  Default ("any") type-check tag.
BOOST_FIXTURE_TEST_CASE(CAExtractDoubleUndefinedAny,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<double> ECLGlobalDoubleArray;

    typedef Opm::GridPropertyAccess::
        Compressed<ECLGlobalDoubleArray> CompressedArray;

    ECLGlobalDoubleArray multpv_glob(ecl, "MULTPV", 1.0);
    CompressedArray multpv(multpv_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_CLOSE(multpv[0], 1.0, reltol);
    BOOST_REQUIRE_CLOSE(multpv[1], 1.0, reltol);
}


// Construct compressed double array based on global, fully specified
// array extracted from input deck.  Type-check tag: NTG.
BOOST_FIXTURE_TEST_CASE(CAExtractDoubleDefinedNTG,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<double> ECLGlobalDoubleArray;

    typedef Opm::GridPropertyAccess::Tag::NTG NTG;

    typedef Opm::GridPropertyAccess::
        Compressed<ECLGlobalDoubleArray, NTG> CompressedArray;

    ECLGlobalDoubleArray ntg_glob(ecl, "NTG", 1.0);
    CompressedArray ntg(ntg_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_CLOSE(ntg[0], 0.2, reltol);
    BOOST_REQUIRE_CLOSE(ntg[1], 0.4, reltol);
}


// Construct compressed integer (int) array based on global, undefined
// (unspecified) array extracted from input deck.  Default ("any")
// type-check tag.
BOOST_FIXTURE_TEST_CASE(CAExtractIntUndefinedAny,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<int> ECLGlobalIntArray;

    typedef Opm::GridPropertyAccess::
        Compressed<ECLGlobalIntArray> CompressedArray;

    ECLGlobalIntArray imbnum_glob(ecl, "IMBNUM", 1);
    CompressedArray imbnum(imbnum_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_EQUAL(imbnum[0], 1);
    BOOST_REQUIRE_EQUAL(imbnum[1], 1);
}


// Construct compressed integer (int) array based on global, fully
// specified array extracted from input deck.  Custom type-check tag.
BOOST_FIXTURE_TEST_CASE(CAExtractIntDefinedCustom,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::ExtractFromDeck<int> ECLGlobalIntArray;

    struct RegionID : public Opm::GridPropertyAccess::Tag::Any {};

    typedef Opm::GridPropertyAccess::
        Compressed<ECLGlobalIntArray, RegionID> CompressedArray;

    ECLGlobalIntArray satnum_glob(ecl, "SATNUM", 1);
    CompressedArray satnum(satnum_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_EQUAL(satnum[0], 3);
    BOOST_REQUIRE_EQUAL(satnum[1], 1);
}


// Construct compressed double array based on global constant value
// for all cells.  Default ("any") type-check tag.
BOOST_FIXTURE_TEST_CASE(CAConstantDoubleAny,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::Constant<double> ConstantDoubleArray;

    typedef Opm::GridPropertyAccess::
        Compressed<ConstantDoubleArray> CompressedArray;

    const double c = 1.234e5;

    ConstantDoubleArray x_glob(c);
    CompressedArray x(x_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_CLOSE(x[0], c, reltol);
    BOOST_REQUIRE_CLOSE(x[1], c, reltol);
}


// Construct compressed double array based on global constant value
// for all cells.  Custom type-check tag.
BOOST_FIXTURE_TEST_CASE(CAConstantDoubleCustom,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::Constant<double> ConstantDoubleArray;

    struct MyTag : public Opm::GridPropertyAccess::Tag::Any {};

    typedef Opm::GridPropertyAccess::
        Compressed<ConstantDoubleArray, MyTag> CompressedArray;

    const double c = 1.234e5;

    ConstantDoubleArray x_glob(c);
    CompressedArray x(x_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_CLOSE(x[0], c, reltol);
    BOOST_REQUIRE_CLOSE(x[1], c, reltol);
}


// Construct compressed integer (int) array based on global constant
// value for all cells.  Custom type-check tag.
BOOST_FIXTURE_TEST_CASE(CAConstantIntAny,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::Constant<int> ConstantIntArray;

    typedef Opm::GridPropertyAccess::
        Compressed<ConstantIntArray> CompressedArray;

    const int i = 12345;

    ConstantIntArray x_glob(i);
    CompressedArray x(x_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_EQUAL(x[0], i);
    BOOST_REQUIRE_EQUAL(x[1], i);
}


// Construct compressed integer (int) array based on global constant
// value for all cells.  Custom type-check tag.
BOOST_FIXTURE_TEST_CASE(CAConstantIntCustom,
                        TestFixture<SetupSimple>)
{
    typedef Opm::GridPropertyAccess::ArrayPolicy
        ::Constant<int> ConstantIntArray;

    struct MyTag : public Opm::GridPropertyAccess::Tag::Any {};

    typedef Opm::GridPropertyAccess::
        Compressed<ConstantIntArray, MyTag> CompressedArray;

    const int i = 12345;

    ConstantIntArray x_glob(i);
    CompressedArray x(x_glob, grid.c_grid()->global_cell);

    BOOST_REQUIRE_EQUAL(x[0], i);
    BOOST_REQUIRE_EQUAL(x[1], i);
}

BOOST_AUTO_TEST_SUITE_END()
