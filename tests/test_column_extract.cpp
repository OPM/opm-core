


#define BOOST_TEST_DYN_LINK
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE ColumnExtractTest
#include <boost/test/unit_test.hpp>
#include <opm/core/utility/ColumnExtract.hpp>
#include <opm/core/GridManager.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>

BOOST_AUTO_TEST_CASE(SingleColumnTest)
{
    using namespace Opm;

    const int size_x = 1, size_y = 1, size_z = 10;
    GridManager manager(size_x, size_y, size_z);

    std::vector<std::vector<int> > columns;
    extractColumn(*manager.c_grid(), columns);

    std::vector<int> correct_answer;
    for (int i = 0; i < 10; ++i) {
        correct_answer.push_back(i);
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(correct_answer.begin(), correct_answer.end(),
                                  columns[0].begin(), columns[0].end());

}


BOOST_AUTO_TEST_CASE(FourByFourColumnTest)
{
    const int size_x = 4, size_y = 4, size_z = 10;
    using namespace Opm;
    GridManager manager(size_x, size_y, size_z);

    std::vector<std::vector<int> > columns;
    extractColumn(*manager.c_grid(), columns);

    std::vector<std::vector<int> > correct_answer;
    correct_answer.resize(size_x * size_y);
    for(int i = 0; i < size_x * size_y; i++) {
        for(int j = 0; j < 10; j++) {
            correct_answer[i].push_back( i + j*size_x*size_y);
        }
    }

    for(int i = 0; i < size_x * size_y; i++) {
        BOOST_CHECK_EQUAL_COLLECTIONS(correct_answer[i].begin(), correct_answer[i].end(),
                                      columns[i].begin(), columns[i].end());
    }

}



BOOST_AUTO_TEST_CASE(DisjointColumn)
{
    std::string grdecl =
"SPECGRID                                        \n"
"3 3 3 1 F                                       \n"
"/                                               \n"
"                                                \n"
"COORD                                           \n"
"0.0 0 0       0.0 0 3.0                         \n"
"1.0 0 0       1.0 0 3.0                         \n"
"2.0 0 0       2.0 0 3.0                         \n"
"3.0 0 0       3.0 0 3.0                         \n"
"                                                \n"
"0.0 1.0 0       0.0 1.0 3.0                     \n"
"1.0 1.0 0       1.0 1.0 3.0                     \n"
"2.0 1.0 0       2.0 1.0 3.0                     \n"
"3.0 1.0 0       3.0 1.0 3.0                     \n"
"                                                \n"
"0.0 2.0 0       0.0 2.0 3.0                     \n"
"1.0 2.0 0       1.0 2.0 3.0                     \n"
"2.0 2.0 0       2.0 2.0 3.0                     \n"
"3.0 2.0 0       3.0 2.0 3.0                     \n"
"                                                \n"
"0.0 3.0 0       0.0 3.0 3.0                     \n"
"1.0 3.0 0       1.0 3.0 3.0                     \n"
"2.0 3.0 0       2.0 3.0 3.0                     \n"
"3.0 3.0 0       3.0 3.0 3.0                     \n"
"/                                               \n"
"                                                \n"
"ZCORN                                           \n"
"36*0.0                                          \n"
"72*1.0                                          \n"
"72*2.0                                          \n"
"36*3.0                                          \n"
"/                                               \n"
"                                                \n"
"ACTNUM                                          \n"
"13*1                                            \n"
"0                                               \n"
"13*1                                            \n"
"/                                               \n"
"\n";
    using namespace Opm;
    std::istringstream is(grdecl);
    EclipseGridParser deck;
    deck.read(is);
    GridManager manager(deck);

    std::vector<std::vector<int> > columns;
    extractColumn(*manager.c_grid(), columns);
    // for (size_t i = 0; i < columns.size(); ++i) {
    //     for (size_t j = 0; j < columns[i].size(); ++j) {
    //         std::cout << columns[i][j] << ' ';
    //     }
    //     std::cout << '\n';
    // }
    // std::cout << std::endl;

    const int correct[10][3] = { { 0, 1, 2 },
                                 { 3, 4, 5 },
                                 { 6, 7, 8 },
                                 { 9, 10, 11 },
                                 { 13, -1, -1 },
                                 { 14, 15, 16 },
                                 { 17, 18, 19 },
                                 { 20, 21, 22 },
                                 { 23, 24, 25 },
                                 { 12, -1, -1 } };
    std::vector<std::vector<int> > correct_answer;
    correct_answer.resize(10);
    for(int i = 0; i < 10; i++) {
        for(int j = 0; j < 3; j++) {
            correct_answer[i].push_back(correct[i][j]);
        }
    }
    correct_answer[4].resize(1);
    correct_answer[9].resize(1);

    BOOST_CHECK_EQUAL(columns.size(), 10);

    for(int i = 0; i < 10; i++) {
        BOOST_CHECK_EQUAL_COLLECTIONS(correct_answer[i].begin(), correct_answer[i].end(),
                                      columns[i].begin(), columns[i].end());
    }

}
