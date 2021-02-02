/**
  * \file test_umat.cpp
  *
  * Tests for umat
  */

#include<umat.h>

#define BOOST_TEST_MODULE test_umat
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( testColumnToRowMajor ){
    /*!
     * Test column to row major conversion function.
     * Uses c++ arrays to avoid adding fortran code to project.
     */

    //Fake a Fortran column major array in memory with a c++ row major vector
    std::vector < double > column_major = { 1, 4,
                                            2, 5,
                                            3, 6 };
    double *column_major_pointer = column_major.data();
    const int width = 3;
    const int height = 2;
    std::vector< std::vector< double > > row_major;
    row_major = columnToRowMajor( column_major_pointer, width, height);
    std::vector< std::vector< double > > answer = { { 1, 2, 3 },
                                                    { 4, 5, 6 } };

    BOOST_CHECK( vectorTools::fuzzyEquals( row_major, answer ) );
}
