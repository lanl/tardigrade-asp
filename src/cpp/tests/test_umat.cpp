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

    //Note that a c++ vector is stored as row major and the column to row major function accesses by pointer.
    //This means that the c++ test is really testing the pointer access of a row-major vector and conversion to a two
    //dimensional array. To test a column major memory access, we must add a fortran executable to the project.
    std::vector < double > column_major = { 1, 2,
                                            3, 4, 
                                            5, 6 };
    double *column_major_pointer = &column_major[0];
    const int width = 3;
    const int height = 2;
    std::vector< std::vector< double > > row_major;
    row_major = columnToRowMajor( column_major_pointer, width, height); 
    std::vector< std::vector< double > > answer = { { 1, 2, 3 },
                                                    { 4, 5, 6 } };

    BOOST_CHECK( row_major == answer );
}
