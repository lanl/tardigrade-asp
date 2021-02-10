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
     * Uses c++ vectors and pointers to avoid adding fortran code to project.
     */

    //Fake a Fortran column major array in memory with a c++ row major vector
    std::vector< double > column_major = { 1, 4,
                                           2, 5,
                                           3, 6 };
    double *column_major_pointer = column_major.data();
    const int height = 2;
    const int width = 3;
    std::vector< std::vector< double > > row_major;
    row_major = columnToRowMajor( column_major_pointer, height, width);
    std::vector< std::vector< double > > answer = { { 1, 2, 3 },
                                                    { 4, 5, 6 } };

    BOOST_CHECK( vectorTools::fuzzyEquals( row_major, answer ) );
}

BOOST_AUTO_TEST_CASE( testRowToColumnMajor ){
    /*!
     * Test row to column major conversion function.
     * Uses c++ vectors and pointers to avoid adding fortran code to project.
     */

    //Test the interface using a c++ vector of vectors
    //Fake a Fortran column major array in memory with a c++ row major vector
    std::vector< double > column_major(6);
    double *column_major_pointer = column_major.data( );
    std::vector< double > expected = { 1, 4,
                                       2, 5,
                                       3, 6 };
    const int height = 2;
    const int width = 3;
    std::vector< std::vector< double > > row_major = { { 1, 2, 3 },
                                                       { 4, 5, 6 } };
    rowToColumnMajor( column_major_pointer, row_major, height, width );

    BOOST_CHECK( vectorTools::fuzzyEquals( column_major, expected ) );

    //Test the interface using a c++ vector saved in row major order
    //Reset fake Fortran column major array
    column_major = { 0, 0, 0, 0, 0, 0 };
    std::vector< double > row_major_vector = { 1, 2, 3,
                                               4, 5, 6 };
    rowToColumnMajor( column_major_pointer, row_major_vector, height, width );

    BOOST_CHECK( vectorTools::fuzzyEquals( column_major, expected ) );

    //Test a single row to single column
    //Create a fake Fortran vector in memory with a c++ row major vector
    std::vector< double > fortran_vector(3);
    double *fortran_vector_pointer = fortran_vector.data( );
    std::vector< double > expected_vector = { 1, 2, 3 };
    std::vector< double > cpp_vector = { 1, 2, 3 };
    rowToColumnMajor( fortran_vector_pointer, cpp_vector, 1, 3 );

    BOOST_CHECK( vectorTools::fuzzyEquals( fortran_vector, expected_vector ) );
}
