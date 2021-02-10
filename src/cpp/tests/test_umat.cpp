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

BOOST_AUTO_TEST_CASE( expandAbaqusStandardStressVector ){
    /*!
     * Test expansion of stress and strain type components to full Abaqus vectors
     *
     * See the Abaqus documentation > Introduction & Spatial Modeling > Conventions chapter > Convention used for stress
     * and strain components.
     *
     * The stress vector components for Abaqus/Standard (UMAT) are
     *
     * \f$ \left { \sigma_{11}, \sigma_{22}, \sigma_{33}, \tau_{12}, \tau_{13}, \tau_{23} } \f$
     *
     * and the strain vector components match as
     *
     * \f$ \left { \epsilon_{11}, \epsilon_{22}, \epsilon_{33}, \gamma_{12}, \gamma_{13}, \gamma_{23} } \f$
     *
     * where components that are zero-valued by definition, e.g. plane stress, are omitted. The shear strain is the
     * engineering shear strain where
     *
     * \f$ \gamma_{ij} = \epsilon_{ij} + \epsilon_{ji} \f$
     */

     //Initialize common test variables
     int NDI;
     int NSHR;
     std::vector< double > vector_expansion(6, -666.);

     //Test full size vector
     std::vector< double > abaqus_full = { 11, 22, 33, 12, 13, 23 };
     std::vector< double > expected_full = { 11, 22, 33, 12, 13, 23 };
     NDI = 3;
     NSHR = 3;

     vector_expansion = expandAbaqusStandardStressVector( abaqus_full, &NDI, &NSHR );

     BOOST_CHECK( vectorTools::FuzzyEquals( vector_expansion, expected_full );

     //Test plane stress vector
     std::fill(vector_expansion.begin(), vector_expansion.end(), -666.);
     std::vector< double > abaqus_plane_stress = { 11, 22, 12 };
     std::vector< double > expected_plane_stress = { 11, 22, 0., 12, 0., 0. };
     NDI = 2;
     NSHR = 1;

     vector_expansion = expandAbaqusStandardStressVector( abaqus_full, NDI, NSHR );

     BOOST_CHECK( vectorTools::FuzzyEquals( vector_expansion, expected_full );
}
