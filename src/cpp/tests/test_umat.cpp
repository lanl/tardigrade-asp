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
//    std::vector < double > column_major = { 1, 4,
//                                            2, 5, 
//                                            3, 6 };
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

    //DEBUGGING
    int column_major_index;
    std::cout << "column_major_pointer" << std::endl;
    for (int row = 0; row < height; row++){
        for (int col = 0; col < width; col++){
            column_major_index = row*width + col;
            std::cout << column_major_pointer[column_major_index] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "row_major" << std::endl;
    for (int row = 0; row < height; row++){
        for (int col = 0; col < width; col++){
            std::cout << row_major[row][col] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    BOOST_CHECK( row_major == answer );
}
