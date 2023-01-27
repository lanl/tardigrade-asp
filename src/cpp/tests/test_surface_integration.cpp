/**
  * \file test_surface_integration.cpp
  *
  * Tests for the surface integration support module
  */

#include<surface_integration.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_surface_integration
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef surfaceIntegration::floatType floatType; //!< Redefinition for the float type
typedef surfaceIntegration::floatVector floatVector; //1< Redefinition for the float vector
typedef surfaceIntegration::floatMatrix floatMatrix; //1< Redefinition for the float matrix

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( test_buildSurfacePoints ){

    floatType x0 =  0.3929383711957233;
    floatType y0 = -0.42772133009924107;
    floatType z0 = -0.5462970928715938;
    floatType dx =  0.10262953816578246;
    floatType dy =  0.43893793957112615;
    unsigned int n_points_edge_x = 5;
    unsigned int n_points_edge_y = 6;

    floatVector result;

    floatVector answer = { 0.39293837, -0.42772133, -0.54629709,  0.49556791, -0.42772133,
                          -0.54629709,  0.59819745, -0.42772133, -0.54629709,  0.70082699,
                          -0.42772133, -0.54629709,  0.80345652, -0.42772133, -0.54629709,
                           0.39293837,  0.01121661, -0.54629709,  0.49556791,  0.01121661,
                          -0.54629709,  0.59819745,  0.01121661, -0.54629709,  0.70082699,
                           0.01121661, -0.54629709,  0.80345652,  0.01121661, -0.54629709,
                           0.39293837,  0.45015455, -0.54629709,  0.49556791,  0.45015455,
                          -0.54629709,  0.59819745,  0.45015455, -0.54629709,  0.70082699,
                           0.45015455, -0.54629709,  0.80345652,  0.45015455, -0.54629709,
                           0.39293837,  0.88909249, -0.54629709,  0.49556791,  0.88909249,
                          -0.54629709,  0.59819745,  0.88909249, -0.54629709,  0.70082699,
                           0.88909249, -0.54629709,  0.80345652,  0.88909249, -0.54629709,
                           0.39293837,  1.32803043, -0.54629709,  0.49556791,  1.32803043,
                          -0.54629709,  0.59819745,  1.32803043, -0.54629709,  0.70082699,
                           1.32803043, -0.54629709,  0.80345652,  1.32803043, -0.54629709,
                           0.39293837,  1.76696837, -0.54629709,  0.49556791,  1.76696837,
                          -0.54629709,  0.59819745,  1.76696837, -0.54629709,  0.70082699,
                           1.76696837, -0.54629709,  0.80345652,  1.76696837, -0.54629709 };

    BOOST_CHECK( !surfaceIntegration::buildSurfacePoints( x0, y0, z0, dx, dy, n_points_edge_x, n_points_edge_y, result ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}

