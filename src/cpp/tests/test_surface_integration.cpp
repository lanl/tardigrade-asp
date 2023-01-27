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

BOOST_AUTO_TEST_CASE( test_rotatePoints ){

    floatVector points = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794,
                          -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                          -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                           0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                           0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                          -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

    floatType thetaX = -0.8157901201098496;
    floatType thetaY = -0.13259765464094353;
    floatType thetaZ = -0.1382744733407124;

    floatVector answer = { 0.29876665, -0.73920054, -0.0103761 ,  0.1824409 ,  0.16523358,
                          -0.40775359,  1.01375449,  0.08664871, -0.16562965, -0.2664906 ,
                           0.15690994,  0.50906204, -0.28995028, -0.71889812,  0.48095937,
                           0.33986109, -0.9639332 ,  0.080416  ,  0.07689573,  0.23098794,
                           0.14498284,  0.77369311,  0.36619883, -0.08082919,  0.36661806,
                          -0.49925749,  0.12665336, -0.60908828, -0.0080854 ,  0.40391106 };

    floatVector result;

    BOOST_CHECK( !surfaceIntegration::rotatePoints( points, thetaX, thetaY, thetaZ, result ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}
