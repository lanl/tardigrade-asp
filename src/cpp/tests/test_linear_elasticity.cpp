/**
  * \file test_linear_elasticity.cpp
  *
  * Tests for the linear elasticity support module
  */

#include<linear_elasticity.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_asp
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef linearElasticity::floatType floatType; //!< Redefinition for the float type
typedef linearElasticity::floatVector floatVector; //1< Redefinition for the float vector

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

BOOST_AUTO_TEST_CASE( test_evaluateEnergy ){

     floatVector chi = { 0.39293837, -0.42772133, -0.54629709,
                         0.10262954,  0.43893794, -0.15378708,
                         0.9615284 ,  0.36965948, -0.0381362 };

     floatVector parameters = { 12.3, 43.4 };

     floatType energy_answer = 43.98356158963631;

     floatType energy;

     BOOST_CHECK( !linearElasticity::evaluateEnergy( chi, parameters, energy ) );

     BOOST_CHECK( vectorTools::fuzzyEquals( energy, energy_answer ) );

     floatVector cauchyStress_answer = { -293.41192005,   43.43903202,   60.23179078,
                                           43.43903204, -224.0643919 ,   11.23726669,
                                           60.23179078,   11.23726668, -135.38556235 };

     floatVector cauchyStress;

     energy = 0;

     BOOST_CHECK( !linearElasticity::evaluateEnergy( chi, parameters, energy, cauchyStress ) );

     BOOST_CHECK( vectorTools::fuzzyEquals( energy, energy_answer ) );

     BOOST_CHECK( vectorTools::fuzzyEquals( cauchyStress, cauchyStress_answer ) );

}
