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
typedef linearElasticity::floatMatrix floatMatrix; //1< Redefinition for the float matrix

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

    unsigned int spatialDimensions = 3;

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

    cauchyStress.clear( );

    energy = 0;

    floatVector dEnergydChi;

    floatMatrix dCauchyStressdChi;

    BOOST_CHECK( !linearElasticity::evaluateEnergy( chi, parameters, energy, cauchyStress, dEnergydChi, dCauchyStressdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy, energy_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cauchyStress, cauchyStress_answer ) );

    // Perform tests of the gradients

    floatType eps = 1e-6;

    floatVector dEnergydChi_answer( chi.size( ), 0 );
    floatMatrix dCauchyStressdChi_answer( chi.size( ), floatVector( chi.size( ), 0 ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] = eps * std::abs( chi[ i ] ) + eps;

        floatType ep, em;

        floatVector cauchyStressp, cauchyStressm;

        BOOST_CHECK( !linearElasticity::evaluateEnergy( chi + delta, parameters, ep ) );

        BOOST_CHECK( !linearElasticity::evaluateEnergy( chi - delta, parameters, em ) );

        dEnergydChi_answer[ i ] = ( ep - em ) / ( 2 * delta[ i ] );

        BOOST_CHECK( !linearElasticity::evaluateEnergy( chi + delta, parameters, ep, cauchyStressp ) );

        BOOST_CHECK( !linearElasticity::evaluateEnergy( chi - delta, parameters, em, cauchyStressm ) );

        for ( unsigned int j = 0; j < chi.size( ); j++ ){

            dCauchyStressdChi_answer[ j ][ i ] = ( cauchyStressp[ j ] - cauchyStressm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dEnergydChi, dEnergydChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCauchyStressdChi, dCauchyStressdChi_answer ) );

}
