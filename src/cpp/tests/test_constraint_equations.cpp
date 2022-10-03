/**
  * \file test_constraint_equations.cpp
  *
  * Tests for the constraint equations support module
  */

#include<constraint_equations.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_constraint_equations
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef constraintEquations::floatType floatType; //!< Redefinition for the float type
typedef constraintEquations::floatVector floatVector; //1< Redefinition for the float vector
typedef constraintEquations::floatMatrix floatMatrix; //1< Redefinition for the float matrix

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

BOOST_AUTO_TEST_CASE( test_tractionConstraint ){
    /*!
     * Test of the traction constraint
     */

    floatVector cauchyStress = { 0.69646919, 0.28613933, 0.22685145,
                                 0.55131477, 0.71946897, 0.42310646,
                                 0.9807642 , 0.68482974, 0.4809319 };

    floatVector normal = { 0.39211752, 0.34317802, 0.72904971 };

    floatVector traction = { 0.43857224, 0.0596779 , 0.39804426 };

    floatType P = 0.7379954057320357;

    floatType C_answer = 0.8992828453968;

    floatType C;

    BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction, P, C ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( C, C_answer ) );

    floatVector dCdCauchyStress, dCdNormal, dCdTraction;

    floatType dCdP;

    BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction, P, C, dCdCauchyStress, dCdNormal, dCdTraction, dCdP ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( C, C_answer ) );

    floatType C_2;

    floatVector dCdCauchyStress_2, dCdNormal_2, dCdTraction_2;

    floatType dCdP_2;

    floatVector d2CdCauchyStressdNormal, d2CdCauchyStressdP, d2CdNormaldP, d2CdTractiondP;

    BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction, P, C_2, dCdCauchyStress_2, dCdNormal_2, dCdTraction_2, dCdP_2,
                                                           d2CdCauchyStressdNormal, d2CdCauchyStressdP, d2CdNormaldP, d2CdTractiondP ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( C, C_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdCauchyStress_2, dCdCauchyStress ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdNormal_2, dCdNormal ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdTraction_2, dCdTraction ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdP_2, dCdP ) );

    floatType eps = 1e-6;

    floatVector dCdCauchyStress_answer( cauchyStress.size( ), 0 );

    floatVector dCdNormal_answer( normal.size( ), 0 );

    floatVector dCdTraction_answer( traction.size( ), 0 );

    floatType   dCdP_answer = 0;

    for ( unsigned int i = 0; i < cauchyStress.size( ); i++ ){

        floatVector delta( cauchyStress.size( ), 0 );

        delta[ i ] += eps * std::abs( cauchyStress[ i ] ) + eps;

        floatType Cp, Cm;

        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress + delta, normal, traction, P, Cp ) );

        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress - delta, normal, traction, P, Cm ) );

        dCdCauchyStress_answer[ i ] = ( Cp - Cm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdCauchyStress, dCdCauchyStress_answer ) );

    floatVector d2CdCauchyStressdNormal_answer( cauchyStress.size( ) * normal.size( ), 0 );

    for ( unsigned int i = 0; i < normal.size( ); i++ ){

        floatVector delta( normal.size( ), 0 );

        delta[ i ] += eps * std::abs( normal[ i ] ) + eps;

        floatType Cp, Cm;

        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal + delta, traction, P, Cp ) );

        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal - delta, traction, P, Cm ) );

        dCdNormal_answer[ i ] = ( Cp - Cm ) / ( 2 * delta[ i ] );

        floatVector dCdCauchyStressp, dCdCauchyStressm;
    
        floatVector dCdNormalp, dCdNormalm;
    
        floatVector dCdTractionp, dCdTractionm;
    
        floatType dCdPp, dCdPm;
    
        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal + delta, traction, P, Cp, dCdCauchyStressp, dCdNormalp, dCdTractionp, dCdPp ) );
    
        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal - delta, traction, P, Cm, dCdCauchyStressm, dCdNormalm, dCdTractionm, dCdPm ) );

        for ( unsigned int j = 0; j < normal.size( ); j++ ){

            for ( unsigned int k = 0; k < normal.size( ); k++ ){

                d2CdCauchyStressdNormal_answer[ normal.size( ) * normal.size( ) * j + normal.size( ) * k + i ] += ( dCdCauchyStressp[ normal.size( ) * j + k ] - dCdCauchyStressm[ normal.size( ) * j + k ] ) / ( 2 * delta[ i ] );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdNormal, dCdNormal_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2CdCauchyStressdNormal, d2CdCauchyStressdNormal_answer ) );

    for ( unsigned int i = 0; i < traction.size( ); i++ ){

        floatVector delta( traction.size( ), 0 );

        delta[ i ] += eps * std::abs( traction[ i ] ) + eps;

        floatType Cp, Cm;

        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction + delta, P, Cp ) );

        BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction - delta, P, Cm ) );

        dCdTraction_answer[ i ] = ( Cp - Cm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdTraction, dCdTraction_answer ) );

    floatType delta = eps * std::abs( C ) + eps;

    floatType Cp, Cm;

    BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction, P + delta, Cp ) );

    BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction, P - delta, Cm ) );

    dCdP_answer = ( Cp - Cm ) / ( 2 * delta );

    floatVector dCdCauchyStressp, dCdCauchyStressm;

    floatVector dCdNormalp, dCdNormalm;

    floatVector dCdTractionp, dCdTractionm;

    floatType dCdPp, dCdPm;

    BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction, P + delta, Cp, dCdCauchyStressp, dCdNormalp, dCdTractionp, dCdPp ) );

    BOOST_CHECK( !constraintEquations::tractionConstraint( cauchyStress, normal, traction, P - delta, Cm, dCdCauchyStressm, dCdNormalm, dCdTractionm, dCdPm ) );

    floatVector d2CdCauchyStressdP_answer = ( dCdCauchyStressp - dCdCauchyStressm ) / ( 2 * delta );

    floatVector d2CdNormaldP_answer = ( dCdNormalp - dCdNormalm ) / ( 2 * delta );

    floatVector d2CdTractiondP_answer = ( dCdTractionp - dCdTractionm ) / ( 2 * delta );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCdP, dCdP_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2CdCauchyStressdP, d2CdCauchyStressdP_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2CdNormaldP, d2CdNormaldP_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2CdTractiondP, d2CdTractiondP_answer ) );

}
