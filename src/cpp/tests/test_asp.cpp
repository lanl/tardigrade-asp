/**
  * \file test_asp.cpp
  *
  * Tests for asp
  */

#include<asp.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_asp
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef asp::floatType floatType; //!< Redefinition of the float type
typedef asp::floatVector floatVector; //!< Redefinition of a vector of floats
typedef asp::floatMatrix floatMatrix; //!< Redefinition of a matrix of floats

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

BOOST_AUTO_TEST_CASE( testSayHello ){
    /*!
     * Test message printed to stdout in sayHello function
     */

    //Setup redirect variables for stdout
    std::stringbuf buffer;
    cout_redirect rd(&buffer);
    boost::test_tools::output_test_stream result;

    //Initialize test variables
    std::string message;
    std::string answer;
    errorOut error = NULL;

    cout_redirect guard( result.rdbuf() );

    //Check normal operation
    message = "World!";
    answer = "Hello World!\n";
    error = asp::sayHello(message);
    BOOST_CHECK( ! error );
    BOOST_CHECK( result.is_equal( answer ) );

    //Reset error code between tests
    error = NULL;

    //Check for "George" error
    message = "George";
    error = asp::sayHello(message);
    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testAbaqusInterface ){
    /*!
     * Test the asp abaqus interface
     */

    double double_scalar = 0.0;
    int int_scalar = 0;

    //Create nominally correct variable holders that match expected Abaqus Fortran interface
    //TODO: fill out nominally correct variable shape and values
    //Strings
    char CMNAME[ ] = "asp";
    //Scalar integers
    int NDI = 3;
    int NSHR = 3;
    int NTENS = 6;
    int NSTATV = 2;
    int NPROPS = 2;
    int NOEL = int_scalar;
    int NPT = int_scalar;
    int LAYER = int_scalar;
    int KSPT = int_scalar;
    int KINC = int_scalar;
    //Scalar doubles
    double SSE = double_scalar;
    double SPD = double_scalar;
    double SCD = double_scalar;
    double RPL = double_scalar;
    double DRPLDT = double_scalar;
    double DTIME = double_scalar;
    double TEMP = double_scalar;
    double DTEMP = double_scalar;
    double PNEWDT = double_scalar;
    double CELENT = double_scalar;
    //Fortan int column major arrays
    std::vector< int > jstep( 4 );
    int* JSTEP  = jstep.data( );
    //Fortan double column major arrays
    std::vector< double > stress( NTENS );
    double* STRESS = stress.data( );
    std::vector< double > statev( NSTATV );
    double* STATEV = statev.data( );
    std::vector< double > ddsdde( NTENS * NTENS );
    double* DDSDDE = ddsdde.data( );
    std::vector< double > ddsddt( NTENS );
    double* DDSDDT = ddsddt.data( );
    std::vector< double > drplde( NTENS );
    double* DRPLDE = drplde.data( );
    std::vector< double > strain( NTENS );
    double* STRAN  = strain.data( );
    std::vector< double > dstrain( NTENS );
    double* DSTRAN = dstrain.data( );
    std::vector< double > time( 2 );
    double* TIME   = time.data( );
    std::vector< double > predef( 1 );
    double* PREDEF = predef.data( );
    std::vector< double > dpred( 1 );
    double* DPRED  = dpred.data( );
    std::vector< double > props( NPROPS );
    double* PROPS  = props.data( );
    //TODO: figure out how to use asp::spatialDimensions here
    std::vector< double > coords( 3 );
    double* COORDS = coords.data( );
    //TODO: figure out how to use asp::spatialDimensions here
    std::vector< double > drot( 3 * 3);
    double* DROT   = drot.data( );
    //TODO: figure out how to use asp::spatialDimensions here
    std::vector< double > dfgrd0( 3 * 3);
    double* DFGRD0 = dfgrd0.data( );
    //TODO: figure out how to use asp::spatialDimensions here
    std::vector< double > dfgrd1( 3 * 3);
    double* DFGRD1 = dfgrd1.data( );

    //Sign of life test. Just run to see if any exceptions are thrown.
    asp::abaqusInterface(
        STRESS, STATEV, DDSDDE, SSE,    SPD,
        SCD,    RPL,    DDSDDT, DRPLDE, DRPLDT,
        STRAN,  DSTRAN, TIME,   DTIME,  TEMP,
        DTEMP,  PREDEF, DPRED,  CMNAME, NDI,
        NSHR,   NTENS,  NSTATV, PROPS,  NPROPS,
        COORDS, DROT,   PNEWDT, CELENT, DFGRD0,
        DFGRD1, NOEL,   NPT,    LAYER,  KSPT,
        JSTEP,  KINC );

    //Check for nStateVariables thrown exception
    std::vector< double > temp = { 1 };
    double* STATEV_incorrect = temp.data( );
    int NSTATV_incorrect = temp.size( );
    BOOST_CHECK_THROW(
        asp::abaqusInterface(
           STRESS, STATEV_incorrect, DDSDDE,           SSE,    SPD,
           SCD,    RPL,              DDSDDT,           DRPLDE, DRPLDT,
           STRAN,  DSTRAN,           TIME,             DTIME,  TEMP,
           DTEMP,  PREDEF,           DPRED,            CMNAME, NDI,
           NSHR,   NTENS,            NSTATV_incorrect, PROPS,  NPROPS,
           COORDS, DROT,             PNEWDT,           CELENT, DFGRD0,
           DFGRD1, NOEL,             NPT,              LAYER,  KSPT,
           JSTEP,  KINC ),
        std::exception );

    //Check for nMaterialParameters thrown exception
    temp = { 1 };
    double* PROPS_incorrect = temp.data( );
    int NPROPS_incorrect = temp.size( );

    BOOST_CHECK_THROW(
        asp::abaqusInterface(
           STRESS, STATEV, DDSDDE, SSE,             SPD,
           SCD,    RPL,    DDSDDT, DRPLDE,          DRPLDT,
           STRAN,  DSTRAN, TIME,   DTIME,           TEMP,
           DTEMP,  PREDEF, DPRED,  CMNAME,          NDI,
           NSHR,   NTENS,  NSTATV, PROPS_incorrect, NPROPS_incorrect,
           COORDS, DROT,   PNEWDT, CELENT,          DFGRD0,
           DFGRD1, NOEL,   NPT,    LAYER,           KSPT,
           JSTEP,  KINC ),
        std::exception );

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeLocalParticleEnergyDensity ){
    /*
     * Test the default implementation of the computation of the local particle's energy density
     */

    floatType previousTime = 1.2;

    floatType deltaTime = 0.4;

    floatVector currentDeformationGradient = { 1, 2, 3,
                                               4, 5, 6,
                                               7, 8, 9 };

    floatVector previousDeformationGradient = { .1, .2, .3,
                                                .4, .5, .6,
                                                .7, .8, .9 };

    floatType currentTemperature = 295.4;

    floatType previousTemperature = 0.43;

    floatVector previousStateVariables;

    floatVector parameters = { 30., 20. };

    floatType energy = 0;

    floatVector cauchyStress;

    asp::aspBase asp;

    BOOST_CHECK( !asp.computeLocalParticleEnergyDensity( previousTime, deltaTime, currentDeformationGradient, previousDeformationGradient,
                                                         currentTemperature, previousTemperature, previousStateVariables, parameters, energy, cauchyStress ) );

    BOOST_CHECK( !vectorTools::fuzzyEquals( energy, 0. ) );

    BOOST_CHECK( cauchyStress.size( ) == currentDeformationGradient.size( ) );

    energy = 0;

    cauchyStress.clear( );

    floatType logProbabilityRatio = 3.4;

    BOOST_CHECK( !asp.computeLocalParticleEnergyDensity( previousTime, deltaTime, currentDeformationGradient, previousDeformationGradient,
                                                         currentTemperature, previousTemperature, previousStateVariables, parameters, energy, cauchyStress, logProbabilityRatio ) );

    BOOST_CHECK( !vectorTools::fuzzyEquals( energy, 0. ) );

    BOOST_CHECK( cauchyStress.size( ) == currentDeformationGradient.size( ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( logProbabilityRatio, 0. ) );

}
