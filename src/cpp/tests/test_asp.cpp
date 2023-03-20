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

// Unit tester to open private members of aspBase for unit testing
namespace asp{

    namespace unit_test{

        class aspBaseTester{

            public:

                static errorOut initializeUnitSphere( asp::aspBase &asp,
                                                      std::pair< bool, floatVector > & unitSpherePoints,
                                                      std::pair< bool, std::vector< unsigned int > > & unitSphereConnectivity ){
    
                    BOOST_CHECK( !asp.initializeUnitSphere( ) );
    
                    unitSpherePoints = asp._unitSpherePoints;
    
                    unitSphereConnectivity = asp._unitSphereConnectivity;
    
                    return NULL;
    
                }

                static errorOut setLocalReferenceRadius( asp::aspBase &asp, std::pair< bool, floatType > &radius ){

                    BOOST_CHECK( !asp.setLocalReferenceRadius( ) );

                    radius = asp._localReferenceRadius;

                    return NULL;

                }

                static errorOut setNonLocalReferenceRadius( asp::aspBase &asp, std::pair< bool, floatType > &radius ){

                    BOOST_CHECK( !asp.setNonLocalReferenceRadius( ) );

                    radius = asp._nonlocalReferenceRadius;

                    return NULL;

                }

                static errorOut setLocalReferenceNormal( asp::aspBase &asp, std::pair< bool, floatVector > &normal ){

                    BOOST_CHECK( !asp.setLocalReferenceNormal( ) );

                    normal = asp._localReferenceNormal;

                    return NULL;

                }

                static errorOut setLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp, std::pair< bool, floatVector > &Xi ){

                    BOOST_CHECK( !asp.setLocalSurfaceReferenceRelativePositionVector( ) );

                    Xi = asp._localSurfaceReferenceRelativePositionVector;

                    return NULL;

                }

                static errorOut setNonLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp, std::pair< bool, floatVector > &Xi ){

                    BOOST_CHECK( !asp.setNonLocalSurfaceReferenceRelativePositionVector( ) );

                    Xi = asp._nonlocalSurfaceReferenceRelativePositionVector;

                    return NULL;

                }

                static errorOut setLocalDeformationGradient( asp::aspBase &asp, std::pair< bool, floatVector > &localDeformationGradient ){

                    BOOST_CHECK( !asp.setLocalDeformationGradient( ) );

                    localDeformationGradient = asp._localDeformationGradient;

                    return NULL;

                }

                static errorOut setLocalMicroDeformation( asp::aspBase &asp, std::pair< bool, floatVector > &localMicroDeformation ){

                    BOOST_CHECK( !asp.setLocalMicroDeformation( ) );

                    localMicroDeformation = asp._localMicroDeformation;

                    return NULL;

                }

                static errorOut setReferenceDistanceVector( asp::aspBase &asp, std::pair< bool, floatVector > &referenceDistanceVector ){

                    BOOST_CHECK( !asp.setReferenceDistanceVector( ) );

                    referenceDistanceVector = asp._referenceDistanceVector;

                    return NULL;

                }

                static errorOut setLocalReferenceParticleSpacing( asp::aspBase &asp,
                                                                  std::pair< bool, floatVector > &localReferenceParticleSpacing ){

                    BOOST_CHECK( !asp.setLocalReferenceParticleSpacing( ) );

                    localReferenceParticleSpacing = asp._localReferenceParticleSpacing;

                    return NULL;

                }

                static errorOut setNonLocalMicroDeformation( asp::aspBase &asp,
                                                             std::pair< bool, floatVector > &nonlocalMicroDeformation ){

                    BOOST_CHECK( !asp.setNonLocalMicroDeformation( ) );

                    nonlocalMicroDeformation = asp._nonlocalMicroDeformation;

                    return NULL;

                }

                static errorOut setLocalCurrentNormal( asp::aspBase &asp,
                                                       std::pair< bool, floatVector > &localCurrentNormal ){

                    BOOST_CHECK( !asp.setLocalCurrentNormal( ) );

                    localCurrentNormal = asp._localCurrentNormal;

                    return NULL;

                }

                static errorOut setCurrentDistance( asp::aspBase &asp,
                                                    std::pair< bool, floatVector > &currentDistance ){

                    BOOST_CHECK( !asp.setCurrentDistance( ) );

                    currentDistance = asp._currentDistanceVector;

                    return NULL;

                }

                static errorOut setSurfaceParameters( asp::aspBase &asp,
                                                      std::pair< bool, floatVector > &currentSurfaceParameters ){

                    BOOST_CHECK( !asp.setSurfaceParameters( ) );

                    currentSurfaceParameters = asp._surfaceParameters;

                    return NULL;

                }

                static errorOut initializeSurfaceIntegrandQuantities( asp::aspBase &asp ){

                    BOOST_CHECK( !asp.initializeSurfaceIntegrandQuantities( ) );

                    return NULL;

                }

                static errorOut computeSurfaceEnergyDensity( asp::aspBase &asp, floatType &surfaceEnergyDensity ){

                    errorOut error = asp.computeSurfaceEnergyDensity( surfaceEnergyDensity );

                    if ( error ){
                        error->print( );
                    }

                    return NULL;

                }

                // Direct write functions for mocking
                static errorOut set_indices( asp::aspBase &asp,
                                            unsigned int localIndex, unsigned int nonlocalIndex, unsigned int localSurfaceNodeIndex ){

                    asp._localIndex = localIndex;

                    asp._nonlocalIndex = nonlocalIndex;

                    asp._localSurfaceNodeIndex = localSurfaceNodeIndex;

                    return NULL;

                }

                static errorOut set_radius( asp::aspBase &asp,
                                            floatType &radius ){

                    asp._radius = radius;

                    return NULL;

                }

                static errorOut set_unitSphere( asp::aspBase &asp,
                                                const floatVector &points, std::vector< unsigned int > &connectivity ){

                    asp._unitSpherePoints.first = true;
                    asp._unitSpherePoints.second = points;

                    asp._unitSphereConnectivity.first = true;
                    asp._unitSphereConnectivity.second = connectivity;

                    return NULL;

                }

                static errorOut set_localReferenceNormal( asp::aspBase &asp,
                                                          const floatVector &normal ){

                    asp._localReferenceNormal.first = true;
                    asp._localReferenceNormal.second = normal;

                    return NULL;

                }

                static errorOut set_localReferenceRadius( asp::aspBase &asp,
                                                          const floatType &radius ){

                    asp._localReferenceRadius.first = true;
                    asp._localReferenceRadius.second = radius;

                    return NULL;

                }

                static errorOut set_nonlocalReferenceRadius( asp::aspBase &asp,
                                                             const floatType &radius ){

                    asp._nonlocalReferenceRadius.first = true;
                    asp._nonlocalReferenceRadius.second = radius;

                    return NULL;

                }

                static errorOut set_deformationGradient( asp::aspBase &asp,
                                                         const floatVector &deformationGradient ){

                    asp._deformationGradient = deformationGradient;

                    return NULL;

                }

                static errorOut set_microDeformation( asp::aspBase &asp,
                                                      const floatVector &microDeformation ){

                    asp._microDeformation = microDeformation;

                    return NULL;

                }

                static errorOut set_localDeformationGradient( asp::aspBase &asp,
                                                              const floatVector &deformationGradient ){

                    asp._localDeformationGradient.first = true;
                    asp._localDeformationGradient.second = deformationGradient;

                    return NULL;

                }

                static errorOut set_localMicroDeformation( asp::aspBase &asp,
                                                           const floatVector &microDeformation ){

                    asp._localMicroDeformation.first = true;
                    asp._localMicroDeformation.second = microDeformation;

                    return NULL;

                }

                static errorOut set_nonlocalMicroDeformation( asp::aspBase &asp,
                                                              const floatVector &nonlocalMicroDeformation ){

                    asp._nonlocalMicroDeformation.first = true;
                    asp._nonlocalMicroDeformation.second = nonlocalMicroDeformation;

                    return NULL;

                }

                static errorOut set_localSurfaceReferenceRelativePositionVector( asp::aspBase &asp,
                                                                                 const floatVector &localXi ){

                    asp._localSurfaceReferenceRelativePositionVector.first = true;
                    asp._localSurfaceReferenceRelativePositionVector.second = localXi;

                    return NULL;

                }

                static errorOut set_nonlocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp,
                                                                                    const floatVector &nonlocalXi ){

                    asp._nonlocalSurfaceReferenceRelativePositionVector.first = true;
                    asp._nonlocalSurfaceReferenceRelativePositionVector.second = nonlocalXi;

                    return NULL;

                }

                static errorOut set_referenceDistanceVector( asp::aspBase &asp, const floatVector &referenceDistanceVector ){

                    asp._referenceDistanceVector.first = true;
                    asp._referenceDistanceVector.second = referenceDistanceVector;

                    return NULL;

                }

                static errorOut set_gradientMicroDeformation( asp::aspBase &asp, const floatVector &gradientMicroDeformation ){

                    asp._gradientMicroDeformation = gradientMicroDeformation;

                    return NULL;

                }

                static errorOut set_localReferenceParticleSpacing( asp::aspBase &asp, const floatVector &referenceParticleSpacing ){

                    asp._localReferenceParticleSpacing.first = true;
                    asp._localReferenceParticleSpacing.second = referenceParticleSpacing;

                    return NULL;

                }

                static errorOut set_surfaceParameters( asp::aspBase &asp, const floatVector &surfaceParameters ){

                    asp._surfaceParameters.first = true;
                    asp._surfaceParameters.second = surfaceParameters;

                    return NULL;

                }

                static errorOut set_currentDistance( asp::aspBase &asp, const floatVector &currentDistance ){

                    asp._currentDistanceVector.first = true;
                    asp._currentDistanceVector.second = currentDistance;

                    return NULL;

                }

                static errorOut set_localCurrentNormal( asp::aspBase &asp, const floatVector &localCurrentNormal ){

                    asp._localCurrentNormal.first = true;
                    asp._localCurrentNormal.second = localCurrentNormal;

                    return NULL;

                }


                // Read functions for checking for errors
                static std::pair< bool, floatVector > getLocalReferenceNormal( asp::aspBase &asp ){

                    return asp._localReferenceNormal;

                }

                static std::pair< bool, floatVector > getLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp ){

                    return asp._localSurfaceReferenceRelativePositionVector;

                }

                static std::pair< bool, floatVector > getNonLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp ){

                    return asp._nonlocalSurfaceReferenceRelativePositionVector;

                }

                static std::pair< bool, floatVector > getReferenceDistanceVector( asp::aspBase &asp ){

                    return asp._referenceDistanceVector;

                }

                static std::pair< bool, floatVector > getLocalDeformationGradient( asp::aspBase &asp ){

                    return asp._localDeformationGradient;

                }

                static std::pair< bool, floatVector > getLocalMicroDeformation( asp::aspBase &asp ){

                    return asp._localMicroDeformation;

                }

                static std::pair< bool, floatVector > getNonLocalMicroDeformation( asp::aspBase &asp ){

                    return asp._nonlocalMicroDeformation;

                }

                static std::pair< bool, floatVector > getCurrentDistance( asp::aspBase &asp ){

                    return asp._currentDistanceVector;

                }

                static std::pair< bool, floatVector > getLocalCurrentNormal( asp::aspBase &asp ){

                    return asp._localCurrentNormal;

                }

                static std::pair< bool, floatVector > getSurfaceParameters( asp::aspBase &asp ){

                    return asp._surfaceParameters;

                }

        };

    }

}

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
    /*!
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

BOOST_AUTO_TEST_CASE( test_aspBase_initializeUnitSphere ){
    /*!
     * Test the initialization of the unit sphere
     */

    asp::aspBase asp;

    std::pair< bool, floatVector > unitSpherePoints;

    std::pair< bool, std::vector< unsigned int > > unitSphereConnectivity;

    asp::unit_test::aspBaseTester::initializeUnitSphere( asp, unitSpherePoints, unitSphereConnectivity );

    BOOST_CHECK( unitSpherePoints.first );

    BOOST_CHECK( unitSphereConnectivity.first );

    unsigned int npoints = unitSpherePoints.second.size( ) / 3;

    BOOST_CHECK( ( unitSpherePoints.second.size( ) % 3 ) == 0 );

    auto it = std::min_element( unitSphereConnectivity.second.begin( ), unitSphereConnectivity.second.end( ) );

    BOOST_CHECK( ( *it ) == 0 );

    it = std::max_element( unitSphereConnectivity.second.begin( ), unitSphereConnectivity.second.end( ) );

    BOOST_CHECK( ( *it ) == ( npoints - 1 ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalReferenceRadius ){

    class aspBaseMock : public asp::aspBase{

        public:

            aspBaseMock( floatType &radius ){
    
                asp::unit_test::aspBaseTester::set_radius( *this, radius );
    
            }

    };

    floatType radiusAnswer = 2.3;

    aspBaseMock asp( radiusAnswer );

    std::pair< bool, floatType > radius;

    asp::unit_test::aspBaseTester::setLocalReferenceRadius( asp, radius );

    BOOST_CHECK( radius.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( radiusAnswer, radius.second ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalReferenceRadius ){

    class aspBaseMock : public asp::aspBase{

        public:

            aspBaseMock( floatType &radius ){
    
                asp::unit_test::aspBaseTester::set_radius( *this, radius );
    
            }

    };

    floatType radiusAnswer = 2.3;

    aspBaseMock asp( radiusAnswer );

    std::pair< bool, floatType > radius;

    asp::unit_test::aspBaseTester::setNonLocalReferenceRadius( asp, radius );

    BOOST_CHECK( radius.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( radiusAnswer, radius.second ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalReferenceNormal ){

    class aspBaseMock : public asp::aspBase{

        errorOut initializeUnitSphere( ){

            floatVector points = {
                                      1, 2, 3,
                                      4, 5, 6,
                                      7, 8, 9
                                  };

            std::vector< unsigned int > connectivity = {
                                                           1, 2, 3
                                                       };

            asp::unit_test::aspBaseTester::set_unitSphere( *this, points, connectivity );

            return NULL;

        }

    };

    aspBaseMock asp;

    floatMatrix normalAnswers = {
                                    { 1, 2, 3},
                                    { 4, 5, 6},
                                    { 7, 8, 9}
                                };

    for ( unsigned int i = 0; i < normalAnswers.size( ); i++ ){

        normalAnswers[ i ] /= vectorTools::l2norm( normalAnswers[ i ] );

    }

    std::pair< bool, floatVector > result;

    unsigned int localIndex = 0;

    unsigned int nonlocalIndex = 0;

    for ( unsigned int i = 0; i < normalAnswers.size( ); i++ ){

        asp::unit_test::aspBaseTester::set_indices( asp, localIndex, nonlocalIndex, i );

        asp::unit_test::aspBaseTester::setLocalReferenceNormal( asp, result );

        BOOST_CHECK( result.first );

        BOOST_CHECK( vectorTools::fuzzyEquals( result.second, normalAnswers[ i ] ) );

    }

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalSurfaceReferenceRelativePositionVector ){

    class aspBaseMock : public asp::aspBase{

        errorOut setLocalReferenceNormal( ){

            floatVector normal = { 1, 2, 3 };

            normal /= vectorTools::l2norm( normal );

            asp::unit_test::aspBaseTester::set_localReferenceNormal( *this, normal );

            return NULL;

        }

        errorOut setLocalReferenceRadius( ){

            floatType radius = 2.45;

            asp::unit_test::aspBaseTester::set_localReferenceRadius( *this, radius );

            return NULL;

        }

    };

    aspBaseMock asp;

    floatVector answer = { 0.65479004, 1.30958009, 1.96437013 };

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setLocalSurfaceReferenceRelativePositionVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalSurfaceReferenceRelativePositionVector ){

    class aspBaseMock : public asp::aspBase{

        errorOut setLocalReferenceNormal( ){

            floatVector normal = { 1, 2, 3 };

            normal /= vectorTools::l2norm( normal );

            asp::unit_test::aspBaseTester::set_localReferenceNormal( *this, normal );

            return NULL;

        }

        errorOut setNonLocalReferenceRadius( ){

            floatType radius = 2.45;

            asp::unit_test::aspBaseTester::set_nonlocalReferenceRadius( *this, radius );

            return NULL;

        }

    };

    aspBaseMock asp;

    floatVector answer = { -0.65479004, -1.30958009, -1.96437013 };

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setNonLocalSurfaceReferenceRelativePositionVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalDeformationGradient ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        
    asp::unit_test::aspBaseTester::set_deformationGradient( asp, answer );

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setLocalDeformationGradient( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalMicroDeformation ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        
    asp::unit_test::aspBaseTester::set_microDeformation( asp, answer );

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setLocalMicroDeformation( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setReferenceDistanceVector ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 0, 0, 0 };

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setReferenceDistanceVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalReferenceParticleSpacing ){

    class aspBaseMock : public asp::aspBase{

        errorOut setLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localSurfaceReferenceRelativePositionVector( *this, value );

            return NULL;

        }

        errorOut setNonLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 4, 5, 6 };

            asp::unit_test::aspBaseTester::set_nonlocalSurfaceReferenceRelativePositionVector( *this, value );

            return NULL;

        }

        errorOut setReferenceDistanceVector( ){

            floatVector value = { 7, 8, 9 };

            asp::unit_test::aspBaseTester::set_referenceDistanceVector( *this, value );

            return NULL;

        }

    };

    aspBaseMock asp;

    floatVector answer = { 4, 5, 6 };

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setLocalReferenceParticleSpacing( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalMicroDeformation ){

    class aspBaseMock : public asp::aspBase{

        errorOut setLocalReferenceParticleSpacing( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localReferenceParticleSpacing( *this, value );

            return NULL;

        }

    };

    aspBaseMock asp;

    floatVector microDeformation = { 4,  5,  6,
                                     7,  8,  9,
                                    10, 11, 12 };

    asp::unit_test::aspBaseTester::set_microDeformation( asp, microDeformation );

    floatVector gradientMicroDeformation = { 13, 14, 15, 16, 17, 18, 19, 20, 21,
                                             22, 23, 24, 25, 26, 27, 28, 29, 30,
                                             31, 32, 33, 34, 35, 36, 37, 38, 39 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp, gradientMicroDeformation );

    floatVector answer = { 90, 109, 128,
                          147, 166, 185, 
                          204, 223, 242 };

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setNonLocalMicroDeformation( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalCurrentNormal ){

    class aspBaseMock : public asp::aspBase{

        errorOut setLocalReferenceNormal( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localReferenceNormal( *this, value );

            return NULL;

        }

        errorOut setLocalMicroDeformation( ){

            floatVector value = { 0.39293837, -0.42772133, -0.54629709,
                                  0.10262954,  0.43893794, -0.15378708,
                                  0.9615284 ,  0.36965948, -0.0381362 };

            asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, value );

            return NULL;

        }

    };

    aspBaseMock asp;

    floatVector answer = { -0.73381014, -0.454509  ,  0.50492004 };

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setLocalCurrentNormal( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setCurrentDistance ){

    class aspBaseMock : public asp::aspBase{

        errorOut setLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localSurfaceReferenceRelativePositionVector( *this, value );

            return NULL;

        }

        errorOut setNonLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 4, 5, 6 };

            asp::unit_test::aspBaseTester::set_nonlocalSurfaceReferenceRelativePositionVector( *this, value );

            return NULL;

        }

        errorOut setReferenceDistanceVector( ){

            floatVector value = { 7, 8, 9 };

            asp::unit_test::aspBaseTester::set_referenceDistanceVector( *this, value );

            return NULL;

        }

        errorOut setLocalDeformationGradient( ){

            floatVector value = { 10, 11, 12,
                                  13, 14, 15,
                                  16, 17, 18 };

            asp::unit_test::aspBaseTester::set_localDeformationGradient( *this, value );

            return NULL;

        }

        errorOut setLocalMicroDeformation( ){

            floatVector value = { 19, 20, 21,
                                  22, 23, 24,
                                  25, 26, 27 };

            asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, value );

            return NULL;

        }

        errorOut setNonLocalMicroDeformation( ){

            floatVector value = { 28, 29, 30,
                                  31, 32, 33,
                                  34, 35, 36 };

            asp::unit_test::aspBaseTester::set_nonlocalMicroDeformation( *this, value );

            return NULL;

        }

    };

    aspBaseMock asp;

    floatVector answer = { 482, 554, 626 };

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setCurrentDistance( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setSurfaceParameters ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer( 0, 0 );

    std::pair< bool, floatVector > result;

    asp::unit_test::aspBaseTester::setSurfaceParameters( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_initializeSurfaceIntegrandQuantities ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector localReferenceNormalAnswer = { 1, 2, 3};

            floatVector localSurfaceReferenceRelativePositionVectorAnswer = { 4, 5, 6 };

            floatVector nonlocalSurfaceReferenceRelativePositionVectorAnswer = { 7, 8, 9 };

            floatVector referenceDistanceVectorAnswer = { 10, 11, 12 };

            floatVector localDeformationGradientAnswer = { 13, 14, 15,
                                                           16, 17, 18,
                                                           19, 20, 21 };

            floatVector localMicroDeformationAnswer = { 22, 23, 24,
                                                        25, 26, 27,
                                                        28, 29, 30 };

            floatVector nonlocalMicroDeformationAnswer = { 31, 32, 33,
                                                           34, 35, 36,
                                                           37, 38, 39 };

            floatVector currentDistanceAnswer = { 40, 41, 42 };

            floatVector localCurrentNormalAnswer = { 43, 44, 45 };

            floatVector surfaceParametersAnswer = { 46, 47 };

        private:

            virtual errorOut setLocalReferenceNormal( ){
    
                asp::unit_test::aspBaseTester::set_localReferenceNormal( *this, localReferenceNormalAnswer );

                return NULL;
    
            }

            virtual errorOut setLocalSurfaceReferenceRelativePositionVector( ){

                asp::unit_test::aspBaseTester::set_localSurfaceReferenceRelativePositionVector( *this, localSurfaceReferenceRelativePositionVectorAnswer );

                return NULL;

            }

            virtual errorOut setNonLocalSurfaceReferenceRelativePositionVector( ){

                asp::unit_test::aspBaseTester::set_nonlocalSurfaceReferenceRelativePositionVector( *this, nonlocalSurfaceReferenceRelativePositionVectorAnswer );

                return NULL;

            }

            virtual errorOut setReferenceDistanceVector( ){

                asp::unit_test::aspBaseTester::set_referenceDistanceVector( *this, referenceDistanceVectorAnswer );

                return NULL;

            }

            virtual errorOut setLocalDeformationGradient( ){

                asp::unit_test::aspBaseTester::set_localDeformationGradient( *this, localDeformationGradientAnswer );

                return NULL;

            }

            virtual errorOut setLocalMicroDeformation( ){

                asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, localMicroDeformationAnswer );

                return NULL;

            }

            virtual errorOut setNonLocalMicroDeformation( ){

                asp::unit_test::aspBaseTester::set_nonlocalMicroDeformation( *this, nonlocalMicroDeformationAnswer );

                return NULL;

            }

            virtual errorOut setCurrentDistance( ){

                asp::unit_test::aspBaseTester::set_currentDistance( *this, currentDistanceAnswer );

                return NULL;

            }

            virtual errorOut setLocalCurrentNormal( ){

                asp::unit_test::aspBaseTester::set_localCurrentNormal( *this, localCurrentNormalAnswer );

                return NULL;

            }

            virtual errorOut setSurfaceParameters( ){

                asp::unit_test::aspBaseTester::set_surfaceParameters( *this, surfaceParametersAnswer );

                return NULL;

            }

    };

    aspBaseMock asp;

    asp::unit_test::aspBaseTester::initializeSurfaceIntegrandQuantities( asp );

    std::pair< bool, floatVector > result;

    result = asp::unit_test::aspBaseTester::getLocalReferenceNormal( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.localReferenceNormalAnswer ) );

    result = asp::unit_test::aspBaseTester::getLocalSurfaceReferenceRelativePositionVector( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.localSurfaceReferenceRelativePositionVectorAnswer ) );

    result = asp::unit_test::aspBaseTester::getNonLocalSurfaceReferenceRelativePositionVector( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.nonlocalSurfaceReferenceRelativePositionVectorAnswer ) );

    result = asp::unit_test::aspBaseTester::getReferenceDistanceVector( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.referenceDistanceVectorAnswer ) );

    result = asp::unit_test::aspBaseTester::getLocalDeformationGradient( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.localDeformationGradientAnswer ) );

    result = asp::unit_test::aspBaseTester::getLocalMicroDeformation( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.localMicroDeformationAnswer ) );

    result = asp::unit_test::aspBaseTester::getNonLocalMicroDeformation( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.nonlocalMicroDeformationAnswer ) );

    result = asp::unit_test::aspBaseTester::getCurrentDistance( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.currentDistanceAnswer ) );

    result = asp::unit_test::aspBaseTester::getLocalCurrentNormal( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.localCurrentNormalAnswer ) );

    result = asp::unit_test::aspBaseTester::getSurfaceParameters( asp );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.surfaceParametersAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceEnergyDensity ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector distanceVector = { 1, 2, 3 };

            floatVector localCurrentNormal = { 0.45584231, 0.56980288, 0.68376346 };

            floatVector surfaceParameters = { 12.3, 45.6 };

            aspBaseMock( ){

                asp::unit_test::aspBaseTester::set_currentDistance( *this, distanceVector );

                asp::unit_test::aspBaseTester::set_localCurrentNormal( *this, localCurrentNormal );

                asp::unit_test::aspBaseTester::set_surfaceParameters( *this, surfaceParameters );

            }

    };
    
    aspBaseMock asp;

    floatType answer = 356.565771656861;

    floatType result;

    asp::unit_test::aspBaseTester::computeSurfaceEnergyDensity( asp, result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}
