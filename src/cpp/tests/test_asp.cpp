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

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer )
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

// Unit tester to open private members of aspBase for unit testing
namespace asp{

    namespace unit_test{

        class aspBaseTester{

            public:

                static void initializeUnitSphere( asp::aspBase &asp,
                                                      asp::dataStorage< floatVector > & unitSpherePoints,
                                                      asp::dataStorage< std::vector< unsigned int > > & unitSphereConnectivity ){
    
                    BOOST_CHECK_NO_THROW( asp.initializeUnitSphere( ) );
    
                    unitSpherePoints = asp._unitSpherePoints;
    
                    unitSphereConnectivity = asp._unitSphereConnectivity;
    
                    return;
    
                }

                static void setLocalReferenceRadius( asp::aspBase &asp, asp::dataStorage< floatType > &radius ){

                    BOOST_CHECK_NO_THROW( asp.setLocalReferenceRadius( ) );

                    radius = asp._localReferenceRadius;

                    return;

                }

                static void setNonLocalReferenceRadius( asp::aspBase &asp, asp::dataStorage< floatType > &radius ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalReferenceRadius( ) );

                    radius = asp._nonLocalReferenceRadius;

                    return;

                }

                static void setLocalReferenceNormal( asp::aspBase &asp, asp::dataStorage< floatVector > &normal ){

                    BOOST_CHECK_NO_THROW( asp.setLocalReferenceNormal( ) );

                    normal = asp._localReferenceNormal;

                    return;

                }

                static void setLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp, asp::dataStorage< floatVector > &Xi ){

                    BOOST_CHECK_NO_THROW( asp.setLocalSurfaceReferenceRelativePositionVector( ) );

                    Xi = asp._localSurfaceReferenceRelativePositionVector;

                    return;

                }

                static void setNonLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp, asp::dataStorage< floatVector > &Xi ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalSurfaceReferenceRelativePositionVector( ) );

                    Xi = asp._nonLocalSurfaceReferenceRelativePositionVector;

                    return;

                }

                static void setLocalDeformationGradient( asp::aspBase &asp, asp::dataStorage< floatVector > &localDeformationGradient ){

                    BOOST_CHECK_NO_THROW( asp.setLocalDeformationGradient( ) );

                    localDeformationGradient = asp._localDeformationGradient;

                    return;

                }

                static void setLocalMicroDeformation( asp::aspBase &asp, asp::dataStorage< floatVector > &localMicroDeformation ){

                    BOOST_CHECK_NO_THROW( asp.setLocalMicroDeformation( ) );

                    localMicroDeformation = asp._localMicroDeformation;

                    return;

                }

                static void setReferenceDistanceVector( asp::aspBase &asp, asp::dataStorage< floatVector > &referenceDistanceVector ){

                    BOOST_CHECK_NO_THROW( asp.setReferenceDistanceVector( ) );

                    referenceDistanceVector = asp._referenceDistanceVector;

                    return;

                }

                static void setLocalReferenceParticleSpacing( asp::aspBase &asp,
                                                                  asp::dataStorage< floatVector > &localReferenceParticleSpacing ){

                    BOOST_CHECK_NO_THROW( asp.setLocalReferenceParticleSpacing( ) );

                    localReferenceParticleSpacing = asp._localReferenceParticleSpacing;

                    return;

                }

                static void setNonLocalMicroDeformation( asp::aspBase &asp,
                                                             asp::dataStorage< floatVector > &nonLocalMicroDeformation ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalMicroDeformation( ) );

                    nonLocalMicroDeformation = asp._nonLocalMicroDeformation;

                    return;

                }

                static void setLocalCurrentNormal( asp::aspBase &asp,
                                                       asp::dataStorage< floatVector > &localCurrentNormal ){

                    BOOST_CHECK_NO_THROW( asp.setLocalCurrentNormal( ) );

                    localCurrentNormal = asp._localCurrentNormal;

                    return;

                }

                static void setCurrentDistanceVector( asp::aspBase &asp,
                                                    asp::dataStorage< floatVector > &currentDistance ){

                    BOOST_CHECK_NO_THROW( asp.setCurrentDistanceVector( ) );

                    currentDistance = asp._currentDistanceVector;

                    return;

                }

                static void setSurfaceParameters( asp::aspBase &asp,
                                                      asp::dataStorage< floatVector > &currentSurfaceParameters ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceParameters( ) );

                    currentSurfaceParameters = asp._surfaceParameters;

                    return;

                }

                static void formBoundingBox( asp::aspBase &asp, floatVector &points, floatMatrix &boundingBox ){

                    BOOST_CHECK_NO_THROW( asp.formBoundingBox( points, boundingBox ) );

                    return;

                }

                static bool pointInBoundingBox( asp::aspBase &asp, floatVector &point, floatMatrix &boundingBox ){

                    return asp.pointInBoundingBox( point, boundingBox );

                }

                static void idBoundingBoxContainedPoints( asp::aspBase &asp, floatVector &points, floatMatrix &boundingBox, std::vector< unsigned int > &containedPoints ){

                    BOOST_CHECK_NO_THROW( asp.idBoundingBoxContainedPoints( points, boundingBox, containedPoints ) );

                    return;

                }

                // Direct write functions for mocking
                static void set_indices( asp::aspBase &asp,
                                            unsigned int localIndex, unsigned int nonLocalIndex, unsigned int localSurfaceNodeIndex ){

                    asp._localIndex = localIndex;

                    asp._nonLocalIndex = nonLocalIndex;

                    asp._localSurfaceNodeIndex = localSurfaceNodeIndex;

                    return;

                }

                static void set_radius( asp::aspBase &asp,
                                            floatType &radius ){

                    asp._radius = radius;

                    return;

                }

                static void set_unitSphere( asp::aspBase &asp,
                                                const floatVector &points, std::vector< unsigned int > &connectivity ){

                    asp._unitSpherePoints.first = true;
                    asp._unitSpherePoints.second = points;

                    asp._unitSphereConnectivity.first = true;
                    asp._unitSphereConnectivity.second = connectivity;

                    return;

                }

                static void set_localReferenceNormal( asp::aspBase &asp,
                                                          const floatVector &normal ){

                    asp._localReferenceNormal.first = true;
                    asp._localReferenceNormal.second = normal;

                    return;

                }

                static void set_localReferenceRadius( asp::aspBase &asp,
                                                          const floatType &radius ){

                    asp._localReferenceRadius.first = true;
                    asp._localReferenceRadius.second = radius;

                    return;

                }

                static void set_nonLocalReferenceRadius( asp::aspBase &asp,
                                                             const floatType &radius ){

                    asp._nonLocalReferenceRadius.first = true;
                    asp._nonLocalReferenceRadius.second = radius;

                    return;

                }

                static void set_deformationGradient( asp::aspBase &asp,
                                                         const floatVector &deformationGradient ){

                    asp._deformationGradient = deformationGradient;

                    return;

                }

                static void set_microDeformation( asp::aspBase &asp,
                                                      const floatVector &microDeformation ){

                    asp._microDeformation = microDeformation;

                    return;

                }

                static void set_localDeformationGradient( asp::aspBase &asp,
                                                              const floatVector &deformationGradient ){

                    asp._localDeformationGradient.first = true;
                    asp._localDeformationGradient.second = deformationGradient;

                    return;

                }

                static void set_localMicroDeformation( asp::aspBase &asp,
                                                           const floatVector &microDeformation ){

                    asp._localMicroDeformation.first = true;
                    asp._localMicroDeformation.second = microDeformation;

                    return;

                }

                static void set_nonLocalMicroDeformation( asp::aspBase &asp,
                                                              const floatVector &nonLocalMicroDeformation ){

                    asp._nonLocalMicroDeformation.first = true;
                    asp._nonLocalMicroDeformation.second = nonLocalMicroDeformation;

                    return;

                }

                static void set_localSurfaceReferenceRelativePositionVector( asp::aspBase &asp,
                                                                                 const floatVector &localXi ){

                    asp._localSurfaceReferenceRelativePositionVector.first = true;
                    asp._localSurfaceReferenceRelativePositionVector.second = localXi;

                    return;

                }

                static void set_nonLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp,
                                                                                    const floatVector &nonLocalXi ){

                    asp._nonLocalSurfaceReferenceRelativePositionVector.first = true;
                    asp._nonLocalSurfaceReferenceRelativePositionVector.second = nonLocalXi;

                    return;

                }

                static void set_referenceDistanceVector( asp::aspBase &asp, const floatVector &referenceDistanceVector ){

                    asp._referenceDistanceVector.first = true;
                    asp._referenceDistanceVector.second = referenceDistanceVector;

                    return;

                }

                static void set_gradientMicroDeformation( asp::aspBase &asp, const floatVector &gradientMicroDeformation ){

                    asp._gradientMicroDeformation = gradientMicroDeformation;

                    return;

                }

                static void set_localReferenceParticleSpacing( asp::aspBase &asp, const floatVector &referenceParticleSpacing ){

                    asp._localReferenceParticleSpacing.first = true;
                    asp._localReferenceParticleSpacing.second = referenceParticleSpacing;

                    return;

                }

                static void set_surfaceParameters( asp::aspBase &asp, const floatVector &surfaceParameters ){

                    asp._surfaceParameters.first = true;
                    asp._surfaceParameters.second = surfaceParameters;

                    return;

                }

                static void set_currentDistance( asp::aspBase &asp, const floatVector &currentDistance ){

                    asp._currentDistanceVector.first = true;
                    asp._currentDistanceVector.second = currentDistance;

                    return;

                }

                static void set_localCurrentNormal( asp::aspBase &asp, const floatVector &localCurrentNormal ){

                    asp._localCurrentNormal.first = true;
                    asp._localCurrentNormal.second = localCurrentNormal;

                    return;

                }

                static void set_localReferenceSurfacePoints( asp::aspBase &asp, const floatVector &localReferenceSurfacePoints ){

                    asp._localReferenceSurfacePoints.first = true;
                    asp._localReferenceSurfacePoints.second = localReferenceSurfacePoints;

                    return;

                }

                static void set_nonLocalReferenceSurfacePoints( asp::aspBase &asp, const floatVector &nonLocalReferenceSurfacePoints ){

                    asp._nonLocalReferenceSurfacePoints.first = true;
                    asp._nonLocalReferenceSurfacePoints.second = nonLocalReferenceSurfacePoints;

                    return;

                }

                static void set_localCurrentSurfacePoints( asp::aspBase &asp, const floatVector &localCurrentSurfacePoints ){

                    asp._localCurrentSurfacePoints.first = true;
                    asp._localCurrentSurfacePoints.second = localCurrentSurfacePoints;

                    return;

                }

                static void set_nonLocalCurrentSurfacePoints( asp::aspBase &asp, const floatVector &nonLocalCurrentSurfacePoints ){

                    asp._nonLocalCurrentSurfacePoints.first = true;
                    asp._nonLocalCurrentSurfacePoints.second = nonLocalCurrentSurfacePoints;

                    return;

                }

                static void set_surfaceOverlapParameters( asp::aspBase &asp, const floatVector &nonLocalCurrentSurfacePoints ){

                    asp._surfaceOverlapParameters.first = true;
                    asp._surfaceOverlapParameters.second = nonLocalCurrentSurfacePoints;

                    return;

                }

                static void set_nonLocalParticleCurrentBoundingBox( asp::aspBase &asp, const floatMatrix &nonLocalParticleBoundingBox ){

                    asp._nonLocalParticleCurrentBoundingBox.first = true;
                    asp._nonLocalParticleCurrentBoundingBox.second = nonLocalParticleBoundingBox;

                    return;

                }

                static void set_localGradientMicroDeformation( asp::aspBase &asp, const floatVector &localGradientMicroDeformation ){

                    asp._localGradientMicroDeformation.first = true;
                    asp._localGradientMicroDeformation.second = localGradientMicroDeformation;

                    return;

                }

                static void set_nonLocalMicroDeformationBase( asp::aspBase &asp, const floatVector &nonLocalMicroDeformationBase ){

                    asp._nonLocalMicroDeformationBase.first = true;
                    asp._nonLocalMicroDeformationBase.second = nonLocalMicroDeformationBase;

                    return;

                }

                static void set_particlePairOverlap( asp::aspBase &asp, const std::unordered_map< unsigned int, floatVector > &particlePairOverlap ){

                    asp._particlePairOverlap.first = true;
                    asp._particlePairOverlap.second = particlePairOverlap;

                    return;

                }

                // Read functions for checking for errors
                static asp::dataStorage< floatVector > getLocalReferenceNormal( asp::aspBase &asp ){

                    return asp._localReferenceNormal;

                }

                static asp::dataStorage< floatVector > getLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp ){

                    return asp._localSurfaceReferenceRelativePositionVector;

                }

                static asp::dataStorage< floatVector > getNonLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp ){

                    return asp._nonLocalSurfaceReferenceRelativePositionVector;

                }

                static asp::dataStorage< floatVector > getReferenceDistanceVector( asp::aspBase &asp ){

                    return asp._referenceDistanceVector;

                }

                static asp::dataStorage< floatVector > getLocalDeformationGradient( asp::aspBase &asp ){

                    return asp._localDeformationGradient;

                }

                static asp::dataStorage< floatVector > getLocalMicroDeformation( asp::aspBase &asp ){

                    return asp._localMicroDeformation;

                }

                static asp::dataStorage< floatVector > getNonLocalMicroDeformation( asp::aspBase &asp ){

                    return asp._nonLocalMicroDeformation;

                }

                static asp::dataStorage< floatVector > getCurrentDistance( asp::aspBase &asp ){

                    return asp._currentDistanceVector;

                }

                static asp::dataStorage< floatVector > getLocalCurrentNormal( asp::aspBase &asp ){

                    return asp._localCurrentNormal;

                }

                static asp::dataStorage< floatVector > getSurfaceParameters( asp::aspBase &asp ){

                    return asp._surfaceParameters;

                }

                static std::vector< asp::dataBase* > getInteractionPairData( asp::aspBase &asp ){

                    return asp._interactionPairData;

                }

                static std::vector< asp::dataBase* > getSurfacePointData( asp::aspBase &asp ){

                    return asp._surfacePointData;

                }

                static std::vector< asp::dataBase* > getLocalParticleData( asp::aspBase &asp ){

                    return asp._localParticleData;

                }

                static void resetInteractionPairData( asp::aspBase & asp ){

                    asp.resetInteractionPairData( );

                }

                static void resetSurfacePointData( asp::aspBase & asp ){

                    asp.resetSurfacePointData( );

                }

                static void resetLocalParticleData( asp::aspBase & asp ){

                    asp.resetLocalParticleData( );

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

    BOOST_CHECK_NO_THROW( asp.computeLocalParticleEnergyDensity( previousTime, deltaTime, currentDeformationGradient, previousDeformationGradient,
                                                                 currentTemperature, previousTemperature, previousStateVariables, parameters, energy, cauchyStress ) );

    BOOST_CHECK_NO_THROW( vectorTools::fuzzyEquals( energy, 0. ) );

    BOOST_CHECK( cauchyStress.size( ) == currentDeformationGradient.size( ) );

    energy = 0;

    cauchyStress.clear( );

    floatType logProbabilityRatio = 3.4;

    BOOST_CHECK_NO_THROW( asp.computeLocalParticleEnergyDensity( previousTime, deltaTime, currentDeformationGradient, previousDeformationGradient,
                                                         currentTemperature, previousTemperature, previousStateVariables, parameters, energy, cauchyStress, logProbabilityRatio ) );

    BOOST_CHECK_NO_THROW( vectorTools::fuzzyEquals( energy, 0. ) );

    BOOST_CHECK( cauchyStress.size( ) == currentDeformationGradient.size( ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( logProbabilityRatio, 0. ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_initializeUnitSphere ){
    /*!
     * Test the initialization of the unit sphere
     */

    asp::aspBase asp;

    asp::dataStorage< floatVector > unitSpherePoints;

    asp::dataStorage< std::vector< unsigned int > > unitSphereConnectivity;

    asp::unit_test::aspBaseTester::initializeUnitSphere( asp, unitSpherePoints, unitSphereConnectivity );

    BOOST_CHECK( unitSpherePoints.first );

    BOOST_CHECK( unitSphereConnectivity.first );

    unsigned int npoints = unitSpherePoints.second.size( ) / 3;

    BOOST_CHECK( ( unitSpherePoints.second.size( ) % 3 ) == 0 );

    auto it = std::min_element( unitSphereConnectivity.second.begin( ), unitSphereConnectivity.second.end( ) );

    BOOST_CHECK( ( *it ) == 0 );

    it = std::max_element( unitSphereConnectivity.second.begin( ), unitSphereConnectivity.second.end( ) );

    BOOST_CHECK( ( *it ) == ( npoints - 1 ) );

    asp::aspBase aspGet1;
    asp::aspBase aspGet2;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet1.getUnitSpherePoints( ), unitSpherePoints.second ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet2.getUnitSphereConnectivity( ), unitSphereConnectivity.second ) );

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

    asp::dataStorage< floatType > radius;

    asp::unit_test::aspBaseTester::setLocalReferenceRadius( asp, radius );

    BOOST_CHECK( radius.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( radiusAnswer, radius.second ) );

    aspBaseMock aspGet( radiusAnswer );

    BOOST_CHECK( vectorTools::fuzzyEquals( radiusAnswer, *aspGet.getLocalReferenceRadius( ) ) );

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

    asp::dataStorage< floatType > radius;

    asp::unit_test::aspBaseTester::setNonLocalReferenceRadius( asp, radius );

    BOOST_CHECK( radius.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( radiusAnswer, radius.second ) );

    aspBaseMock aspGet( radiusAnswer );

    BOOST_CHECK( vectorTools::fuzzyEquals( radiusAnswer, *aspGet.getNonLocalReferenceRadius( ) ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalReferenceNormal ){

    class aspBaseMock : public asp::aspBase{

        void initializeUnitSphere( ){

            floatVector points = {
                                      1, 2, 3,
                                      4, 5, 6,
                                      7, 8, 9
                                  };

            std::vector< unsigned int > connectivity = {
                                                           1, 2, 3
                                                       };

            asp::unit_test::aspBaseTester::set_unitSphere( *this, points, connectivity );

            return;

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

    asp::dataStorage< floatVector > result;

    unsigned int localIndex = 0;

    unsigned int nonLocalIndex = 0;

    for ( unsigned int i = 0; i < normalAnswers.size( ); i++ ){

        asp::unit_test::aspBaseTester::set_indices( asp, localIndex, nonLocalIndex, i );

        asp::unit_test::aspBaseTester::setLocalReferenceNormal( asp, result );

        BOOST_CHECK( result.first );

        BOOST_CHECK( vectorTools::fuzzyEquals( result.second, normalAnswers[ i ] ) );

        aspBaseMock aspGet;

        asp::unit_test::aspBaseTester::set_indices( aspGet, localIndex, nonLocalIndex, i );

        BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalReferenceNormal( ), normalAnswers[ i ] ) );

    }

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalSurfaceReferenceRelativePositionVector ){

    class aspBaseMock : public asp::aspBase{

        void setLocalReferenceNormal( ){

            floatVector normal = { 1, 2, 3 };

            normal /= vectorTools::l2norm( normal );

            asp::unit_test::aspBaseTester::set_localReferenceNormal( *this, normal );

            return;

        }

        void setLocalReferenceRadius( ){

            floatType radius = 2.45;

            asp::unit_test::aspBaseTester::set_localReferenceRadius( *this, radius );

            return;

        }

    };

    aspBaseMock asp;

    floatVector answer = { 0.65479004, 1.30958009, 1.96437013 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalSurfaceReferenceRelativePositionVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalSurfaceReferenceRelativePositionVector( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalSurfaceReferenceRelativePositionVector ){

    class aspBaseMock : public asp::aspBase{

        void setLocalReferenceNormal( ){

            floatVector normal = { 1, 2, 3 };

            normal /= vectorTools::l2norm( normal );

            asp::unit_test::aspBaseTester::set_localReferenceNormal( *this, normal );

            return;

        }

        void setNonLocalReferenceRadius( ){

            floatType radius = 2.45;

            asp::unit_test::aspBaseTester::set_nonLocalReferenceRadius( *this, radius );

            return;

        }

    };

    aspBaseMock asp;

    floatVector answer = { -0.65479004, -1.30958009, -1.96437013 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setNonLocalSurfaceReferenceRelativePositionVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getNonLocalSurfaceReferenceRelativePositionVector( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalDeformationGradient ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        
    asp::unit_test::aspBaseTester::set_deformationGradient( asp, answer );

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalDeformationGradient( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    asp::unit_test::aspBaseTester::set_deformationGradient( aspGet, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalDeformationGradient( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalMicroDeformation ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        
    asp::unit_test::aspBaseTester::set_microDeformation( asp, answer );

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalMicroDeformation( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    asp::unit_test::aspBaseTester::set_microDeformation( aspGet, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalMicroDeformation( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setReferenceDistanceVector ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 0, 0, 0 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setReferenceDistanceVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getReferenceDistanceVector( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalReferenceParticleSpacing ){

    class aspBaseMock : public asp::aspBase{

        void setLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localSurfaceReferenceRelativePositionVector( *this, value );

            return;

        }

        void setNonLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 4, 5, 6 };

            asp::unit_test::aspBaseTester::set_nonLocalSurfaceReferenceRelativePositionVector( *this, value );

            return;

        }

        void setReferenceDistanceVector( ){

            floatVector value = { 7, 8, 9 };

            asp::unit_test::aspBaseTester::set_referenceDistanceVector( *this, value );

            return;

        }

    };

    aspBaseMock asp;

    floatVector answer = { 4, 5, 6 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalReferenceParticleSpacing( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalReferenceParticleSpacing( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalMicroDeformation ){

    class aspBaseMock : public asp::aspBase{

        void setLocalReferenceParticleSpacing( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localReferenceParticleSpacing( *this, value );

            return;

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

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setNonLocalMicroDeformation( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    asp::unit_test::aspBaseTester::set_microDeformation( aspGet, microDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspGet, gradientMicroDeformation );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getNonLocalMicroDeformation( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalCurrentNormal ){

    class aspBaseMock : public asp::aspBase{

        void setLocalReferenceNormal( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localReferenceNormal( *this, value );

            return;

        }

        void setLocalMicroDeformation( ){

            floatVector value = { 0.39293837, -0.42772133, -0.54629709,
                                  0.10262954,  0.43893794, -0.15378708,
                                  0.9615284 ,  0.36965948, -0.0381362 };

            asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, value );

            return;

        }

    };

    aspBaseMock asp;

    floatVector answer = { -0.73381014, -0.454509  ,  0.50492004 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalCurrentNormal( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalCurrentNormal( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setCurrentDistanceVector ){

    class aspBaseMock : public asp::aspBase{

        void setLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localSurfaceReferenceRelativePositionVector( *this, value );

            return;

        }

        void setNonLocalSurfaceReferenceRelativePositionVector( ){

            floatVector value = { 4, 5, 6 };

            asp::unit_test::aspBaseTester::set_nonLocalSurfaceReferenceRelativePositionVector( *this, value );

            return;

        }

        void setReferenceDistanceVector( ){

            floatVector value = { 7, 8, 9 };

            asp::unit_test::aspBaseTester::set_referenceDistanceVector( *this, value );

            return;

        }

        void setLocalDeformationGradient( ){

            floatVector value = { 10, 11, 12,
                                  13, 14, 15,
                                  16, 17, 18 };

            asp::unit_test::aspBaseTester::set_localDeformationGradient( *this, value );

            return;

        }

        void setLocalMicroDeformation( ){

            floatVector value = { 19, 20, 21,
                                  22, 23, 24,
                                  25, 26, 27 };

            asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, value );

            return;

        }

        void setNonLocalMicroDeformation( ){

            floatVector value = { 28, 29, 30,
                                  31, 32, 33,
                                  34, 35, 36 };

            asp::unit_test::aspBaseTester::set_nonLocalMicroDeformation( *this, value );

            return;

        }

    };

    aspBaseMock asp;

    floatVector answer = { 482, 554, 626 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setCurrentDistanceVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getCurrentDistanceVector( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setSurfaceParameters ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector surfaceParametersAnswer = { 1, 2, 3 };

        private:

            virtual void setSurfaceParameters( ){

                asp::unit_test::aspBaseTester::set_surfaceParameters( *this, surfaceParametersAnswer );

            }
    };

    aspBaseMock asp;

    floatVector answer = { 1, 2, 3 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setSurfaceParameters( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getSurfaceParameters( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceAdhesionEnergyDensity ){

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

    floatType answer = 178.2828858284305;

    floatType result;

    asp.computeSurfaceAdhesionEnergyDensity( result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    result = *asp.getSurfaceAdhesionEnergyDensity( );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceAdhesionTraction ){

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
    
    aspBaseMock asp, aspGet;

    floatVector dn = vectorTools::dot( asp.distanceVector, asp.localCurrentNormal ) * asp.localCurrentNormal;

    floatVector dt = ( asp.distanceVector - dn );

    floatVector answer = asp.surfaceParameters[ 0 ] * dn + asp.surfaceParameters[ 1 ] * dt;

    floatVector result;

    asp.computeSurfaceAdhesionTraction( result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    result = *aspGet.getSurfaceAdhesionTraction( );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getSurfaceAdhesionEnergyDensity_error ){

    class aspBaseMock : public asp::aspBase{

        private:

            virtual void computeSurfaceAdhesionEnergyDensity( floatType &surfaceAdhesionEnergyDensity ){

                throw std::runtime_error( "This should throw an error" );

            }

    };

    aspBaseMock asp;

    //Setup redirect variables for stderr
    std::stringbuf buffer;
    cerr_redirect rd( &buffer );

    BOOST_REQUIRE_THROW(
        try{
            asp.getSurfaceAdhesionEnergyDensity( );
        }
        catch(std::exception &e){
            errorTools::printNestedExceptions( e );
            throw;
        }
    , std::exception );

    BOOST_CHECK( buffer.str( ).find( "This should throw an error" ) != std::string::npos );

    BOOST_CHECK( buffer.str( ).find( "getSurfaceAdhesionEnergyDensity" ) != std::string::npos );

    BOOST_CHECK( buffer.str( ).find( "setSurfaceAdhesionEnergyDensity" ) != std::string::npos );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getSurfaceAdhesionEnergyDensity ){

    class aspBaseMock : public asp::aspBase{

        private:

            virtual void computeSurfaceAdhesionEnergyDensity( floatType &surfaceAdhesionEnergyDensity ){

                surfaceAdhesionEnergyDensity = 123;

                return;

            }

    };

    aspBaseMock asp;

    //Setup redirect variables for stderr
    std::stringbuf buffer;
    cerr_redirect rd( &buffer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getSurfaceAdhesionEnergyDensity( ), 123. ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalReferenceSurfacePoints ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector points = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            std::vector< unsigned int > connectivity = { 1, 2, 3 };

            floatType radius = 1.3;

        private:

            virtual void initializeUnitSphere( ){

                asp::unit_test::aspBaseTester::set_unitSphere( *this, points, connectivity );

            }

            virtual void setLocalReferenceRadius( ){

                asp::unit_test::aspBaseTester::set_localReferenceRadius( *this, radius );

            }

    };

    aspBaseMock asp;

    floatVector answer = asp.points * asp.radius;

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getLocalReferenceSurfacePoints( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getNonLocalReferenceSurfacePoints ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector points = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            std::vector< unsigned int > connectivity = { 1, 2, 3 };

            floatType radius = 1.3;

        private:

            virtual void initializeUnitSphere( ){

                asp::unit_test::aspBaseTester::set_unitSphere( *this, points, connectivity );

            }

            virtual void setNonLocalReferenceRadius( ){

                asp::unit_test::aspBaseTester::set_nonLocalReferenceRadius( *this, radius );

            }

    };

    aspBaseMock asp;

    floatVector answer = asp.points * asp.radius;

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getNonLocalReferenceSurfacePoints( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalCurrentSurfacePoints ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector points = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

            floatVector microDeformation = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };

        private:

            virtual void setLocalReferenceSurfacePoints( ){

                asp::unit_test::aspBaseTester::set_localReferenceSurfacePoints( *this, points );

            }

            virtual void setLocalMicroDeformation( ){

                asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, microDeformation );

            }

    };

    aspBaseMock asp;

    floatVector answer = { 1.4,  3.2,  5. ,
                           3.2,  7.7, 12.2,
                           5. , 12.2, 19.4,
                           6.8, 16.7, 26.6};

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getLocalCurrentSurfacePoints( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getNonLocalCurrentSurfacePoints ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector points = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

            floatVector microDeformation = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };

        private:

            virtual void setNonLocalReferenceSurfacePoints( ){

                asp::unit_test::aspBaseTester::set_nonLocalReferenceSurfacePoints( *this, points );

            }

            virtual void setNonLocalMicroDeformation( ){

                asp::unit_test::aspBaseTester::set_nonLocalMicroDeformation( *this, microDeformation );

            }

    };

    aspBaseMock asp;

    floatVector answer = { 1.4,  3.2,  5. ,
                           3.2,  7.7, 12.2,
                           5. , 12.2, 19.4,
                           6.8, 16.7, 26.6};

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getNonLocalCurrentSurfacePoints( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalParticleCurrentBoundingBox ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector points = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

        private:

            virtual void setLocalCurrentSurfacePoints( ){

                asp::unit_test::aspBaseTester::set_localCurrentSurfacePoints( *this, points );

            }

    };

    aspBaseMock asp;

    floatMatrix answer = { { 1, 10 },
                           { 2, 11 },
                           { 3, 12 } };

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getLocalParticleCurrentBoundingBox( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getNonLocalParticleCurrentBoundingBox ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector points = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

        private:

            virtual void setNonLocalCurrentSurfacePoints( ){

                asp::unit_test::aspBaseTester::set_nonLocalCurrentSurfacePoints( *this, points );

            }

    };

    aspBaseMock asp;

    floatMatrix answer = { { 1, 10 },
                           { 2, 11 },
                           { 3, 12 } };

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getNonLocalParticleCurrentBoundingBox( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_formBoundingBox ){

    asp::aspBase asp;

    floatVector points = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    floatMatrix answer = { { 1, 7 },
                           { 2, 8 },
                           { 3, 9 } };

    floatMatrix result;

    asp::unit_test::aspBaseTester::formBoundingBox( asp, points, result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_pointInBoundingBox ){

    asp::aspBase asp;

    floatVector point = { 1, 2 };

    floatMatrix boundingBox = { { 0, 2 }, { 1, 3 } };

    BOOST_CHECK( asp::unit_test::aspBaseTester::pointInBoundingBox( asp, point, boundingBox ) );

    point = { 1, -2 };

    BOOST_CHECK( !asp::unit_test::aspBaseTester::pointInBoundingBox( asp, point, boundingBox ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_idBoundingBoxContainedPoints ){

    asp::aspBase asp;

    floatVector points = { 1, 2, 3, 1, -2, 4, 0.5, 1.1, 4.9 };

    floatMatrix boundingBox = { { 0, 2 }, { 1, 3 }, { 2, 5 } };

    std::vector< unsigned int > answer = { 0, 2 };

    std::vector< unsigned int > result;

    asp::unit_test::aspBaseTester::idBoundingBoxContainedPoints( asp, points, boundingBox, result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getMicroDeformation ){

    asp::aspBase asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6 };

    asp::unit_test::aspBaseTester::set_microDeformation( asp, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getMicroDeformation( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getGradientMicroDeformation ){

    asp::aspBase asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getGradientMicroDeformation( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalGradientMicroDeformation ){

    asp::aspBase asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getLocalGradientMicroDeformation( ), answer ) );

}


BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalMicroDeformationBase ){

    asp::aspBase asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6 };

    asp::unit_test::aspBaseTester::set_microDeformation( asp, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getNonLocalMicroDeformationBase( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setSurfaceOverlapParameters ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector overlapParameters = { 1, 2, 3};

        private:

            virtual void setSurfaceOverlapParameters( ){

                asp::unit_test::aspBaseTester::set_surfaceOverlapParameters( *this, overlapParameters );

            }

    };

    aspBaseMock asp;

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp.getSurfaceOverlapParameters( ), asp.overlapParameters ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setParticlePairOverlap ){
    /*!
     * Test setting the particle pair overlap values
     */

    class aspBaseMock : public asp::aspBase{

        public:

            floatMatrix nonLocalParticleCurrentBoundingBox = { { 0, 2 }, { 0, 2 }, { 0, 2 } };

            floatVector localReferenceSurfacePoints = { 0, 0, 1, 1, 0, 0, 1.1, 0, 0 };

            floatVector localCurrentSurfacePoints = { 0, 0, 1, 1.1, 0, 0, 1.1 * 1.1, 0, 0 };

            floatVector localDeformationGradient = { 1.1, 0, 0, 0, 1, 0, 0, 0, 1 };

            floatType nonLocalReferenceRadius = 2.3;

            floatVector localMicroDeformation = { 1.5, 0, 0, 0, 1, 0, 0, 0, 1 };

            floatVector localReferenceParticleSpacing = { 2. / 1.1, 0, 0 };

            floatVector localGradientMicroDeformation = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0 };

            floatVector nonLocalMicroDeformationBase = { 1 / 2.3, 0, 0, 0, 1, 0, 0, 0, 1 };

        private:

            virtual void setNonLocalParticleCurrentBoundingBox( ){

                asp::unit_test::aspBaseTester::set_nonLocalParticleCurrentBoundingBox( *this, nonLocalParticleCurrentBoundingBox );

            }

            virtual void setLocalCurrentSurfacePoints( ){

                asp::unit_test::aspBaseTester::set_localCurrentSurfacePoints( *this, localCurrentSurfacePoints );

            }

            virtual void setLocalReferenceSurfacePoints( ){

                asp::unit_test::aspBaseTester::set_localReferenceSurfacePoints( *this, localReferenceSurfacePoints );

            }

            virtual void setLocalDeformationGradient( ){

                asp::unit_test::aspBaseTester::set_localDeformationGradient( *this, localDeformationGradient );

            }

            virtual void setNonLocalReferenceRadius( ){

                asp::unit_test::aspBaseTester::set_nonLocalReferenceRadius( *this, nonLocalReferenceRadius );

            }

            virtual void setLocalMicroDeformation( ){

                asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, localMicroDeformation );

            }

            virtual void setLocalReferenceParticleSpacing( ){

                asp::unit_test::aspBaseTester::set_localReferenceParticleSpacing( *this, localReferenceParticleSpacing );

            }

            virtual void setLocalGradientMicroDeformation( ){

                asp::unit_test::aspBaseTester::set_localGradientMicroDeformation( *this, localGradientMicroDeformation );

            }

            virtual void setNonLocalMicroDeformationBase( ){

                asp::unit_test::aspBaseTester::set_nonLocalMicroDeformationBase( *this, nonLocalMicroDeformationBase );

            }

    };

    aspBaseMock asp;

    std::unordered_map< unsigned int, floatVector > answer = { { 0, { 0, 0, 0 } }, { 1, { -0.5, 0, 0 } }, { 2, { -0.65, 0, 0 } } };

    const std::unordered_map< unsigned int, floatVector > *result = asp.getParticlePairOverlap( );

    for ( auto p = answer.begin( ); p != answer.end( ); p++ ){

        auto search = result->find( p->first );

        BOOST_CHECK( search != result->end( ) );

        if ( search != result->end( ) ){

            BOOST_CHECK( vectorTools::fuzzyEquals( p->second, search->second ) );

        }

    }

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceOverlapEnergyDensity ){

    class aspBaseMock : public asp::aspBase{

        public:
   
            std::unordered_map< unsigned int, floatVector > particlePairOverlap = { { 0, { -0.5, 0, 0 } }, { 4, { 2, -1, 4 } } };

            floatVector surfaceOverlapParameters = { 2.3 };

        private:

            virtual void setParticlePairOverlap( ){

                asp::unit_test::aspBaseTester::set_particlePairOverlap( *this, particlePairOverlap );

            }

            virtual void setSurfaceOverlapParameters( ){

                asp::unit_test::aspBaseTester::set_surfaceOverlapParameters( *this, surfaceOverlapParameters );

            }

    };

    aspBaseMock asp, aspGet;

    std::unordered_map< unsigned int, floatType > result;

    BOOST_CHECK_NO_THROW( asp.computeSurfaceOverlapEnergyDensity( result ) );

    const std::unordered_map< unsigned int, floatType > *resultGet;

    resultGet = aspGet.getSurfaceOverlapEnergyDensity( );

    for ( auto p = asp.particlePairOverlap.begin( ); p != asp.particlePairOverlap.end( ); p++ ){

        auto overlapEnergy = result.find( p->first );

        BOOST_CHECK( overlapEnergy != result.end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapEnergy->second, 0.5 * asp.surfaceOverlapParameters[ 0 ] * vectorTools::dot( p->second, p->second ) ) );

        auto overlapEnergyGet = resultGet->find( p->first );

        BOOST_CHECK( overlapEnergyGet != result.end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapEnergyGet->second, 0.5 * asp.surfaceOverlapParameters[ 0 ] * vectorTools::dot( p->second, p->second ) ) );

    } 

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceOverlapTraction ){

    class aspBaseMock : public asp::aspBase{

        public:
   
            std::unordered_map< unsigned int, floatVector > particlePairOverlap = { { 0, { -0.5, 0, 0 } }, { 4, { 2, -1, 4 } } };

            floatVector surfaceOverlapParameters = { 2.3 };

        private:

            virtual void setParticlePairOverlap( ){

                asp::unit_test::aspBaseTester::set_particlePairOverlap( *this, particlePairOverlap );

            }

            virtual void setSurfaceOverlapParameters( ){

                asp::unit_test::aspBaseTester::set_surfaceOverlapParameters( *this, surfaceOverlapParameters );

            }

    };

    aspBaseMock asp, aspGet;

    std::unordered_map< unsigned int, floatVector > result;

    BOOST_CHECK_NO_THROW( asp.computeSurfaceOverlapTraction( result ) );

    const std::unordered_map< unsigned int, floatVector > *resultGet;

    resultGet = aspGet.getSurfaceOverlapTraction( );

    for ( auto p = asp.particlePairOverlap.begin( ); p != asp.particlePairOverlap.end( ); p++ ){

        auto overlapTraction = result.find( p->first );

        BOOST_CHECK( overlapTraction != result.end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapTraction->second, asp.surfaceOverlapParameters[ 0 ] * p->second ) );

        auto overlapTractionGet = resultGet->find( p->first );

        BOOST_CHECK( overlapTractionGet != result.end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapTraction->second, asp.surfaceOverlapParameters[ 0 ] * p->second ) );

    } 

}

BOOST_AUTO_TEST_CASE( test_aspBase_resetInteractionPairData ){

    class aspBaseMock : public asp::aspBase{

        public:

            asp::dataStorage< int > interactionData1;
        
            asp::dataStorage< floatVector > interactionData2;
    
            asp::dataStorage< floatType > surfaceData1;
    
            asp::dataStorage< floatMatrix > surfaceData2;
    
            asp::dataStorage< std::unordered_map< int, floatType > > particleData1;

            aspBaseMock( ){

                interactionData1 = asp::dataStorage< int >( true, 1 );
        
                interactionData2 = asp::dataStorage< floatVector >( true, { 1, 2, 3 } );
    
                surfaceData1 = asp::dataStorage< floatType >( true, 2. );
    
                surfaceData2 = asp::dataStorage< floatMatrix >( true, { { 1, 2, 3 }, { 4, 5, 6 } } );
    
                particleData1 = asp::dataStorage< std::unordered_map< int, floatType > >( true, { { 1, 1 }, { 2, 3 } } );

                addInteractionPairData( &interactionData1 );

                addInteractionPairData( &interactionData2 );

                addSurfacePointData( &surfaceData1 );

                addSurfacePointData( &surfaceData2 );

                addLocalParticleData( &particleData1 );

            }

    };

    aspBaseMock asp;

    BOOST_CHECK( asp::unit_test::aspBaseTester::getInteractionPairData( asp ).size( ) == 2 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getSurfacePointData( asp ).size( ) == 2 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getLocalParticleData( asp ).size( ) == 1 );

    asp::unit_test::aspBaseTester::resetInteractionPairData( asp );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getInteractionPairData( asp ).size( ) == 0 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getSurfacePointData( asp ).size( ) == 2 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getLocalParticleData( asp ).size( ) == 1 );

    BOOST_CHECK( !asp.interactionData1.first );

    BOOST_CHECK( !asp.interactionData2.first );

    BOOST_CHECK( asp.surfaceData1.first );

    BOOST_CHECK( asp.surfaceData2.first );

    BOOST_CHECK( asp.particleData1.first );

}

BOOST_AUTO_TEST_CASE( test_aspBase_resetSurfacePointData ){

    class aspBaseMock : public asp::aspBase{

        public:

            asp::dataStorage< int > interactionData1;
        
            asp::dataStorage< floatVector > interactionData2;
    
            asp::dataStorage< floatType > surfaceData1;
    
            asp::dataStorage< floatMatrix > surfaceData2;
    
            asp::dataStorage< std::unordered_map< int, floatType > > particleData1;

            aspBaseMock( ){

                interactionData1 = asp::dataStorage< int >( true, 1 );
        
                interactionData2 = asp::dataStorage< floatVector >( true, { 1, 2, 3 } );
    
                surfaceData1 = asp::dataStorage< floatType >( true, 2. );
    
                surfaceData2 = asp::dataStorage< floatMatrix >( true, { { 1, 2, 3 }, { 4, 5, 6 } } );
    
                particleData1 = asp::dataStorage< std::unordered_map< int, floatType > >( true, { { 1, 1 }, { 2, 3 } } );

                addInteractionPairData( &interactionData1 );

                addInteractionPairData( &interactionData2 );

                addSurfacePointData( &surfaceData1 );

                addSurfacePointData( &surfaceData2 );

                addLocalParticleData( &particleData1 );

            }

    };

    aspBaseMock asp;

    BOOST_CHECK( asp::unit_test::aspBaseTester::getInteractionPairData( asp ).size( ) == 2 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getSurfacePointData( asp ).size( ) == 2 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getLocalParticleData( asp ).size( ) == 1 );

    asp::unit_test::aspBaseTester::resetSurfacePointData( asp );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getInteractionPairData( asp ).size( ) == 0 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getSurfacePointData( asp ).size( ) == 0 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getLocalParticleData( asp ).size( ) == 1 );

    BOOST_CHECK( !asp.interactionData1.first );

    BOOST_CHECK( !asp.interactionData2.first );

    BOOST_CHECK( !asp.surfaceData1.first );

    BOOST_CHECK( !asp.surfaceData2.first );

    BOOST_CHECK( asp.particleData1.first );

}

BOOST_AUTO_TEST_CASE( test_aspBase_resetLocalParticleData ){

    class aspBaseMock : public asp::aspBase{

        public:

            asp::dataStorage< int > interactionData1;
        
            asp::dataStorage< floatVector > interactionData2;
    
            asp::dataStorage< floatType > surfaceData1;
    
            asp::dataStorage< floatMatrix > surfaceData2;
    
            asp::dataStorage< std::unordered_map< int, floatType > > particleData1;

            aspBaseMock( ){

                interactionData1 = asp::dataStorage< int >( true, 1 );
        
                interactionData2 = asp::dataStorage< floatVector >( true, { 1, 2, 3 } );
    
                surfaceData1 = asp::dataStorage< floatType >( true, 2. );
    
                surfaceData2 = asp::dataStorage< floatMatrix >( true, { { 1, 2, 3 }, { 4, 5, 6 } } );
    
                particleData1 = asp::dataStorage< std::unordered_map< int, floatType > >( true, { { 1, 1 }, { 2, 3 } } );

                addInteractionPairData( &interactionData1 );

                addInteractionPairData( &interactionData2 );

                addSurfacePointData( &surfaceData1 );

                addSurfacePointData( &surfaceData2 );

                addLocalParticleData( &particleData1 );

            }

    };

    aspBaseMock asp;

    BOOST_CHECK( asp::unit_test::aspBaseTester::getInteractionPairData( asp ).size( ) == 2 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getSurfacePointData( asp ).size( ) == 2 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getLocalParticleData( asp ).size( ) == 1 );

    asp::unit_test::aspBaseTester::resetLocalParticleData( asp );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getInteractionPairData( asp ).size( ) == 0 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getSurfacePointData( asp ).size( ) == 0 );

    BOOST_CHECK( asp::unit_test::aspBaseTester::getLocalParticleData( asp ).size( ) == 0 );

    BOOST_CHECK( !asp.interactionData1.first );

    BOOST_CHECK( !asp.interactionData2.first );

    BOOST_CHECK( !asp.surfaceData1.first );

    BOOST_CHECK( !asp.surfaceData2.first );

    BOOST_CHECK( !asp.particleData1.first );

}
