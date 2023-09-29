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

                static bool searchLocalParticleData( const asp::aspBase &asp, const asp::dataBase *value ){

                    return std::find( asp._localParticleData.begin( ), asp._localParticleData.end( ), value ) != asp._localParticleData.end( );

                }

                static bool searchSurfacePointData( const asp::aspBase &asp, const asp::dataBase *value ){

                    return std::find( asp._surfacePointData.begin( ), asp._surfacePointData.end( ), value ) != asp._surfacePointData.end( );

                }

                static bool searchInteractionPairData( const asp::aspBase &asp, const asp::dataBase *value ){

                    return std::find( asp._interactionPairData.begin( ), asp._interactionPairData.end( ), value ) != asp._interactionPairData.end( );

                }

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

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localReferenceRadius ) );

                    return;

                }

                static void setNonLocalReferenceRadius( asp::aspBase &asp, asp::dataStorage< floatType > &radius ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalReferenceRadius( ) );

                    radius = asp._nonLocalReferenceRadius;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._nonLocalReferenceRadius ) );

                    return;

                }

                static void setLocalReferenceNormal( asp::aspBase &asp, asp::dataStorage< floatVector > &normal ){

                    BOOST_CHECK_NO_THROW( asp.setLocalReferenceNormal( ) );

                    normal = asp._localReferenceNormal;

                    BOOST_CHECK( searchSurfacePointData( asp, &asp._localReferenceNormal ) );

                    return;

                }

                static void setLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp, asp::dataStorage< floatVector > &Xi ){

                    BOOST_CHECK_NO_THROW( asp.setLocalSurfaceReferenceRelativePositionVector( ) );

                    Xi = asp._localSurfaceReferenceRelativePositionVector;

                    BOOST_CHECK( searchSurfacePointData( asp, &asp._localSurfaceReferenceRelativePositionVector ) );

                    return;

                }

                static void setNonLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp, asp::dataStorage< floatVector > &Xi ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalSurfaceReferenceRelativePositionVector( ) );

                    Xi = asp._nonLocalSurfaceReferenceRelativePositionVector;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._nonLocalSurfaceReferenceRelativePositionVector ) );

                    return;

                }

                static void setLocalDeformationGradient( asp::aspBase &asp, asp::dataStorage< floatVector > &localDeformationGradient ){

                    BOOST_CHECK_NO_THROW( asp.setLocalDeformationGradient( ) );

                    localDeformationGradient = asp._localDeformationGradient;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localDeformationGradient ) );

                    return;

                }

                static void setPreviousLocalDeformationGradient( asp::aspBase &asp, asp::dataStorage< floatVector > &previousLocalDeformationGradient ){

                    BOOST_CHECK_NO_THROW( asp.setPreviousLocalDeformationGradient( ) );

                    previousLocalDeformationGradient = asp._previousLocalDeformationGradient;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._previousLocalDeformationGradient ) );

                    return;

                }

                static void setLocalMicroDeformation( asp::aspBase &asp, asp::dataStorage< floatVector > &localMicroDeformation ){

                    BOOST_CHECK_NO_THROW( asp.setLocalMicroDeformation( ) );

                    localMicroDeformation = asp._localMicroDeformation;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localMicroDeformation ) );

                    return;

                }

                static void setPreviousLocalMicroDeformation( asp::aspBase &asp, asp::dataStorage< floatVector > &previousLocalMicroDeformation ){

                    BOOST_CHECK_NO_THROW( asp.setPreviousLocalMicroDeformation( ) );

                    previousLocalMicroDeformation = asp._previousLocalMicroDeformation;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._previousLocalMicroDeformation ) );

                    return;

                }

                static void setReferenceDistanceVector( asp::aspBase &asp, asp::dataStorage< floatVector > &referenceDistanceVector ){

                    BOOST_CHECK_NO_THROW( asp.setReferenceDistanceVector( ) );

                    referenceDistanceVector = asp._referenceDistanceVector;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._referenceDistanceVector ) );

                    return;

                }

                static void setLocalReferenceParticleSpacingVector( asp::aspBase &asp,
                                                                  asp::dataStorage< floatVector > &localReferenceParticleSpacing ){

                    BOOST_CHECK_NO_THROW( asp.setLocalReferenceParticleSpacingVector( ) );

                    localReferenceParticleSpacing = asp._localReferenceParticleSpacing;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._localReferenceParticleSpacing ) );

                    return;

                }

                static void setNonLocalMicroDeformation( asp::aspBase &asp,
                                                             asp::dataStorage< floatVector > &nonLocalMicroDeformation ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalMicroDeformation( ) );

                    nonLocalMicroDeformation = asp._nonLocalMicroDeformation;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._nonLocalMicroDeformation ) );

                    return;

                }

                static void setLocalCurrentNormal( asp::aspBase &asp,
                                                       asp::dataStorage< floatVector > &localCurrentNormal ){

                    BOOST_CHECK_NO_THROW( asp.setLocalCurrentNormal( ) );

                    localCurrentNormal = asp._localCurrentNormal;

                    BOOST_CHECK( searchSurfacePointData( asp, &asp._localCurrentNormal ) );

                    return;

                }

                static void setCurrentDistanceVector( asp::aspBase &asp,
                                                    asp::dataStorage< floatVector > &currentDistance ){

                    BOOST_CHECK_NO_THROW( asp.setCurrentDistanceVector( ) );

                    currentDistance = asp._currentDistanceVector;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._currentDistanceVector ) );

                    return;

                }

                static void setSurfaceParameters( asp::aspBase &asp,
                                                      asp::dataStorage< floatVector > &currentSurfaceParameters ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceParameters( ) );

                    currentSurfaceParameters = asp._surfaceParameters;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._surfaceParameters ) );

                    return;

                }

                static void setSurfaceAdhesionEnergyDensity( asp::aspBase &asp,
                                                             asp::dataStorage< floatType > &surfaceAdhesionEnergyDensity ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceAdhesionEnergyDensity( ) );

                    surfaceAdhesionEnergyDensity = asp._surfaceAdhesionEnergyDensity;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._surfaceAdhesionEnergyDensity ) );

                    return;

                }

                static void setSurfaceAdhesionTraction( asp::aspBase &asp,
                                                        asp::dataStorage< floatVector > &surfaceAdhesionTraction ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceAdhesionTraction( ) );

                    surfaceAdhesionTraction = asp._surfaceAdhesionTraction;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._surfaceAdhesionTraction ) );

                    return;

                }

                static void setLocalReferenceSurfacePoints( asp::aspBase &asp,
                                                            asp::dataStorage< floatVector > &localReferenceSurfacePoints ){

                    BOOST_CHECK_NO_THROW( asp.setLocalReferenceSurfacePoints( ) );

                    localReferenceSurfacePoints = asp._localReferenceSurfacePoints;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localReferenceSurfacePoints ) );

                    return;

                }

                static void setNonLocalReferenceSurfacePoints( asp::aspBase &asp,
                                                               asp::dataStorage< floatVector > &nonLocalReferenceSurfacePoints ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalReferenceSurfacePoints( ) );

                    nonLocalReferenceSurfacePoints = asp._nonLocalReferenceSurfacePoints;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._nonLocalReferenceSurfacePoints ) );

                    return;

                }

                static void setLocalCurrentSurfacePoints( asp::aspBase &asp,
                                                          asp::dataStorage< floatVector > &localCurrentSurfacePoints ){

                    BOOST_CHECK_NO_THROW( asp.setLocalCurrentSurfacePoints( ) );

                    localCurrentSurfacePoints = asp._localCurrentSurfacePoints;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localCurrentSurfacePoints ) );

                    return;

                }

                static void setNonLocalCurrentSurfacePoints( asp::aspBase &asp,
                                                             asp::dataStorage< floatVector > &nonLocalCurrentSurfacePoints ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalCurrentSurfacePoints( ) );

                    nonLocalCurrentSurfacePoints = asp._nonLocalCurrentSurfacePoints;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._nonLocalCurrentSurfacePoints ) );

                    return;

                }

                static void setLocalParticleCurrentBoundingBox( asp::aspBase &asp,
                                                                asp::dataStorage< floatMatrix > &localParticleCurrentBoundingBox ){

                    BOOST_CHECK_NO_THROW( asp.setLocalParticleCurrentBoundingBox( ) );

                    localParticleCurrentBoundingBox = asp._localParticleCurrentBoundingBox;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleCurrentBoundingBox ) );

                    return;

                }

                static void setNonLocalParticleCurrentBoundingBox( asp::aspBase &asp,
                                                                   asp::dataStorage< floatMatrix > &nonLocalParticleCurrentBoundingBox ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalParticleCurrentBoundingBox( ) );

                    nonLocalParticleCurrentBoundingBox = asp._nonLocalParticleCurrentBoundingBox;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._nonLocalParticleCurrentBoundingBox ) );

                    return;

                }

                static void setLocalGradientMicroDeformation( asp::aspBase &asp,
                                                              asp::dataStorage< floatVector > &localGradientMicroDeformation ){

                    BOOST_CHECK_NO_THROW( asp.setLocalGradientMicroDeformation( ) );

                    localGradientMicroDeformation = asp._localGradientMicroDeformation;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localGradientMicroDeformation ) );

                    return;

                }

                static void setNonLocalMicroDeformationBase( asp::aspBase &asp,
                                                             asp::dataStorage< floatVector > &nonLocalMicroDeformationBase ){

                    BOOST_CHECK_NO_THROW( asp.setNonLocalMicroDeformationBase( ) );

                    nonLocalMicroDeformationBase = asp._nonLocalMicroDeformationBase;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._nonLocalMicroDeformationBase ) );

                    return;

                }

                static void setSurfaceOverlapParameters( asp::aspBase &asp,
                                                         asp::dataStorage< floatVector > &surfaceOverlapParameters ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceOverlapParameters( ) );

                    surfaceOverlapParameters = asp._surfaceOverlapParameters;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._surfaceOverlapParameters ) );

                    return;

                }

                static void setParticlePairOverlap( asp::aspBase &asp,
                                                    asp::dataStorage< std::unordered_map< unsigned int, floatVector > > &particlePairOverlap ){

                    BOOST_CHECK_NO_THROW( asp.setParticlePairOverlap( ) );

                    particlePairOverlap = asp._particlePairOverlap;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._particlePairOverlap ) );

                    return;

                }

                static void setSurfaceOverlapEnergyDensity( asp::aspBase &asp,
                                                            asp::dataStorage< std::unordered_map< unsigned int, floatType > > &surfaceOverlapEnergyDensity ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceOverlapEnergyDensity( ) );

                    surfaceOverlapEnergyDensity = asp._surfaceOverlapEnergyDensity;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._surfaceOverlapEnergyDensity ) );

                    return;

                }

                static void setLocalParticleReferenceVolume( asp::aspBase &asp,
                                                             asp::dataStorage< floatType > &result ){

                    BOOST_CHECK_NO_THROW( asp.setLocalParticleReferenceVolume( ) );

                    result = asp._localParticleReferenceVolume;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleReferenceVolume ) );

                }

                static void setLocalParticleCurrentVolume( asp::aspBase &asp,
                                                           asp::dataStorage< floatType > &result ){

                    BOOST_CHECK_NO_THROW( asp.setLocalParticleCurrentVolume( ) );

                    result = asp._localParticleCurrentVolume;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleCurrentVolume ) );

                }

                static void setLocalParticleEnergy( asp::aspBase &asp,
                                                    asp::dataStorage< floatType > &result ){

                    try{
                        asp.setLocalParticleEnergy( );
                    }
                    catch(std::exception &e){
                        errorTools::printNestedExceptions( e );
                    }

                    BOOST_CHECK_NO_THROW( asp.setLocalParticleEnergy( ) );

                    result = asp._localParticleEnergy;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleEnergy ) );

                }

                static void setLocalParticleQuantities( asp::aspBase &asp,
                                                        asp::dataStorage< floatType > &energyResult,
                                                        asp::dataStorage< floatVector > &stressResult,
                                                        asp::dataStorage< floatVector > &sdvResult,
                                                        asp::dataStorage< floatType > &probabilityResult ){

                    BOOST_CHECK_NO_THROW( asp.setLocalParticleQuantities( ) );

                    energyResult = asp._localParticleEnergyDensity;

                    stressResult = asp._localParticleMicroCauchyStress;

                    sdvResult = asp._localParticleStateVariables;

                    probabilityResult = asp._localParticleLogProbabilityRatio;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleEnergyDensity       ) );

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleMicroCauchyStress   ) );

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleStateVariables      ) );

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleLogProbabilityRatio ) );

                }

                static void setLocalParticleParameters( asp::aspBase &asp,
                                                        asp::dataStorage< floatVector > &result ){

                    BOOST_CHECK_NO_THROW( asp.setLocalParticleParameters( ) );

                    result = asp._localParticleParameters;

                    BOOST_CHECK( searchLocalParticleData( asp, &asp._localParticleParameters ) );

                }

                static void setSurfaceAdhesionThickness( asp::aspBase &asp,
                                                         asp::dataStorage< floatType > &result ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceAdhesionThickness( ) );

                    result = asp._surfaceAdhesionThickness;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._surfaceAdhesionThickness ) );

                }

                static void setSurfaceOverlapThickness( asp::aspBase &asp,
                                                        asp::dataStorage< std::unordered_map< unsigned int, floatType > > &result ){

                    BOOST_CHECK_NO_THROW( asp.setSurfaceOverlapThickness( ) );

                    result = asp._surfaceOverlapThickness;

                    BOOST_CHECK( searchInteractionPairData( asp, &asp._surfaceOverlapThickness ) );

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

                    asp.addSurfacePointData( &asp._localReferenceNormal );

                    return;

                }

                static void set_localReferenceRadius( asp::aspBase &asp,
                                                          const floatType &radius ){

                    asp._localReferenceRadius.first = true;
                    asp._localReferenceRadius.second = radius;

                    asp.addLocalParticleData( &asp._localReferenceRadius );

                    return;

                }

                static void set_nonLocalReferenceRadius( asp::aspBase &asp,
                                                             const floatType &radius ){

                    asp._nonLocalReferenceRadius.first = true;
                    asp._nonLocalReferenceRadius.second = radius;

                    asp.addInteractionPairData( &asp._nonLocalReferenceRadius );

                    return;

                }

                static void set_deformationGradient( asp::aspBase &asp,
                                                         const floatVector &deformationGradient ){

                    asp._deformationGradient = deformationGradient;

                    return;

                }

                static void set_previousDeformationGradient( asp::aspBase &asp,
                                                             const floatVector &previousDeformationGradient ){

                    asp._previousDeformationGradient = previousDeformationGradient;

                    return;

                }

                static void set_microDeformation( asp::aspBase &asp,
                                                      const floatVector &microDeformation ){

                    asp._microDeformation = microDeformation;

                    return;

                }

                static void set_previousMicroDeformation( asp::aspBase &asp,
                                                          const floatVector &previousMicroDeformation ){

                    asp._previousMicroDeformation = previousMicroDeformation;

                    return;

                }

                static void set_localDeformationGradient( asp::aspBase &asp,
                                                              const floatVector &deformationGradient ){

                    asp._localDeformationGradient.first = true;
                    asp._localDeformationGradient.second = deformationGradient;

                    asp.addLocalParticleData( &asp._localDeformationGradient );

                    return;

                }

                static void set_localMicroDeformation( asp::aspBase &asp,
                                                           const floatVector &microDeformation ){

                    asp._localMicroDeformation.first = true;
                    asp._localMicroDeformation.second = microDeformation;

                    asp.addLocalParticleData( &asp._localMicroDeformation );

                    return;

                }

                static void set_nonLocalMicroDeformation( asp::aspBase &asp,
                                                              const floatVector &nonLocalMicroDeformation ){

                    asp._nonLocalMicroDeformation.first = true;
                    asp._nonLocalMicroDeformation.second = nonLocalMicroDeformation;

                    asp.addInteractionPairData( &asp._nonLocalMicroDeformation );

                    return;

                }

                static void set_localSurfaceReferenceRelativePositionVector( asp::aspBase &asp,
                                                                                 const floatVector &localXi ){

                    asp._localSurfaceReferenceRelativePositionVector.first = true;
                    asp._localSurfaceReferenceRelativePositionVector.second = localXi;

                    asp.addSurfacePointData( &asp._localSurfaceReferenceRelativePositionVector );

                    return;

                }

                static void set_nonLocalSurfaceReferenceRelativePositionVector( asp::aspBase &asp,
                                                                                    const floatVector &nonLocalXi ){

                    asp._nonLocalSurfaceReferenceRelativePositionVector.first = true;
                    asp._nonLocalSurfaceReferenceRelativePositionVector.second = nonLocalXi;

                    asp.addInteractionPairData( &asp._nonLocalSurfaceReferenceRelativePositionVector );

                    return;

                }

                static void set_referenceDistanceVector( asp::aspBase &asp, const floatVector &referenceDistanceVector ){

                    asp._referenceDistanceVector.first = true;
                    asp._referenceDistanceVector.second = referenceDistanceVector;

                    asp.addInteractionPairData( &asp._referenceDistanceVector );

                    return;

                }

                static void set_gradientMicroDeformation( asp::aspBase &asp, const floatVector &gradientMicroDeformation ){

                    asp._gradientMicroDeformation = gradientMicroDeformation;

                    return;

                }

                static void set_localReferenceParticleSpacing( asp::aspBase &asp, const floatVector &referenceParticleSpacing ){

                    asp._localReferenceParticleSpacing.first = true;
                    asp._localReferenceParticleSpacing.second = referenceParticleSpacing;

                    asp.addInteractionPairData( &asp._localReferenceParticleSpacing );

                    return;

                }

                static void set_surfaceParameters( asp::aspBase &asp, const floatVector &surfaceParameters ){

                    asp._surfaceParameters.first = true;
                    asp._surfaceParameters.second = surfaceParameters;

                    asp.addInteractionPairData( &asp._surfaceParameters );

                    return;

                }

                static void set_currentDistanceVector( asp::aspBase &asp, const floatVector &currentDistance ){

                    asp._currentDistanceVector.first = true;
                    asp._currentDistanceVector.second = currentDistance;

                    asp.addInteractionPairData( &asp._currentDistanceVector );

                    return;

                }

                static void set_localCurrentNormal( asp::aspBase &asp, const floatVector &localCurrentNormal ){

                    asp._localCurrentNormal.first = true;
                    asp._localCurrentNormal.second = localCurrentNormal;

                    asp.addSurfacePointData( &asp._localCurrentNormal );

                    return;

                }

                static void set_localReferenceSurfacePoints( asp::aspBase &asp, const floatVector &localReferenceSurfacePoints ){

                    asp._localReferenceSurfacePoints.first = true;
                    asp._localReferenceSurfacePoints.second = localReferenceSurfacePoints;

                    asp.addLocalParticleData( &asp._localReferenceSurfacePoints );

                    return;

                }

                static void set_nonLocalReferenceSurfacePoints( asp::aspBase &asp, const floatVector &nonLocalReferenceSurfacePoints ){

                    asp._nonLocalReferenceSurfacePoints.first = true;
                    asp._nonLocalReferenceSurfacePoints.second = nonLocalReferenceSurfacePoints;

                    asp.addInteractionPairData( &asp._nonLocalReferenceSurfacePoints );

                    return;

                }

                static void set_localCurrentSurfacePoints( asp::aspBase &asp, const floatVector &localCurrentSurfacePoints ){

                    asp._localCurrentSurfacePoints.first = true;
                    asp._localCurrentSurfacePoints.second = localCurrentSurfacePoints;

                    asp.addLocalParticleData( &asp._localCurrentSurfacePoints );

                    return;

                }

                static void set_nonLocalCurrentSurfacePoints( asp::aspBase &asp, const floatVector &nonLocalCurrentSurfacePoints ){

                    asp._nonLocalCurrentSurfacePoints.first = true;
                    asp._nonLocalCurrentSurfacePoints.second = nonLocalCurrentSurfacePoints;

                    asp.addInteractionPairData( &asp._nonLocalCurrentSurfacePoints );

                    return;

                }

                static void set_surfaceOverlapParameters( asp::aspBase &asp, const floatVector &nonLocalCurrentSurfacePoints ){

                    asp._surfaceOverlapParameters.first = true;
                    asp._surfaceOverlapParameters.second = nonLocalCurrentSurfacePoints;

                    asp.addInteractionPairData( &asp._surfaceOverlapParameters );

                    return;

                }

                static void set_nonLocalParticleCurrentBoundingBox( asp::aspBase &asp, const floatMatrix &nonLocalParticleBoundingBox ){

                    asp._nonLocalParticleCurrentBoundingBox.first = true;
                    asp._nonLocalParticleCurrentBoundingBox.second = nonLocalParticleBoundingBox;

                    asp.addInteractionPairData( &asp._nonLocalParticleCurrentBoundingBox );

                    return;

                }

                static void set_localGradientMicroDeformation( asp::aspBase &asp, const floatVector &localGradientMicroDeformation ){

                    asp._localGradientMicroDeformation.first = true;
                    asp._localGradientMicroDeformation.second = localGradientMicroDeformation;

                    asp.addLocalParticleData( &asp._localGradientMicroDeformation );

                    return;

                }

                static void set_nonLocalMicroDeformationBase( asp::aspBase &asp, const floatVector &nonLocalMicroDeformationBase ){

                    asp._nonLocalMicroDeformationBase.first = true;
                    asp._nonLocalMicroDeformationBase.second = nonLocalMicroDeformationBase;

                    asp.addInteractionPairData( &asp._nonLocalMicroDeformationBase );

                    return;

                }

                static void set_particlePairOverlap( asp::aspBase &asp, const std::unordered_map< unsigned int, floatVector > &particlePairOverlap ){

                    asp._particlePairOverlap.first = true;
                    asp._particlePairOverlap.second = particlePairOverlap;

                    asp.addInteractionPairData( &asp._particlePairOverlap );

                    return;

                }

                static void set_localParticleEnergyDensity( asp::aspBase &asp, floatType &value ){

                    asp._localParticleEnergyDensity.first = true;
                    asp._localParticleEnergyDensity.second = value;

                    asp.addLocalParticleData( &asp._localParticleEnergyDensity );

                }

                static void set_localParticleEnergy( asp::aspBase &asp, floatType &value ){

                    asp._localParticleEnergy.first = true;
                    asp._localParticleEnergy.second = value;

                    asp.addLocalParticleData( &asp._localParticleEnergy );

                }

                static void set_localParticleReferenceVolume( asp::aspBase &asp, floatType &referenceVolume ){

                    asp._localParticleReferenceVolume.first = true;
                    asp._localParticleReferenceVolume.second = referenceVolume;

                    asp.addLocalParticleData( &asp._localParticleReferenceVolume );

                }

                static void set_localParticleCurrentVolume( asp::aspBase &asp, floatType &currentVolume ){

                    asp._localParticleCurrentVolume.first = true;
                    asp._localParticleCurrentVolume.second = currentVolume;

                    asp.addLocalParticleData( &asp._localParticleCurrentVolume );

                }

                static void set_particleParameters( asp::aspBase &asp, floatVector &value ){

                    asp._particleParameters = value;

                }

                static void set_localParticleParameters( asp::aspBase &asp, floatVector &value ){

                    asp._localParticleParameters.first = true;
                    asp._localParticleParameters.second = value;

                    asp.addLocalParticleData( &asp._localParticleParameters );

                }

                static void set_localParticleMicroCauchyStress( asp::aspBase &asp, floatVector &value ){

                    asp._localParticleMicroCauchyStress.first = true;
                    asp._localParticleMicroCauchyStress.second = value;

                    asp.addLocalParticleData( &asp._localParticleMicroCauchyStress );

                }

                static void set_localParticleStateVariables( asp::aspBase &asp, floatVector &value ){

                    asp._localParticleStateVariables.first = true;
                    asp._localParticleStateVariables.second = value;

                    asp.addLocalParticleData( &asp._localParticleStateVariables );

                }

                static void set_localParticleLogProbabilityRatio( asp::aspBase &asp, floatType &value ){

                    asp._localParticleLogProbabilityRatio.first = true;
                    asp._localParticleLogProbabilityRatio.second = value;

                    asp.addLocalParticleData( &asp._localParticleLogProbabilityRatio );

                }

                static void set_numLocalParticles( asp::aspBase &asp, unsigned int &value ){

                    asp._numLocalParticles = value;

                }

                static void set_surfaceAdhesionEnergyDensity( asp::aspBase &asp, floatType &value ){

                    asp._surfaceAdhesionEnergyDensity.first = true;
                    asp._surfaceAdhesionEnergyDensity.second = value;

                    asp.addInteractionPairData( &asp._surfaceAdhesionEnergyDensity );

                }

                static void set_surfaceAdhesionTraction( asp::aspBase &asp, floatVector &value ){

                    asp._surfaceAdhesionTraction.first = true;
                    asp._surfaceAdhesionTraction.second = value;

                    asp.addInteractionPairData( &asp._surfaceAdhesionTraction );

                }

                static void set_surfaceAdhesionThickness( asp::aspBase &asp, floatType &value ){

                    asp._surfaceAdhesionThickness.first = true;
                    asp._surfaceAdhesionThickness.second = value;

                    asp.addInteractionPairData( &asp._surfaceAdhesionThickness );

                }

                static void set_surfaceOverlapEnergyDensity( asp::aspBase &asp, mapFloatType &value ){

                    asp._surfaceOverlapEnergyDensity.first = true;
                    asp._surfaceOverlapEnergyDensity.second = value;

                    asp.addInteractionPairData( &asp._surfaceOverlapEnergyDensity );

                }

                static void set_surfaceOverlapTraction( asp::aspBase &asp, mapFloatVector &value ){

                    asp._surfaceOverlapTraction.first = true;
                    asp._surfaceOverlapTraction.second = value;

                    asp.addInteractionPairData( &asp._surfaceOverlapTraction );

                }

                static void set_surfaceOverlapThickness( asp::aspBase &asp, mapFloatType &value ){

                    asp._surfaceOverlapThickness.first = true;
                    asp._surfaceOverlapThickness.second = value;

                    asp.addInteractionPairData( &asp._surfaceOverlapThickness );

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

                static void resetInteractionPairData( asp::aspBase &asp ){

                    asp.resetInteractionPairData( );

                }

                static void resetSurfacePointData( asp::aspBase &asp ){

                    asp.resetSurfacePointData( );

                }

                static void resetLocalParticleData( asp::aspBase &asp ){

                    asp.resetLocalParticleData( );

                }

                static void checkNumLocalParticles( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._numLocalParticles == asp.getNumLocalParticles( ) );

                }

                static void checkPreviousTime( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._previousTime == asp.getPreviousTime( ) );

                }

                static void checkDeltaTime( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._deltaTime == asp.getDeltaTime( ) );

                }

                static void checkPreviousTemperature( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._previousTemperature == asp.getPreviousTemperature( ) );

                }

                static void checkTemperature( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._temperature == asp.getTemperature( ) );

                }

                static void checkPreviousDeformationGradient( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._previousDeformationGradient == asp.getPreviousDeformationGradient( ) );

                }

                static void checkPreviousMicroDeformation( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._previousMicroDeformation == asp.getPreviousMicroDeformation( ) );

                }

                static void checkPreviousStateVariables( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._previousStateVariables == asp.getPreviousStateVariables( ) );

                }

                static void checkPreviousLocalStateVariables( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._previousStateVariables == asp.getPreviousLocalStateVariables( ) );

                }

                static void checkDeformationGradient( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._deformationGradient == asp.getDeformationGradient( ) );

                }

                static void checkParticleParameters( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._particleParameters == asp.getParticleParameters( ) );

                }

                static void checkLocalIndex( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._localIndex == asp.getLocalIndex( ) );

                }

                static void checkNonLocalIndex( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._nonLocalIndex == asp.getNonLocalIndex( ) );

                }

                static void checkLocalSurfaceNodeIndex( asp::aspBase &asp ){

                    BOOST_CHECK( &asp._localSurfaceNodeIndex == asp.getLocalSurfaceNodeIndex( ) );

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

    class aspBaseMock : public asp::aspBase{

        public:

            aspBaseMock( floatType &radius ){

                asp::unit_test::aspBaseTester::set_radius( *this, radius );

            }

    };

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

    floatVector currentStateVariables;

    floatVector parameters = { 30., 20. };

    floatType radius = 2.4;

    floatType energy = 0;

    floatVector cauchyStress;

    aspBaseMock asp( radius );

    BOOST_CHECK_NO_THROW( asp.computeLocalParticleEnergyDensity( previousTime, deltaTime, currentDeformationGradient, previousDeformationGradient,
                                                                 currentTemperature, previousTemperature, previousStateVariables, parameters, energy, cauchyStress, currentStateVariables ) );

    BOOST_CHECK( !vectorTools::fuzzyEquals( energy, 0. ) );

    BOOST_CHECK( cauchyStress.size( ) == currentDeformationGradient.size( ) );

    energy = 0;

    cauchyStress.clear( );

    floatType logProbabilityRatio = 3.4;

    BOOST_CHECK_NO_THROW( asp.computeLocalParticleEnergyDensity( previousTime, deltaTime, currentDeformationGradient, previousDeformationGradient,
                                                                 currentTemperature, previousTemperature, previousStateVariables, parameters, energy, cauchyStress, currentStateVariables, logProbabilityRatio ) );

    BOOST_CHECK( !vectorTools::fuzzyEquals( energy, 0. ) );

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

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalReferenceNormal ){

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

    floatVector result;

    for ( unsigned int i = 0; i < normalAnswers.size( ); i++ ){

        asp.getLocalReferenceNormal( i, result );

        BOOST_CHECK( vectorTools::fuzzyEquals( result, normalAnswers[ i ] ) );

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

BOOST_AUTO_TEST_CASE( test_aspBase_setPreviousLocalDeformationGradient ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        
    asp::unit_test::aspBaseTester::set_previousDeformationGradient( asp, answer );

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setPreviousLocalDeformationGradient( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    asp::unit_test::aspBaseTester::set_previousDeformationGradient( aspGet, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getPreviousLocalDeformationGradient( ), answer ) );

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

BOOST_AUTO_TEST_CASE( test_aspBase_setPreviousLocalMicroDeformation ){

    class aspBaseMock : public asp::aspBase{

    };

    aspBaseMock asp;

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        
    asp::unit_test::aspBaseTester::set_previousMicroDeformation( asp, answer );

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setPreviousLocalMicroDeformation( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    asp::unit_test::aspBaseTester::set_previousMicroDeformation( aspGet, answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getPreviousLocalMicroDeformation( ), answer ) );

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

    asp::unit_test::aspBaseTester::setLocalReferenceParticleSpacingVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalReferenceParticleSpacingVector( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalMicroDeformation ){

    class aspBaseMock : public asp::aspBase{

        floatVector microDeformationBase = { 4,  5,  6,
                                             7,  8,  9,
                                            10, 11, 12 };

        void setLocalReferenceParticleSpacingVector( ){

            floatVector value = { 1, 2, 3 };

            asp::unit_test::aspBaseTester::set_localReferenceParticleSpacing( *this, value );

            return;

        }

        void setNonLocalMicroDeformationBase( ){

            asp::aspBase::setNonLocalMicroDeformationBase( microDeformationBase );

        }

    };

    aspBaseMock asp;

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

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspGet, gradientMicroDeformation );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getNonLocalMicroDeformation( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalMicroDeformationDerivatives ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector microDeformationBase = { 4,  5,  6,
                                                 7,  8,  9,
                                                10, 11, 12 };
    
            floatVector Xi = { 1, 2, 3 };
    
            floatVector D = { 3, 4, 5 };
    
            floatVector XiNL = { 6, 7, 8 };

        private:

            void setLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setLocalSurfaceReferenceRelativePositionVector( Xi );
    
                return;
    
            }
    
            void setReferenceDistanceVector( ) override {
    
                asp::aspBase::setReferenceDistanceVector( D );
    
                return;
    
            }
    
            void setNonLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setNonLocalSurfaceReferenceRelativePositionVector( XiNL );
    
                return;
    
            }
    
            void setNonLocalMicroDeformationBase( ) override{
    
                asp::aspBase::setNonLocalMicroDeformationBase( microDeformationBase );
    
            }

    };

    aspBaseMock asp1, asp2, asp3, asp4, asp5;

    floatVector gradientMicroDeformation = { 13, 14, 15, 16, 17, 18, 19, 20, 21,
                                             22, 23, 24, 25, 26, 27, 28, 29, 30,
                                             31, 32, 33, 34, 35, 36, 37, 38, 39 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp1, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp2, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp3, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp4, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp5, gradientMicroDeformation );

    floatMatrix dChiNLdChiNLBase( 9, floatVector( 9, 0 ) );

    floatMatrix dChiNLdGradChi( 9, floatVector( 27, 0 ) );

    floatMatrix dChiNLdXi( 9, floatVector( 3, 0 ) );

    floatMatrix dChiNLdD( 9, floatVector( 3, 0 ) );

    floatMatrix dChiNLdXiNL( 9, floatVector( 3, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < asp1.microDeformationBase.size( ); i++ ){

        floatVector deltas( asp1.microDeformationBase.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp1.microDeformationBase[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.microDeformationBase += deltas;

        aspm.microDeformationBase -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 9; j++ ){

            dChiNLdChiNLBase[ j ][ i ] = ( ( *aspp.getNonLocalMicroDeformation( ) )[ j ] - ( *aspm.getNonLocalMicroDeformation( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector deltas( gradientMicroDeformation.size( ), 0 );

        deltas[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation + deltas );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation - deltas );

        for ( unsigned int j = 0; j < 9; j++ ){

            dChiNLdGradChi[ j ][ i ] = ( ( *aspp.getNonLocalMicroDeformation( ) )[ j ] - ( *aspm.getNonLocalMicroDeformation( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp3.Xi.size( ); i++ ){

        floatVector deltas( asp3.Xi.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp3.Xi[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.Xi += deltas;

        aspm.Xi -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 9; j++ ){

            dChiNLdXi[ j ][ i ] = ( ( *aspp.getNonLocalMicroDeformation( ) )[ j ] - ( *aspm.getNonLocalMicroDeformation( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp4.D.size( ); i++ ){

        floatVector deltas( asp4.D.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp4.D[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.D += deltas;

        aspm.D -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 9; j++ ){

            dChiNLdD[ j ][ i ] = ( ( *aspp.getNonLocalMicroDeformation( ) )[ j ] - ( *aspm.getNonLocalMicroDeformation( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp5.XiNL.size( ); i++ ){

        floatVector deltas( asp5.XiNL.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp5.XiNL[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.XiNL += deltas;

        aspm.XiNL -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 9; j++ ){

            dChiNLdXiNL[ j ][ i ] = ( ( *aspp.getNonLocalMicroDeformation( ) )[ j ] - ( *aspm.getNonLocalMicroDeformation( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dChiNLdChiNLBase, *asp1.getdNonLocalMicroDeformationdNonLocalMicroDeformationBase( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dChiNLdGradChi, *asp2.getdNonLocalMicroDeformationdGradientMicroDeformation( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dChiNLdXi, *asp3.getdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dChiNLdD, *asp4.getdNonLocalMicroDeformationdLocalReferenceDistanceVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dChiNLdXiNL, *asp5.getdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( ) ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalMicroDeformationSecondDerivatives ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector microDeformationBase = { 4,  5,  6,
                                                 7,  8,  9,
                                                10, 11, 12 };
    
            floatVector Xi = { 1, 2, 3 };
    
            floatVector D = { 3, 4, 5 };
    
            floatVector XiNL = { 6, 7, 8 };

        private:

            void setLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setLocalSurfaceReferenceRelativePositionVector( Xi );
    
                return;
    
            }
    
            void setReferenceDistanceVector( ) override {
    
                asp::aspBase::setReferenceDistanceVector( D );
    
                return;
    
            }
    
            void setNonLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setNonLocalSurfaceReferenceRelativePositionVector( XiNL );
    
                return;
    
            }
    
            void setNonLocalMicroDeformationBase( ) override{
    
                asp::aspBase::setNonLocalMicroDeformationBase( microDeformationBase );
    
            }

    };

    aspBaseMock asp1, asp2, asp3, asp4, asp5;

    floatVector gradientMicroDeformation = { 13, 14, 15, 16, 17, 18, 19, 20, 21,
                                             22, 23, 24, 25, 26, 27, 28, 29, 30,
                                             31, 32, 33, 34, 35, 36, 37, 38, 39 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp1, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp2, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp3, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp4, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp5, gradientMicroDeformation );

    floatMatrix d2ChiNLdXidGradChi( 9, floatVector( 3 * 27, 0 ) );

    floatMatrix d2ChiNLdDdGradChi( 9, floatVector( 3 * 27, 0 ) );

    floatMatrix d2ChiNLdXiNLdGradChi( 9, floatVector( 3 * 27, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector deltas( gradientMicroDeformation.size( ), 0 );

        deltas[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation + deltas );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation - deltas );

        for ( unsigned int j = 0; j < 9; j++ ){

            for ( unsigned int k = 0; k < 3; k++ ){

                d2ChiNLdXidGradChi[ j ][ 27 * k + i ] = ( ( *aspp.getdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ) )[ j ][ k ] - ( *aspm.getdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

                d2ChiNLdDdGradChi[ j ][ 27 * k + i ] = ( ( *aspp.getdNonLocalMicroDeformationdLocalReferenceDistanceVector( ) )[ j ][ k ] - ( *aspm.getdNonLocalMicroDeformationdLocalReferenceDistanceVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

                d2ChiNLdXiNLdGradChi[ j ][ 27 * k + i ] = ( ( *aspp.getdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( ) )[ j ][ k ] - ( *aspm.getdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ChiNLdXidGradChi, *asp1.getd2NonLocalMicroDeformationdLocalReferenceRelativePositionVectordGradientMicroDeformation( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ChiNLdDdGradChi, *asp2.getd2NonLocalMicroDeformationdLocalReferenceDistanceVectordGradientMicroDeformation( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ChiNLdXiNLdGradChi, *asp3.getd2NonLocalMicroDeformationdNonLocalReferenceRelativePositionVectordGradientMicroDeformation( ) ) );

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

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalCurrentNormalGradients ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector N = { 1, 2, 3 };
    
            floatVector chi = { 0.39293837, -0.42772133, -0.54629709,
                                0.10262954,  0.43893794, -0.15378708,
                                0.9615284 ,  0.36965948, -0.0381362 };

        private:

            void setLocalReferenceNormal( ){
    
                asp::aspBase::setLocalReferenceNormal( N );
    
                return;
    
            }
    
            void setLocalMicroDeformation( ){
    
                asp::aspBase::setLocalMicroDeformation( chi );
    
            }

    };

    aspBaseMock asp1, asp2;

    floatMatrix dndChi( 3, floatVector( 9, 0 ) );

    floatMatrix dndN( 3, floatVector( 3, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < asp1.chi.size( ); i++ ){

        floatVector deltas( asp1.chi.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp1.chi[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.chi += deltas;

        aspm.chi -= deltas;

        for ( unsigned int j = 0; j < 3; j++ ){

            dndChi[ j ][ i ] = ( ( *aspp.getLocalCurrentNormal( ) )[ j ] - ( *aspm.getLocalCurrentNormal( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp1.N.size( ); i++ ){

        floatVector deltas( asp1.N.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp1.N[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.N += deltas;

        aspm.N -= deltas;

        for ( unsigned int j = 0; j < 3; j++ ){

            dndN[ j ][ i ] = ( ( *aspp.getLocalCurrentNormal( ) )[ j ] - ( *aspm.getLocalCurrentNormal( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dndChi, *asp1.getdLocalCurrentNormaldLocalMicroDeformation( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dndN, *asp1.getdLocalCurrentNormaldLocalReferenceNormal( ) ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalCurrentNormal ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatMatrix values = { { 1, 2, 3 }, { 4, 5, 6 } };


        private:

            void getLocalReferenceNormal( const unsigned int &index, floatVector &value ){
    
                value = values[ index ];
    
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

    floatMatrix answers = {
                              { -0.73381014, -0.454509  ,  0.50492004 },
                              { -0.68612755, -0.39784135,  0.60905767 },
                          };

    floatVector result;

    for ( unsigned int i = 0; i < answers.size( ); i++ ){

        asp.getLocalCurrentNormal( i, result );

        BOOST_CHECK( vectorTools::fuzzyEquals( result, answers[ i ] ) );

    }

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

    floatVector gradientMicroDeformation( 27, 0 );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp, gradientMicroDeformation );

    floatVector answer = { 482, 554, 626 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setCurrentDistanceVector( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    aspBaseMock aspGet;

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspGet, gradientMicroDeformation );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getCurrentDistanceVector( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setCurrentDistanceVectorGradients ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector Xi = { 1, 2, 3 };
    
            floatVector XiNL = { 4, 5, 6 };
    
            floatVector D = { .7, .8, .9 };
    
            floatVector F = { 10, 11, 12,
                              13, 14, 15,
                              16, 17, 18 };
    
            floatVector chi = { 19, 20, 21,
                                22, 23, 24,
                                25, 26, 27 };
    
            floatVector chiNLBase = { 28, 29, 30,
                                      31, 32, 33,
                                      34, 35, 36 };

        private:

            void setLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setLocalSurfaceReferenceRelativePositionVector( Xi );
    
                return;
    
            }
    
            void setNonLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setNonLocalSurfaceReferenceRelativePositionVector( XiNL );
    
                return;
    
            }
    
            void setReferenceDistanceVector( ) override {
    
                asp::aspBase::setReferenceDistanceVector( D );
    
                return;
    
            }
    
            void setLocalDeformationGradient( ) override {
    
                asp::aspBase::setLocalDeformationGradient( F );
    
                return;
    
            }
    
            void setLocalMicroDeformation( ) override {
    
                asp::aspBase::setLocalMicroDeformation( chi );
    
                return;
    
            }
    
            void setNonLocalMicroDeformationBase( ) override {
    
                asp::aspBase::setNonLocalMicroDeformationBase( chiNLBase );
    
                return;
    
            }

    };

    aspBaseMock asp1, asp2, asp3, asp4, asp5, asp6, asp7;

    floatVector gradientMicroDeformation = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                             1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
                                             1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp1, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp2, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp3, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp4, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp5, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp6, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp7, gradientMicroDeformation );

    floatMatrix dddXi( 3, floatVector( 3, 0 ) );

    floatMatrix dddD( 3, floatVector( 3, 0 ) );

    floatMatrix dddXiNL( 3, floatVector( 3, 0 ) );

    floatMatrix dddF( 3, floatVector( 9, 0 ) );

    floatMatrix dddChi( 3, floatVector( 9, 0 ) );

    floatMatrix dddChiNLBase( 3, floatVector( 9, 0 ) );

    floatMatrix dddGradChi( 3, floatVector( 27, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < asp1.Xi.size( ); i++ ){

        floatVector deltas( asp1.Xi.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp1.Xi[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.Xi += deltas;

        aspm.Xi -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            dddXi[ j ][ i ] = ( ( *aspp.getCurrentDistanceVector( ) )[ j ] - ( *aspm.getCurrentDistanceVector( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp2.D.size( ); i++ ){

        floatVector deltas( asp2.D.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp1.D[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.D += deltas;

        aspm.D -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            dddD[ j ][ i ] = ( ( *aspp.getCurrentDistanceVector( ) )[ j ] - ( *aspm.getCurrentDistanceVector( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp3.XiNL.size( ); i++ ){

        floatVector deltas( asp3.XiNL.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp3.XiNL[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.XiNL += deltas;

        aspm.XiNL -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            dddXiNL[ j ][ i ] = ( ( *aspp.getCurrentDistanceVector( ) )[ j ] - ( *aspm.getCurrentDistanceVector( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp4.F.size( ); i++ ){

        floatVector deltas( asp4.F.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp4.F[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.F += deltas;

        aspm.F -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            dddF[ j ][ i ] = ( ( *aspp.getCurrentDistanceVector( ) )[ j ] - ( *aspm.getCurrentDistanceVector( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp5.chi.size( ); i++ ){

        floatVector deltas( asp5.chi.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp5.chi[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.chi += deltas;

        aspm.chi -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            dddChi[ j ][ i ] = ( ( *aspp.getCurrentDistanceVector( ) )[ j ] - ( *aspm.getCurrentDistanceVector( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < asp6.chiNLBase.size( ); i++ ){

        floatVector deltas( asp6.chiNLBase.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp6.chiNLBase[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.chiNLBase += deltas;

        aspm.chiNLBase -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            dddChiNLBase[ j ][ i ] = ( ( *aspp.getCurrentDistanceVector( ) )[ j ] - ( *aspm.getCurrentDistanceVector( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector deltas( gradientMicroDeformation.size( ), 0 );

        deltas[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation + deltas );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation - deltas );

        for ( unsigned int j = 0; j < 3; j++ ){

            dddGradChi[ j ][ i ] = ( ( *aspp.getCurrentDistanceVector( ) )[ j ] - ( *aspm.getCurrentDistanceVector( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi, *asp1.getdCurrentDistanceVectordLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddD, *asp2.getdCurrentDistanceVectordLocalReferenceDistanceVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXiNL, *asp3.getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddF, *asp4.getdCurrentDistanceVectordLocalDeformationGradient( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddChi, *asp5.getdCurrentDistanceVectordLocalMicroDeformation( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddChiNLBase, *asp6.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddGradChi, *asp7.getdCurrentDistanceVectordGradientMicroDeformation( ) ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setCurrentDistanceVectorSecondGradients ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector Xi = { 1, 2, 3 };
    
            floatVector XiNL = { 4, 5, 6 };
    
            floatVector D = { .7, .8, .9 };
    
            floatVector F = { 10, 11, 12,
                              13, 14, 15,
                              16, 17, 18 };
    
            floatVector chi = { 19, 20, 21,
                                22, 23, 24,
                                25, 26, 27 };
    
            floatVector chiNLBase = { 28, 29, 30,
                                      31, 32, 33,
                                      34, 35, 36 };

        private:

            void setLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setLocalSurfaceReferenceRelativePositionVector( Xi );
    
                return;
    
            }
    
            void setNonLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setNonLocalSurfaceReferenceRelativePositionVector( XiNL );
    
                return;
    
            }
    
            void setReferenceDistanceVector( ) override {
    
                asp::aspBase::setReferenceDistanceVector( D );
    
                return;
    
            }
    
            void setLocalDeformationGradient( ) override {
    
                asp::aspBase::setLocalDeformationGradient( F );
    
                return;
    
            }
    
            void setLocalMicroDeformation( ) override {
    
                asp::aspBase::setLocalMicroDeformation( chi );
    
                return;
    
            }
    
            void setNonLocalMicroDeformationBase( ) override {
    
                asp::aspBase::setNonLocalMicroDeformationBase( chiNLBase );
    
                return;
    
            }

    };

    aspBaseMock asp1, asp2, asp3, asp4, asp5, asp6, asp7;

    floatVector gradientMicroDeformation = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                             1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
                                             1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp1, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp2, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp3, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp4, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp5, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp6, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp7, gradientMicroDeformation );

    floatMatrix d2ddXidXi_answer( 3, floatVector( 3 * 3, 0 ) );
    floatMatrix d2ddXidXi = d2ddXidXi_answer;

    floatMatrix d2ddDdXi_answer( 3, floatVector( 3 * 3, 0 ) );
    floatMatrix d2ddDdXi = d2ddDdXi_answer;

    floatMatrix d2ddXiNLdXi( 3, floatVector( 3 * 3, 0 ) );

    floatMatrix d2ddFdXi( 3, floatVector( 9 * 3, 0 ) );

    floatMatrix d2ddChidXi( 3, floatVector( 9 * 3, 0 ) );

    floatMatrix d2ddChiNLBasedXi_answer( 3, floatVector( 9 * 3, 0 ) );
    floatMatrix d2ddChiNLBasedXi = d2ddChiNLBasedXi_answer;

    floatMatrix d2ddGradChidXi( 3, floatVector( 27 * 3, 0 ) );

    floatMatrix d2ddDdD_answer( 3, floatVector( 3 * 3, 0 ) );
    floatMatrix d2ddDdD = d2ddDdD_answer;

    floatMatrix d2ddXiNLdD( 3, floatVector( 3 * 3, 0 ) );

    floatMatrix d2ddFdD( 3, floatVector( 9 * 3, 0 ) );

    floatMatrix d2ddChidD_answer( 3, floatVector( 9 * 3, 0 ) );
    floatMatrix d2ddChidD = d2ddChidD_answer;

    floatMatrix d2ddChiNLBasedD_answer( 3, floatVector( 9 * 3, 0 ) );
    floatMatrix d2ddChiNLBasedD = d2ddChiNLBasedD_answer;

    floatMatrix d2ddGradChidD( 3, floatVector( 27 * 3, 0 ) );

    floatMatrix d2ddXiNLdXiNL( 3, floatVector( 3 * 3, 0 ) );

    floatMatrix d2ddFdXiNL( 3, floatVector( 9 * 3, 0 ) );

    floatMatrix d2ddChidXiNL_answer( 3, floatVector( 9 * 3, 0 ) );
    floatMatrix d2ddChidXiNL = d2ddChidXiNL_answer;

    floatMatrix d2ddChiNLBasedXiNL( 3, floatVector( 9 * 3, 0 ) );

    floatMatrix d2ddGradChidXiNL( 3, floatVector( 27 * 3, 0 ) );

    floatMatrix d2ddFdF_answer( 3, floatVector( 9 * 9, 0 ) );
    floatMatrix d2ddFdF( 3, floatVector( 9 * 9, 0 ) );

    floatMatrix d2ddChidF_answer( 3, floatVector( 9 * 9, 0 ) );
    floatMatrix d2ddChidF( 3, floatVector( 9 * 9, 0 ) );

    floatMatrix d2ddChiNLBasedF_answer( 3, floatVector( 9 * 9, 0 ) );
    floatMatrix d2ddChiNLBasedF( 3, floatVector( 9 * 9, 0 ) );

    floatMatrix d2ddGradChidF_answer( 3, floatVector( 27 * 9, 0 ) );
    floatMatrix d2ddGradChidF( 3, floatVector( 27 * 9, 0 ) );

    floatMatrix d2ddChidChi_answer( 3, floatVector( 9 * 9, 0 ) );
    floatMatrix d2ddChidChi( 3, floatVector( 9 * 9, 0 ) );

    floatMatrix d2ddChiNLBasedChi_answer( 3, floatVector( 9 * 9, 0 ) );
    floatMatrix d2ddChiNLBasedChi( 3, floatVector( 9 * 9, 0 ) );

    floatMatrix d2ddGradChidChi_answer( 3, floatVector( 27 * 9, 0 ) );
    floatMatrix d2ddGradChidChi( 3, floatVector( 27 * 9, 0 ) );

    floatMatrix d2ddChiNLBasedChiNLBase_answer( 3, floatVector( 9 * 9, 0 ) );
    floatMatrix d2ddChiNLBasedChiNLBase( 3, floatVector( 9 * 9, 0 ) );

    floatMatrix d2ddGradChidChiNLBase_answer( 3, floatVector( 27 * 9, 0 ) );
    floatMatrix d2ddGradChidChiNLBase( 3, floatVector( 27 * 9, 0 ) );

    floatMatrix d2ddGradChidGradChi_answer( 3, floatVector( 27 * 27, 0 ) );
    floatMatrix d2ddGradChidGradChi( 3, floatVector( 27 * 27, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < asp1.Xi.size( ); i++ ){

        floatVector deltas( asp1.Xi.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp1.Xi[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.Xi += deltas;

        aspm.Xi -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            for ( unsigned int k = 0; k < 3; k++ ){

                d2ddXidXi[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalReferenceRelativePositionVector( ) )[ j ][ k ]
                                              - ( *aspm.getdCurrentDistanceVectordLocalReferenceRelativePositionVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 3; k++ ){

                d2ddDdXi[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalReferenceDistanceVector( ) )[ j ][ k ]
                                             - ( *aspm.getdCurrentDistanceVectordLocalReferenceDistanceVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 3; k++ ){

                d2ddXiNLdXi[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) )[ j ][ k ]
                                                - ( *aspm.getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddFdXi[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ]
                                             - ( *aspm.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChidXi[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ]
                                               - ( *aspm.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChiNLBasedXi[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ]
                                                     - ( *aspm.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 27; k++ ){

                d2ddGradChidXi[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ]
                                                   - ( *aspm.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }

    }

    for ( unsigned int i = 0; i < asp2.D.size( ); i++ ){

        floatVector deltas( asp2.D.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp1.D[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.D += deltas;

        aspm.D -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            for ( unsigned int k = 0; k < 3; k++ ){

                d2ddDdD[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalReferenceDistanceVector( ) )[ j ][ k ]
                                            - ( *aspm.getdCurrentDistanceVectordLocalReferenceDistanceVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 3; k++ ){

                d2ddXiNLdD[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) )[ j ][ k ]
                                               - ( *aspm.getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddFdD[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ]
                                            - ( *aspm.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChidD[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ]
                                              - ( *aspm.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChiNLBasedD[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ]
                                                    - ( *aspm.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 27; k++ ){

                d2ddGradChidD[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ]
                                                  - ( *aspm.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }

    }

    for ( unsigned int i = 0; i < asp3.XiNL.size( ); i++ ){

        floatVector deltas( asp3.XiNL.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp3.XiNL[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.XiNL += deltas;

        aspm.XiNL -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            for ( unsigned int k = 0; k < 3; k++ ){

                d2ddXiNLdXiNL[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) )[ j ][ k ]
                                                  - ( *aspm.getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddFdXiNL[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ]
                                               - ( *aspm.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChidXiNL[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ]
                                                 - ( *aspm.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChiNLBasedXiNL[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ]
                                                       - ( *aspm.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 27; k++ ){

                d2ddGradChidXiNL[ j ][ 3 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ]
                                                     - ( *aspm.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }

    }

    for ( unsigned int i = 0; i < asp4.F.size( ); i++ ){

        floatVector deltas( asp4.F.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp4.F[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.F += deltas;

        aspm.F -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddFdF[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ]
                                            - ( *aspm.getdCurrentDistanceVectordLocalDeformationGradient( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChidF[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ]
                                              - ( *aspm.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChiNLBasedF[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ]
                                                    - ( *aspm.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 27; k++ ){

                d2ddGradChidF[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ]
                                                  - ( *aspm.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }
    }

    for ( unsigned int i = 0; i < asp5.chi.size( ); i++ ){

        floatVector deltas( asp5.chi.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp5.chi[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.chi += deltas;

        aspm.chi -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChidChi[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ]
                                                - ( *aspm.getdCurrentDistanceVectordLocalMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChiNLBasedChi[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ]
                                                      - ( *aspm.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 27; k++ ){

                d2ddGradChidChi[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ]
                                                    - ( *aspm.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }

    }

    for ( unsigned int i = 0; i < asp6.chiNLBase.size( ); i++ ){

        floatVector deltas( asp6.chiNLBase.size( ), 0 );

        deltas[ i ] = eps * std::fabs( asp6.chiNLBase[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.chiNLBase += deltas;

        aspm.chiNLBase -= deltas;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        for ( unsigned int j = 0; j < 3; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                d2ddChiNLBasedChiNLBase[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ]
                                                            - ( *aspm.getdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

            for ( unsigned int k = 0; k < 27; k++ ){

                d2ddGradChidChiNLBase[ j ][ 9 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ]
                                                          - ( *aspm.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }

    }

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector deltas( gradientMicroDeformation.size( ), 0 );

        deltas[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation + deltas );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation - deltas );

        for ( unsigned int j = 0; j < 3; j++ ){

            for ( unsigned int k = 0; k < 27; k++ ){

                d2ddGradChidGradChi[ j ][ 27 * k + i ] = ( ( *aspp.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ]
                                                         - ( *aspm.getdCurrentDistanceVectordGradientMicroDeformation( ) )[ j ][ k ] ) / ( 2 * deltas[ i ] );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddXidXi, d2ddXidXi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddDdXi, d2ddDdXi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddXiNLdXi, *asp1.getd2CurrentDistanceVectordNonLocalReferenceRelativePositionVectordLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXi, *asp1.getd2CurrentDistanceVectordLocalDeformationGradientdLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidXi, *asp1.getd2CurrentDistanceVectordLocalMicroDeformationdLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChiNLBasedXi, d2ddChiNLBasedXi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidXi, *asp1.getd2CurrentDistanceVectordGradientMicroDeformationdLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddDdD, d2ddDdD_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddXiNLdD, *asp1.getd2CurrentDistanceVectordNonLocalReferenceRelativePositionVectordLocalReferenceDistanceVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdD, *asp1.getd2CurrentDistanceVectordLocalDeformationGradientdLocalReferenceDistanceVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidD, d2ddChidD_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChiNLBasedD, d2ddChiNLBasedD_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidD, *asp1.getd2CurrentDistanceVectordGradientMicroDeformationdLocalReferenceDistanceVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddXiNLdXiNL, *asp1.getd2CurrentDistanceVectordNonLocalReferenceRelativePositionVectordNonLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXiNL, *asp1.getd2CurrentDistanceVectordLocalDeformationGradientdNonLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidXiNL, d2ddChidXiNL_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChiNLBasedXiNL, *asp1.getd2CurrentDistanceVectordNonLocalMicroDeformationBasedNonLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidXiNL, *asp1.getd2CurrentDistanceVectordGradientMicroDeformationdNonLocalReferenceRelativePositionVector( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdF, d2ddFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidF, d2ddChidF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChiNLBasedF, d2ddChiNLBasedF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidF, d2ddGradChidF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidChi, d2ddChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChiNLBasedChi, d2ddChiNLBasedChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidChi, d2ddGradChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChiNLBasedChiNLBase, d2ddChiNLBasedChiNLBase_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidChiNLBase, d2ddGradChidChiNLBase_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidGradChi, d2ddGradChidGradChi_answer ) );

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

                asp::unit_test::aspBaseTester::set_currentDistanceVector( *this, distanceVector );

                asp::unit_test::aspBaseTester::set_localCurrentNormal( *this, localCurrentNormal );

                asp::unit_test::aspBaseTester::set_surfaceParameters( *this, surfaceParameters );

            }

    };
    
    aspBaseMock asp, aspGet;

    floatType answer = 178.2828858284305;

    floatType result;

    asp.computeSurfaceAdhesionEnergyDensity( result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    asp::dataStorage< floatType > result2;

    asp::unit_test::aspBaseTester::setSurfaceAdhesionEnergyDensity( asp, result2 );

    BOOST_CHECK( result2.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result2.second, answer ) );

    result = *aspGet.getSurfaceAdhesionEnergyDensity( );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceAdhesionEnergyDensity2 ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector Xi = { .1, .2, .3 };
    
            floatVector XiNL = { .4, .5, .6 };
    
            floatVector D = { .7, .8, .9 };
    
            floatVector F = { 1.1, 0.2, 0.3,
                              0.4, 1.5, 0.6,
                              0.7, 0.8, 1.9 };
    
            floatVector chi = { 1.19, 0.20, 0.21,
                                0.22, 1.23, 0.24,
                                0.25, 0.26, 0.27 };
    
            floatVector chiNLBase = { 1.28, 0.29, 0.30,
                                      0.31, 1.32, 0.33,
                                      0.34, 0.35, 1.36 };

            floatVector localReferenceNormal = { 0.4, .8, -0.2 };

            floatVector surfaceParameters = { 12.3, 45.6 };

        private:

            void setLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setLocalSurfaceReferenceRelativePositionVector( Xi );
    
                return;
    
            }
    
            void setNonLocalSurfaceReferenceRelativePositionVector( ) override {
    
                asp::aspBase::setNonLocalSurfaceReferenceRelativePositionVector( XiNL );
    
                return;
    
            }
    
            void setReferenceDistanceVector( ) override {
    
                asp::aspBase::setReferenceDistanceVector( D );
    
                return;
    
            }
    
            void setLocalDeformationGradient( ) override {
    
                asp::aspBase::setLocalDeformationGradient( F );
    
                return;
    
            }
    
            void setLocalMicroDeformation( ) override {
    
                asp::aspBase::setLocalMicroDeformation( chi );
    
                return;
    
            }
    
            void setNonLocalMicroDeformationBase( ) override {
    
                asp::aspBase::setNonLocalMicroDeformationBase( chiNLBase );
    
                return;
    
            }

            void setSurfaceParameters( ) override {

                asp::aspBase::setSurfaceParameters( surfaceParameters );

            }

            void setLocalReferenceNormal( ) override {

                asp::aspBase::setLocalReferenceNormal( localReferenceNormal );

            }

    };
    
    aspBaseMock asp, aspGet, asp1, asp2, asp3;

    floatVector gradientMicroDeformation = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                             1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
                                             1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspGet, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp1, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp2, gradientMicroDeformation );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp3, gradientMicroDeformation );

    floatType answer;

    floatType result;

    asp.computeSurfaceAdhesionEnergyDensity( answer );

    floatVector _dSurfaceAdhesionEnergyDensitydLocalDeformationGradient,
                _dSurfaceAdhesionEnergyDensitydLocalMicroDeformation,
                _dSurfaceAdhesionEnergyDensitydGradientMicroDeformation;

    asp.computeSurfaceAdhesionEnergyDensity( result,
                                             _dSurfaceAdhesionEnergyDensitydLocalDeformationGradient,
                                             _dSurfaceAdhesionEnergyDensitydLocalMicroDeformation,
                                             _dSurfaceAdhesionEnergyDensitydGradientMicroDeformation );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    aspGet.getdSurfaceAdhesionEnergyDensitydLocalDeformationGradient( );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, *aspGet.getSurfaceAdhesionEnergyDensity( ) ) );

    floatVector dSurfaceAdhesionEnergyDensitydLocalDeformationGradient( 9, 0 );

    floatVector dSurfaceAdhesionEnergyDensitydLocalMicroDeformation( 9, 0 );

    floatVector dSurfaceAdhesionEnergyDensitydGradientMicroDeformation( 27, 0 );
    
    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < asp1.F.size( ); i++ ){

        floatVector delta( asp1.F.size( ), 0 );

        delta[ i ] = eps * std::fabs( asp1.F[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.F += delta;

        aspm.F -= delta;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        dSurfaceAdhesionEnergyDensitydLocalDeformationGradient[ i ] = (  ( *aspp.getSurfaceAdhesionEnergyDensity( ) ) - ( *aspm.getSurfaceAdhesionEnergyDensity( ) ) ) / ( 2 * delta[ i ] );

    }

    for ( unsigned int i = 0; i < asp2.chi.size( ); i++ ){

        floatVector delta( asp2.chi.size( ), 0 );

        delta[ i ] = eps * std::fabs( asp2.chi[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        aspp.chi += delta;

        aspm.chi -= delta;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation );

        dSurfaceAdhesionEnergyDensitydLocalMicroDeformation[ i ] = (  ( *aspp.getSurfaceAdhesionEnergyDensity( ) ) - ( *aspm.getSurfaceAdhesionEnergyDensity( ) ) ) / ( 2 * delta[ i ] );

    }

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector delta( gradientMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        aspBaseMock aspp, aspm;

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspp, gradientMicroDeformation + delta );

        asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspm, gradientMicroDeformation - delta );

        dSurfaceAdhesionEnergyDensitydGradientMicroDeformation[ i ] = (  ( *aspp.getSurfaceAdhesionEnergyDensity( ) ) - ( *aspm.getSurfaceAdhesionEnergyDensity( ) ) ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dSurfaceAdhesionEnergyDensitydLocalDeformationGradient, *asp1.getdSurfaceAdhesionEnergyDensitydLocalDeformationGradient( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dSurfaceAdhesionEnergyDensitydLocalMicroDeformation, *asp2.getdSurfaceAdhesionEnergyDensitydLocalMicroDeformation( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dSurfaceAdhesionEnergyDensitydGradientMicroDeformation, *asp3.getdSurfaceAdhesionEnergyDensitydGradientMicroDeformation( ) ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceAdhesionTraction ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector distanceVector = { 1, 2, 3 };

            floatVector localCurrentNormal = { 0.45584231, 0.56980288, 0.68376346 };

            floatVector surfaceParameters = { 12.3, 45.6 };

            aspBaseMock( ){

                asp::unit_test::aspBaseTester::set_currentDistanceVector( *this, distanceVector );

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

    asp::dataStorage< floatVector > result2;

    asp::unit_test::aspBaseTester::setSurfaceAdhesionTraction( asp, result2 );

    BOOST_CHECK( result2.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result2.second, answer ) );

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

    aspBaseMock asp, aspGet;

    floatVector answer = asp.points * asp.radius;

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalReferenceSurfacePoints( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalReferenceSurfacePoints( ), answer ) );

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

    aspBaseMock asp, aspGet;

    floatVector answer = asp.points * asp.radius;

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setNonLocalReferenceSurfacePoints( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getNonLocalReferenceSurfacePoints( ), answer ) );

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

    aspBaseMock asp, aspGet;

    floatVector answer = { 1.4,  3.2,  5. ,
                           3.2,  7.7, 12.2,
                           5. , 12.2, 19.4,
                           6.8, 16.7, 26.6};

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalCurrentSurfacePoints( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalCurrentSurfacePoints( ), answer ) );

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

    aspBaseMock asp, aspGet;

    floatVector answer = { 1.4,  3.2,  5. ,
                           3.2,  7.7, 12.2,
                           5. , 12.2, 19.4,
                           6.8, 16.7, 26.6};

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setNonLocalCurrentSurfacePoints( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getNonLocalCurrentSurfacePoints( ), answer ) );

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

    aspBaseMock asp, aspGet;

    floatMatrix answer = { { 1, 10 },
                           { 2, 11 },
                           { 3, 12 } };

    asp::dataStorage< floatMatrix > result;

    asp::unit_test::aspBaseTester::setLocalParticleCurrentBoundingBox( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalParticleCurrentBoundingBox( ), answer ) );

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

    aspBaseMock asp, aspGet;

    floatMatrix answer = { { 1, 10 },
                           { 2, 11 },
                           { 3, 12 } };

    asp::dataStorage< floatMatrix > result;

    asp::unit_test::aspBaseTester::setNonLocalParticleCurrentBoundingBox( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getNonLocalParticleCurrentBoundingBox( ), answer ) );

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

    asp::aspBase asp, aspGet;

    floatVector answer = { 1, 2, 3, 4, 5, 6 };

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( asp, answer );

    asp::unit_test::aspBaseTester::set_gradientMicroDeformation( aspGet, answer );

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setLocalGradientMicroDeformation( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalGradientMicroDeformation( ), answer ) );

}


BOOST_AUTO_TEST_CASE( test_aspBase_setNonLocalMicroDeformationBase ){

    asp::aspBase asp, aspGet;

    floatVector answer = { 1, 2, 3, 4, 5, 6 };

    asp::unit_test::aspBaseTester::set_microDeformation( asp, answer );

    asp::unit_test::aspBaseTester::set_microDeformation( aspGet, answer );

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setNonLocalMicroDeformationBase( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getNonLocalMicroDeformationBase( ), answer ) );

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

    aspBaseMock asp, aspGet;

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::setSurfaceOverlapParameters( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, asp.overlapParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getSurfaceOverlapParameters( ), asp.overlapParameters ) );

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

            virtual void setLocalReferenceParticleSpacingVector( ){

                asp::unit_test::aspBaseTester::set_localReferenceParticleSpacing( *this, localReferenceParticleSpacing );

            }

            virtual void setLocalGradientMicroDeformation( ){

                asp::unit_test::aspBaseTester::set_localGradientMicroDeformation( *this, localGradientMicroDeformation );

            }

            virtual void setNonLocalMicroDeformationBase( ){

                asp::unit_test::aspBaseTester::set_nonLocalMicroDeformationBase( *this, nonLocalMicroDeformationBase );

            }

    };

    aspBaseMock asp, aspGet;

    std::unordered_map< unsigned int, floatVector > answer = { { 0, { 0, 0, 0 } }, { 1, { -0.5, 0, 0 } }, { 2, { -0.65, 0, 0 } } };

    asp::dataStorage< std::unordered_map< unsigned int, floatVector > > result;

    asp::unit_test::aspBaseTester::setParticlePairOverlap( asp, result );

    BOOST_CHECK( result.first );

    const std::unordered_map< unsigned int, floatVector > *result2 = aspGet.getParticlePairOverlap( );

    for ( auto p = answer.begin( ); p != answer.end( ); p++ ){

        auto search = result.second.find( p->first );

        BOOST_CHECK( search != result.second.end( ) );

        if ( search != result.second.end( ) ){

            BOOST_CHECK( vectorTools::fuzzyEquals( p->second, search->second ) );

        }

        auto search2 = result2->find( p->first );

        BOOST_CHECK( search2 != result2->end( ) );

        if ( search2 != result2->end( ) ){

            BOOST_CHECK( vectorTools::fuzzyEquals( p->second, search2->second ) );

        }

    }

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceOverlapEnergyDensity ){

    class aspBaseMock : public asp::aspBase{

        public:
   
            std::unordered_map< unsigned int, floatVector > particlePairOverlap = { { 0, { -0.5, 0, 0 } }, { 4, { 2, -1, 4 } } };

            floatVector surfaceOverlapParameters = { 2.3 };

            floatMatrix currentNormals = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 }, { 10, 11, 12 }, { 13, 14, 15 } };

        private:

            virtual void setParticlePairOverlap( ){

                asp::unit_test::aspBaseTester::set_particlePairOverlap( *this, particlePairOverlap );

            }

            virtual void setSurfaceOverlapParameters( ){

                asp::unit_test::aspBaseTester::set_surfaceOverlapParameters( *this, surfaceOverlapParameters );

            }

            virtual void getLocalCurrentNormal( const unsigned int &index, floatVector &normal ){

                normal = currentNormals[ index ];

            }

    };

    aspBaseMock asp, aspGet;

    std::unordered_map< unsigned int, floatType > result;

    asp::dataStorage< std::unordered_map< unsigned int, floatType > > resultSet;

    BOOST_CHECK_NO_THROW( asp.computeSurfaceOverlapEnergyDensity( result ) );

    asp::unit_test::aspBaseTester::setSurfaceOverlapEnergyDensity( asp, resultSet );

    BOOST_CHECK( resultSet.first );

    const std::unordered_map< unsigned int, floatType > *resultGet;

    resultGet = aspGet.getSurfaceOverlapEnergyDensity( );

    for ( auto p = asp.particlePairOverlap.begin( ); p != asp.particlePairOverlap.end( ); p++ ){

        auto overlapEnergy = result.find( p->first );

        BOOST_CHECK( overlapEnergy != result.end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapEnergy->second, 0.5 * asp.surfaceOverlapParameters[ 0 ] * vectorTools::dot( p->second, p->second ) * std::fabs( vectorTools::dot( asp.currentNormals[ p->first ], p->second ) ) ) );

        auto overlapEnergySet = resultSet.second.find( p->first );

        BOOST_CHECK( overlapEnergySet != resultSet.second.end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapEnergySet->second, 0.5 * asp.surfaceOverlapParameters[ 0 ] * vectorTools::dot( p->second, p->second ) * std::fabs( vectorTools::dot( asp.currentNormals[ p->first ], p->second ) ) ) );

        auto overlapEnergyGet = resultGet->find( p->first );

        BOOST_CHECK( overlapEnergyGet != resultGet->end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapEnergyGet->second, 0.5 * asp.surfaceOverlapParameters[ 0 ] * vectorTools::dot( p->second, p->second ) * std::fabs( vectorTools::dot( asp.currentNormals[ p->first ], p->second ) ) ) );

    } 

}

BOOST_AUTO_TEST_CASE( test_aspBase_computeSurfaceOverlapThickness ){

    class aspBaseMock : public asp::aspBase{

        public:
   
            std::unordered_map< unsigned int, floatVector > particlePairOverlap = { { 0, { -0.5, 0, 0 } }, { 4, { 2, -1, 4 } } };

            floatMatrix currentNormals = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 }, { 10, 11, 12 }, { 13, 14, 15 } };

        private:

            virtual void setParticlePairOverlap( ){

                asp::unit_test::aspBaseTester::set_particlePairOverlap( *this, particlePairOverlap );

            }

            virtual void getLocalCurrentNormal( const unsigned int &index, floatVector &normal ){

                normal = currentNormals[ index ];

            }

    };

    aspBaseMock asp, aspGet;

    asp::dataStorage< std::unordered_map< unsigned int, floatType > > resultSet;

    asp::unit_test::aspBaseTester::setSurfaceOverlapThickness( asp, resultSet );

    BOOST_CHECK( resultSet.first );

    const std::unordered_map< unsigned int, floatType > *resultGet;

    resultGet = aspGet.getSurfaceOverlapThickness( );

    for ( auto p = asp.particlePairOverlap.begin( ); p != asp.particlePairOverlap.end( ); p++ ){

        auto overlapThicknessSet = resultSet.second.find( p->first );

        BOOST_CHECK( overlapThicknessSet != resultSet.second.end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapThicknessSet->second, std::fabs( vectorTools::dot( asp.currentNormals[ p->first ], p->second ) ) ) );

        auto overlapThicknessGet = resultGet->find( p->first );

        BOOST_CHECK( overlapThicknessGet != resultGet->end( ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( overlapThicknessGet->second, std::fabs( vectorTools::dot( asp.currentNormals[ p->first ], p->second ) ) ) );

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

BOOST_AUTO_TEST_CASE( test_aspBase_getNumLocalParticles ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkNumLocalParticles( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getPreviousTime ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkPreviousTime( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getDeltaTime ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkDeltaTime( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getPreviousTemperature ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkPreviousTemperature( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getTemperature ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkTemperature( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getDeformationGradient ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkDeformationGradient( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getPreviousDeformationGradient ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkPreviousDeformationGradient( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getPreviousMicroDeformation ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkPreviousMicroDeformation( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getPreviousStateVariables ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkPreviousStateVariables( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getPreviousLocalStateVariables ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkPreviousLocalStateVariables( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getParticleParameters ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkParticleParameters( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalIndex ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkLocalIndex( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getNonLocalIndex ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkNonLocalIndex( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalSurfaceNodeIndex ){

    asp::aspBase asp;

    asp::unit_test::aspBaseTester::checkLocalSurfaceNodeIndex( asp );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalParticleEnergy ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatType volume = 2.4;
    
            floatType energyDensity = 3.5;
    
            floatVector stress = { 1, 2, 3 };
    
            floatVector stateVariables = { -1, -2 };
    
            floatType probabilityRatio = .34;

        private:

            void setLocalParticleCurrentVolume( ){
    
                asp::unit_test::aspBaseTester::set_localParticleCurrentVolume( *this, volume );
    
            }

            void setLocalParticleQuantities( ){
    
                asp::unit_test::aspBaseTester::set_localParticleEnergyDensity( *this, energyDensity );

                asp::unit_test::aspBaseTester::set_localParticleMicroCauchyStress( *this, stress );

                asp::unit_test::aspBaseTester::set_localParticleStateVariables( *this, stateVariables );

                asp::unit_test::aspBaseTester::set_localParticleLogProbabilityRatio( *this, probabilityRatio );
    
            }

    };

    aspBaseMock asp, aspGet;

    floatType answer = 8.4;

    asp::dataStorage< floatType > result;

    asp::unit_test::aspBaseTester::setLocalParticleEnergy( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalParticleEnergy( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalParticleEnergyDensity ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatType volume = 2.4;
    
            floatType energyDensity = 3.5;

            floatVector stress = { 1, 2, 3 };
    
            floatVector stateVariables = { -1, -2 };
    
            floatType probabilityRatio = .34;

        private:

            void setLocalParticleQuantities( ){
    
                asp::unit_test::aspBaseTester::set_localParticleEnergyDensity( *this, energyDensity );
    
                asp::unit_test::aspBaseTester::set_localParticleMicroCauchyStress( *this, stress );

                asp::unit_test::aspBaseTester::set_localParticleStateVariables( *this, stateVariables );

                asp::unit_test::aspBaseTester::set_localParticleLogProbabilityRatio( *this, probabilityRatio );
    
            }

    };

    aspBaseMock asp, aspGet;

    floatType answer = asp.energyDensity;

    asp::unit_test::aspBaseTester::set_localParticleEnergyDensity( asp, asp.energyDensity );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalParticleEnergyDensity( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalParticleLogProbabilityRatio ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatType volume = 2.4;
    
            floatType energyDensity = 3.5;

            floatVector stress = { 1, 2, 3 };
    
            floatVector stateVariables = { -1, -2 };
    
            floatType probabilityRatio = .34;

        private:

            void setLocalParticleQuantities( ){
    
                asp::unit_test::aspBaseTester::set_localParticleEnergyDensity( *this, energyDensity );
    
                asp::unit_test::aspBaseTester::set_localParticleMicroCauchyStress( *this, stress );

                asp::unit_test::aspBaseTester::set_localParticleStateVariables( *this, stateVariables );

                asp::unit_test::aspBaseTester::set_localParticleLogProbabilityRatio( *this, probabilityRatio );
    
            }

    };

    aspBaseMock asp, aspGet;

    floatType answer = asp.probabilityRatio;

    asp::unit_test::aspBaseTester::set_localParticleLogProbabilityRatio( asp, asp.probabilityRatio );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalParticleLogProbabilityRatio( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalParticleQuantities ){

    class aspBaseMock : public asp::aspBase{

        private:

            virtual void computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
                                                            const floatVector &currentMicroDeformation, const floatVector &previousMicroDeformation,
                                                            const floatType &currentTemperature, const floatType &previousTemperature,
                                                            const floatVector &previousStateVariables,
                                                            const floatVector &parameters,
                                                            floatType &energyDensity, floatVector &cauchyStress, floatVector &stateVariables, floatType &logProbabilityRatio ){

                energyDensity = 1.2;

                cauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

                stateVariables = { -1, -2 };

                logProbabilityRatio = 0.47;

                return;

            }

    };

    aspBaseMock asp, aspGet1, aspGet2, aspGet3, aspGet4;

    floatType energyDensityAnswer = 1.2;

    floatVector cauchyStressAnswer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    floatVector stateVariablesAnswer = { -1, -2 };

    floatType logProbabilityRatioAnswer = 0.47;

    asp::dataStorage< floatType > energyDensityResult;

    asp::dataStorage< floatVector > cauchyStressResult;

    asp::dataStorage< floatVector > stateVariableResult;

    asp::dataStorage< floatType > logProbabilityRatioResult;

    asp::unit_test::aspBaseTester::setLocalParticleQuantities( asp, energyDensityResult, cauchyStressResult, stateVariableResult, logProbabilityRatioResult );

    BOOST_CHECK( energyDensityResult.first );

    BOOST_CHECK( cauchyStressResult.first );

    BOOST_CHECK( stateVariableResult.first );

    BOOST_CHECK( logProbabilityRatioResult.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( energyDensityResult.second, energyDensityAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cauchyStressResult.second, cauchyStressAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( stateVariableResult.second, stateVariablesAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( logProbabilityRatioResult.second, logProbabilityRatioAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet1.getLocalParticleEnergyDensity( ), energyDensityAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet2.getLocalParticleMicroCauchyStress( ), cauchyStressAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet3.getLocalParticleStateVariables( ), stateVariablesAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet4.getLocalParticleLogProbabilityRatio( ), logProbabilityRatioAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalParticleReferenceVolume ){

    class aspBaseMock : public asp::aspBase{

        floatType radius = 2.4;

        private:

            void setLocalReferenceRadius( ){
    
                asp::unit_test::aspBaseTester::set_localReferenceRadius( *this, radius );
    
            }

    };

    aspBaseMock asp, aspGet;

    floatType answer = 57.90583579096705;

    asp::dataStorage< floatType > result;

    asp::unit_test::aspBaseTester::setLocalParticleReferenceVolume( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalParticleReferenceVolume( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_getLocalParticleCurrentVolume ){

    class aspBaseMock : public asp::aspBase{

        floatType referenceVolume = 2.4;

        floatVector microDeformation = { 1, 2, 3, 4, 5, 6, 7, 8, 2 };

        private:

            void setLocalParticleReferenceVolume( ){
    
                asp::unit_test::aspBaseTester::set_localParticleReferenceVolume( *this, referenceVolume );
    
            }

            void setLocalMicroDeformation( ){
    
                asp::unit_test::aspBaseTester::set_localMicroDeformation( *this, microDeformation );
    
            }

    };

    aspBaseMock asp, aspGet;

    floatType answer = 50.39999999999997;

    asp::dataStorage< floatType > result;

    asp::unit_test::aspBaseTester::setLocalParticleCurrentVolume( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalParticleCurrentVolume( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setLocalParticleParameters ){

    class aspBaseMock : public asp::aspBase{

        floatVector particleParameters = { 1, 2, 3 };

    };

    aspBaseMock asp, aspGet;

    floatVector answer = { 1, 2, 3 };

    asp::dataStorage< floatVector > result;

    asp::unit_test::aspBaseTester::set_particleParameters( asp, answer );

    asp::unit_test::aspBaseTester::set_particleParameters( aspGet, answer );

    asp::unit_test::aspBaseTester::setLocalParticleParameters( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getLocalParticleParameters( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_assembleLocalParticle ){

    class aspBaseMock : public asp::aspBase{

        public:

            unsigned int numLocalParticles = 4;
    
            floatVector energies = { 1, 2, 3, 4 };

            floatVector energyDensities = { 10, 20, 30, 40 };
    
            floatMatrix microCauchyStresses = { { 1, 2, 3 }, { 4, 5 }, { 6, 7, 8 }, { 9, 10, 11, 12 } };
    
            floatVector localParticleVolumes = { 0.1, 0.2, 0.3, 0.4 };
    
            floatVector probabilityRatios = { -0.2, -0.3, -0.4, -0.5 };
    
            aspBaseMock( ) : aspBase( ){
    
                asp::unit_test::aspBaseTester::set_numLocalParticles( *this, numLocalParticles );
    
            }

        private:

            virtual void setLocalParticleEnergy( ){

                const unsigned int* localIndex = getLocalIndex( );

                asp::unit_test::aspBaseTester::set_localParticleEnergy( *this, energies[ *localIndex ] );

            }

            virtual void setLocalParticleQuantities( ){
    
                const unsigned int* localIndex = getLocalIndex( );
    
                asp::unit_test::aspBaseTester::set_localParticleEnergyDensity( *this, energies[ *localIndex ] );

                asp::unit_test::aspBaseTester::set_localParticleMicroCauchyStress( *this, microCauchyStresses[ *localIndex ] );
    
                asp::unit_test::aspBaseTester::set_localParticleCurrentVolume( *this, localParticleVolumes[ *localIndex ] );
    
                asp::unit_test::aspBaseTester::set_localParticleLogProbabilityRatio( *this, probabilityRatios[ *localIndex ] );
    
            }

    };

    aspBaseMock asp1, asp2, asp3, asp4;

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp1.getAssembledLocalParticleEnergies( ), asp1.energies ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp2.getAssembledLocalParticleMicroCauchyStresses( ), asp2.microCauchyStresses ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp3.getAssembledLocalParticleVolumes( ), asp3.localParticleVolumes ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp4.getAssembledLocalParticleLogProbabilityRatios( ), asp4.probabilityRatios ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_setSurfaceAdhesionThickness ){

    class aspBaseMock : public asp::aspBase{

        public:

            floatVector normalVector = { 1, 2, 3 };

            floatVector distanceVector = { 4, 5, 6 };

            aspBaseMock( ){

                normalVector /= vectorTools::l2norm( normalVector );

            }

        private:

            virtual void setCurrentDistanceVector( ){

                asp::unit_test::aspBaseTester::set_currentDistanceVector( *this, distanceVector );

            }

            virtual void setLocalCurrentNormal( ){

                asp::unit_test::aspBaseTester::set_localCurrentNormal( *this, normalVector );

            }

    };

    aspBaseMock asp, aspGet;

    floatType answer = std::fabs( vectorTools::dot( asp.normalVector, asp.distanceVector ) );

    asp::dataStorage< floatType > result;

    asp::unit_test::aspBaseTester::setSurfaceAdhesionThickness( asp, result );

    BOOST_CHECK( result.first );

    BOOST_CHECK( vectorTools::fuzzyEquals( result.second, answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *aspGet.getSurfaceAdhesionThickness( ), answer ) );

}

BOOST_AUTO_TEST_CASE( test_aspBase_assembleSurfaceEnergies ){

    class aspBaseMock : public asp::aspBase{

        public:

            unsigned int numLocalParticles = 4;

            unsigned int numSurfacePoints = 2;

            floatVector unitSpherePoints = { 1, 2, 3, 4, 5, 6 };

            std::vector< unsigned int > unitSphereConnectivity = { 10, 11, 12, 13 };

            std::vector< std::vector< std::vector< std::vector< unsigned int > > > > neighboringOverlapParticles
                {
                    {
                        {
                            { 1, 7 },
                            { 0, },
                            { },
                            { 13, 8, 17},
                        },
                        {
                            { 11, 4 },
                            { 9, 34 },
                            { 10 },
                            { 8, 9, 10, 11},
                        },
                    },
                    {
                        {
                            { 8, 13},
                            { 2 },
                            { 0, 1, 2 },
                            { 13 },
                        },
                        {
                            { },
                            { },
                            { 4, 7 },
                            { 45 },
                        },
                    },
                    {
                        {
                            { 6, 7, 8, 9 },
                            { 1, },
                            { 7, 6 },
                            { 35 },
                        },
                        {
                            { 1, 2 },
                            { },
                            { 6 },
                            { },
                        },
                    },
                    {
                        {
                            { 78, 22, 45, 13 },
                            { 9, 87 },
                            { },
                            { 22 },
                        },
                        {
                            { },
                            { },
                            { },
                            { 35, 6, 7 },
                        },
                    },
                };

            std::vector< std::vector< floatVector > > adhesionEnergyDensities =
                {
                    {
                        {
                            1, 2, 3, 4
                        },
                        {
                            5, 6, 7, 8
                        },
                    },
                    {
                        {
                            9, 10, 11, 12
                        },
                        {
                            13, 14, 15, 16
                        },
                    },
                    {
                        {
                            17, 18, 19, 20
                        },
                        {
                            21, 22, 23, 24
                        },
                    },
                    {
                        {
                            25, 26, 27, 28
                        },
                        {
                            29, 30, 31, 32
                        },
                    },
                };

            std::vector< std::vector< floatVector > > adhesionThicknesses =
            {
                {
                    {
                        33, 34, 35, 36
                    },
                    {
                        37, 38, 39, 40
                    },
                },
                {
                    {
                        41, 42, 43, 44
                    },
                    {
                        45, 46, 47, 48
                    },
                },
                {
                    {
                        49, 50, 51, 52
                    },
                    {
                        53, 54, 55, 56
                    },
                },
                {
                    {
                        57, 58, 59, 60
                    },
                    {
                        61, 62, 63, 64
                    },
                },
            };

            std::vector< std::vector< floatMatrix > > adhesionTractions =
            {
                {
                    {
                        {  1,  2,  3 },
                        {  4,  5,  6 },
                        {  7,  8,  9 },
                        { 10, 11, 12 },
                    },
                    {
                        { 13, 14, 15 },
                        { 16, 17, 18 },
                        { 19, 20, 21 },
                        { 22, 23, 24 },
                    },
                },
                {
                    {
                        { 25, 26, 27 },
                        { 28, 29, 30 },
                        { 31, 32, 33 },
                        { 34, 35, 36 },
                    },
                    {
                        { 37, 38, 39 },
                        { 40, 41, 42 },
                        { 43, 44, 45 },
                        { 46, 47, 48 },
                    },
                },
                {
                    {
                        { 49, 50, 51 },
                        { 52, 53, 54 },
                        { 55, 56, 57 },
                        { 58, 59, 60 },
                    },
                    {
                        { 61, 62, 63 },
                        { 64, 65, 66 },
                        { 67, 68, 69 },
                        { 70, 71, 72 },
                    },
                },
                {
                    {
                        { 73, 74, 75 },
                        { 76, 77, 78 },
                        { 79, 80, 81 },
                        { 82, 83, 84 },
                    },
                    {
                        { 85, 86, 87 },
                        { 88, 89, 90 },
                        { 91, 92, 93 },
                        { 94, 95, 96 },
                    },
                },
            };

            std::vector< std::vector< std::vector< std::unordered_map< unsigned int, floatType > > > > overlapEnergyDensities;

            std::vector< std::vector< std::vector< std::unordered_map< unsigned int, floatType > > > > overlapThicknesses;

            std::vector< std::vector< std::vector< std::unordered_map< unsigned int, floatVector > > > > overlapTractions;
    
            aspBaseMock( ) : aspBase( ){
    
                asp::unit_test::aspBaseTester::set_numLocalParticles( *this, numLocalParticles );
   
                initializeOverlapQuantitiesTest( );
 
            }

        private:

            void initializeOverlapQuantitiesTest( ){

                overlapEnergyDensities.resize( numLocalParticles );

                overlapThicknesses.resize( numLocalParticles );

                overlapTractions.resize( numLocalParticles );

                for ( unsigned int i = 0; i < numLocalParticles; i++ ){

                    overlapEnergyDensities[ i ].resize( numSurfacePoints );
    
                    overlapThicknesses[ i ].resize( numSurfacePoints );
    
                    overlapTractions[ i ].resize( numSurfacePoints );

                    for ( unsigned int j = 0; j < numSurfacePoints; j++ ){

                        overlapEnergyDensities[ i ][ j ].resize( numLocalParticles );

                        overlapThicknesses[ i ][ j ].resize( numLocalParticles );

                        overlapTractions[ i ][ j ].resize( numLocalParticles );

                        for ( unsigned int k = 0; k < numLocalParticles; k++ ){

                            for ( auto index = neighboringOverlapParticles[ i ][ j ][ k ].begin( );
                                       index != neighboringOverlapParticles[ i ][ j ][ k ].end( ); index++ ){

                                overlapEnergyDensities[ i ][ j ][ k ].emplace( *index, adhesionEnergyDensities[ i ][ j ][ k ] * ( *index ) * .142 );

                                overlapThicknesses[ i ][ j ][ k ].emplace( *index, adhesionThicknesses[ i ][ j ][ k ] * ( ( *index ) - 3 ) * .89 );

                                overlapTractions[ i ][ j ][ k ].emplace( *index, adhesionTractions[ i ][ j ][ k ] * ( ( *index ) - .78 ) * .09 );

                            }

                        }

                    }

                }

            }

            virtual void initializeUnitSphere( ){

                asp::unit_test::aspBaseTester::set_unitSphere( *this, unitSpherePoints, unitSphereConnectivity );

            }

            virtual void setSurfaceAdhesionEnergyDensity( ){
    
                const unsigned int* localIndex = getLocalIndex( );

                const unsigned int* localSurfaceNodeIndex = getLocalSurfaceNodeIndex( );

                const unsigned int* nonLocalIndex = getNonLocalIndex( );
    
                asp::unit_test::aspBaseTester::set_surfaceAdhesionEnergyDensity( *this, adhesionEnergyDensities[ *localIndex ][ *localSurfaceNodeIndex ][ *nonLocalIndex ] );
    
            }

            virtual void setSurfaceAdhesionTraction( ){
    
                const unsigned int* localIndex = getLocalIndex( );

                const unsigned int* localSurfaceNodeIndex = getLocalSurfaceNodeIndex( );

                const unsigned int* nonLocalIndex = getNonLocalIndex( );
    
                asp::unit_test::aspBaseTester::set_surfaceAdhesionTraction( *this, adhesionTractions[ *localIndex ][ *localSurfaceNodeIndex ][ *nonLocalIndex ] );
    
            }

            virtual void setSurfaceAdhesionThickness( ){
    
                const unsigned int* localIndex = getLocalIndex( );

                const unsigned int* localSurfaceNodeIndex = getLocalSurfaceNodeIndex( );

                const unsigned int* nonLocalIndex = getNonLocalIndex( );
    
                asp::unit_test::aspBaseTester::set_surfaceAdhesionThickness( *this, adhesionThicknesses[ *localIndex ][ *localSurfaceNodeIndex ][ *nonLocalIndex ] );
    
            }

            virtual void setSurfaceOverlapEnergyDensity( ){
    
                const unsigned int* localIndex = getLocalIndex( );

                const unsigned int* localSurfaceNodeIndex = getLocalSurfaceNodeIndex( );

                const unsigned int* nonLocalIndex = getNonLocalIndex( );
    
                asp::unit_test::aspBaseTester::set_surfaceOverlapEnergyDensity( *this, overlapEnergyDensities[ *localIndex ][ *localSurfaceNodeIndex ][ *nonLocalIndex ] );
    
            }

            virtual void setSurfaceOverlapTraction( ){
    
                const unsigned int* localIndex = getLocalIndex( );

                const unsigned int* localSurfaceNodeIndex = getLocalSurfaceNodeIndex( );

                const unsigned int* nonLocalIndex = getNonLocalIndex( );
    
                asp::unit_test::aspBaseTester::set_surfaceOverlapTraction( *this, overlapTractions[ *localIndex ][ *localSurfaceNodeIndex ][ *nonLocalIndex ] );
    
            }

            virtual void setSurfaceOverlapThickness( ){
    
                const unsigned int* localIndex = getLocalIndex( );

                const unsigned int* localSurfaceNodeIndex = getLocalSurfaceNodeIndex( );

                const unsigned int* nonLocalIndex = getNonLocalIndex( );
    
                asp::unit_test::aspBaseTester::set_surfaceOverlapThickness( *this, overlapThicknesses[ *localIndex ][ *localSurfaceNodeIndex ][ *nonLocalIndex ] );
    
            }

    };

    aspBaseMock asp1, asp2, asp3, asp4, asp5, asp6;

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp1.getAssembledSurfaceAdhesionEnergyDensities( ), asp1.adhesionEnergyDensities ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp2.getAssembledSurfaceAdhesionTractions( ), asp2.adhesionTractions ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *asp3.getAssembledSurfaceAdhesionThicknesses( ), asp3.adhesionThicknesses ) );

    for ( unsigned int i = 0; i < asp1.numLocalParticles; i++ ){

        for ( unsigned int j = 0; j < asp1.numSurfacePoints; j++ ){

            for ( unsigned int k = 0; k < asp1.numLocalParticles; k++ ){

                for ( auto index = asp1.neighboringOverlapParticles[ i ][ j ][ k ].begin( );
                           index != asp1.neighboringOverlapParticles[ i ][ j ][ k ].end( );
                           index++ ){

                    auto energy    = ( *asp4.getAssembledSurfaceOverlapEnergyDensities( ) )[ i ][ j ][ k ].find( *index );

                    auto thickness = ( *asp5.getAssembledSurfaceOverlapThicknesses( ) )[ i ][ j ][ k ].find( *index );

                    auto traction  = ( *asp6.getAssembledSurfaceOverlapTractions( ) )[ i ][ j ][ k ].find( *index );

                    BOOST_CHECK( energy != ( *asp4.getAssembledSurfaceOverlapEnergyDensities( ) )[ i ][ j ][ k ].end( ) );

                    BOOST_CHECK( thickness != ( *asp5.getAssembledSurfaceOverlapThicknesses( ) )[ i ][ j ][ k ].end( ) );

                    BOOST_CHECK( traction != ( *asp6.getAssembledSurfaceOverlapTractions( ) )[ i ][ j ][ k ].end( ) );

                    BOOST_CHECK( vectorTools::fuzzyEquals( energy->second, asp4.overlapEnergyDensities[ i ][ j ][ k ][ *index ] ) );

                    BOOST_CHECK( vectorTools::fuzzyEquals( thickness->second, asp5.overlapThicknesses[ i ][ j ][ k ][ *index ] ) );

                    BOOST_CHECK( vectorTools::fuzzyEquals( traction->second, asp6.overlapTractions[ i ][ j ][ k ][ *index ] ) );

                }

            }

        }

    }

}

