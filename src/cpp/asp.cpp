/**
  ******************************************************************************
  * \file asp.cpp
  ******************************************************************************
  * A C++ library for printing messages to stdout. Used as a stub repo example.
  ******************************************************************************
  */

#include<asp.h>
#include<stress_tools.h>
#include<traction_separation.h>
#include<surface_integration.h>

namespace asp{

    //Define asp global constants in a place that Doxygen can pick up for documentation
    /** \brief Define the expected number of tensor spatial dimensions for the Abaqus interface. */
    const int spatialDimensions = 3;

    aspBase::aspBase( ){
        /*!
         * The default constructor for ASP
         */

        // Initialize surface integral pairs
        _localReferenceRadius = std::make_pair( false, 0 );

        _nonLocalReferenceRadius = std::make_pair( false, 0 );

        _unitSpherePoints = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localReferenceNormal = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _nonLocalSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _referenceDistanceVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localReferenceParticleSpacing = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localDeformationGradient = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localMicroDeformation = std::make_pair( false, floatVector( _dimension, 0 ) );

        _nonLocalMicroDeformation = std::make_pair( false, floatVector( _dimension, 0 ) );

        _currentDistanceVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localCurrentNormal = std::make_pair( false, floatVector( _dimension, 0 ) );

    }

    void aspBase::computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
                                                     const floatVector &currentMicroDeformation, const floatVector &previousMicroDeformation,
                                                     const floatType &currentTemperature, const floatType &previousTemperature,
                                                     const floatVector &previousStateVariables,
                                                     const floatVector &parameters,
                                                     floatType &energy, floatVector &cauchyStress ){
        /*!
         * Default function for computing the local particle's energy per unit current volume
         * Defaults to a linear-elastic solid with parameters defined in the reference configuration.
         * 
         * \param &previousTime: The previous value of the time
         * \param &deltaTime: The change in time
         * \param &currentMicroDeformation: The current value of the micro-deformation tensor \f$\chi\f$
         * \param &previousMicroDeformation: The previous value of the micro-deformation tensor \f$\chi\f$
         * \param &currentTemperature: The current value of the temperature
         * \param &previousTemperature: The previous value of the temperature
         * \param &parameters: The parameters for the model
         * \param &energy: The current value of the Helmholtz free energy per unit current volume
         * \param &cauchyStress: The current value of the Cauchy stress
         */

        ERROR_TOOLS_CATCH_NODE_POINTER( stressTools::linearElasticity::evaluateEnergy( currentMicroDeformation, parameters, energy, cauchyStress ) );

        return;

    }

    void aspBase::computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
                                                     const floatVector &currentMicroDeformation, const floatVector &previousMicroDeformation,
                                                     const floatType &currentTemperature, const floatType &previousTemperature,
                                                     const floatVector &previousStateVariables,
                                                     const floatVector &parameters,
                                                     floatType &energy, floatVector &cauchyStress, floatType &logProbabilityRatio ){
        /*!
         * Default function for computing the local particle's energy per unit current volume
         * Defaults to a linear-elastic solid with parameters defined in the reference configuration.
         * 
         * \param &previousTime: The previous value of the time
         * \param &deltaTime: The change in time
         * \param &currentMicroDeformation: The current value of the micro-deformation tensor \f$\chi\f$
         * \param &previousMicroDeformation: The previous value of the micro-deformation tensor \f$\chi\f$
         * \param &currentTemperature: The current value of the temperature
         * \param &previousTemperature: The previous value of the temperature
         * \param &parameters: The parameters for the model
         * \param &energy: The current value of the Helmholtz free energy per unit current volume
         * \param &cauchyStress: The current value of the Cauchy stress
         * \param &logProbabilityRatio: The natural logarithm of the probability ratio between the new
         *     probability density and the previous value.
         */

        ERROR_TOOLS_CATCH_NODE_POINTER( stressTools::linearElasticity::evaluateEnergy( currentMicroDeformation, parameters, energy, cauchyStress ) );

        logProbabilityRatio = 0.;

        return;

    }

    void aspBase::initializeUnitSphere( ){
        /*!
         * Initialize the unit sphere (i.e. a sphere of radius 1) to integrate over
         */

        surfaceIntegration::decomposeSphere( 1.0, _surfaceElementCount, _unitSpherePoints.second, _unitSphereConnectivity.second );

        _unitSpherePoints.first = true;

        _unitSphereConnectivity.first = true;

        return;

    }

    const floatVector* aspBase::getUnitSpherePoints( ){
        /*!
         * Get the points on the unit sphere
         */

        if ( !_unitSpherePoints.first ){

            ERROR_TOOLS_CATCH( initializeUnitSphere( ) );

        }

        return &_unitSpherePoints.second;

    }

    const std::vector< unsigned int >* aspBase::getUnitSphereConnectivity( ){
        /*!
         * Get the points on the unit sphere
         */

        if ( !_unitSphereConnectivity.first ){

            ERROR_TOOLS_CATCH( initializeUnitSphere( ) );

        }

        return &_unitSphereConnectivity.second;

    }

    void aspBase::setLocalReferenceNormal( ){
        /*!
         * Set the local reference normal vector
         * 
         * Because we have a unit sphere centered on zero as a basis the normal vector is the same
         * as any position vector on the surface of the sphere.
         */

        const floatVector* unitSpherePoints;
        ERROR_TOOLS_CATCH( unitSpherePoints = getUnitSpherePoints( ) );

        if ( _dimension * ( _localIndex + 1 ) > unitSpherePoints->size( ) ){

            std::string message = "The requested index is greater than the number of points available on the unit sphere.\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";
            message            += "  number of points:      " + std::to_string( unitSpherePoints->size( ) / _dimension );

            ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        _localReferenceNormal.second = floatVector( unitSpherePoints->begin( ) + _dimension * _localSurfaceNodeIndex,
                                                    unitSpherePoints->begin( ) + _dimension * ( _localSurfaceNodeIndex + 1 ) );

        _localReferenceNormal.second /= vectorTools::l2norm( _localReferenceNormal.second );

        _localReferenceNormal.first = true;

        return;

    }

    const floatVector* aspBase::getLocalReferenceNormal( ){
        /*!
         * Get the local reference normal vector
         * 
         * Get the local reference normal vector or set it and return it if necessary
         */

        if ( !_localReferenceNormal.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceNormal( ) );

        }

        return &_localReferenceNormal.second;

    }

    void aspBase::setLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Set the local surface relative position vector
         * 
         * \f$\Xi_I = R^{local} N_I\f$
         * 
         * where \f$R^{local}\f$ is the local radius and \f$N_I\f$ is the local reference normal.
         */

        const floatVector* localReferenceNormal;
        ERROR_TOOLS_CATCH( localReferenceNormal = getLocalReferenceNormal( ) );

        const floatType* localReferenceRadius;
        ERROR_TOOLS_CATCH( localReferenceRadius = getLocalReferenceRadius( ) );

        _localSurfaceReferenceRelativePositionVector.second = ( *localReferenceRadius ) * ( *localReferenceNormal );

        _localSurfaceReferenceRelativePositionVector.first = true;

        return;

    }

    const floatVector* aspBase::getLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Get the value of the local surface reference relative position vector
         */

        if ( !_localSurfaceReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setLocalSurfaceReferenceRelativePositionVector( ) );

        }

        return &_localSurfaceReferenceRelativePositionVector.second;

    }

    void aspBase::setNonLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Set the non-local relative position vector
         * 
         * \f$\Xi_I^{non-local} = -R^{non-local} N_I\f$
         * 
         * where \f$R^{non-local}\f$ is the non-local radius and \f$N_I\f$ is the local reference normal.
         */

        const floatVector* localReferenceNormal;
        ERROR_TOOLS_CATCH( localReferenceNormal = getLocalReferenceNormal( ) );

        const floatType* nonLocalReferenceRadius;
        ERROR_TOOLS_CATCH( nonLocalReferenceRadius = getNonLocalReferenceRadius( ) );

        _nonLocalSurfaceReferenceRelativePositionVector.second = -( *nonLocalReferenceRadius ) * ( *localReferenceNormal );

        _nonLocalSurfaceReferenceRelativePositionVector.first = true;

        return;

    }

    const floatVector* aspBase::getNonLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Get the non-local surface reference relative position vector
         */

        if ( !_nonLocalSurfaceReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setNonLocalSurfaceReferenceRelativePositionVector( ) );

        }

        return &_nonLocalSurfaceReferenceRelativePositionVector.second;

    }

    void aspBase::setLocalReferenceRadius( ){
        /*!
         * Set the local reference radius
         */

        _localReferenceRadius.second = _radius;

        _localReferenceRadius.first = true;

        return;

    }

    const floatType* aspBase::getLocalReferenceRadius( ){
        /*!
         * Get the local reference radius
         */

        if ( !_localReferenceRadius.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceRadius( ) );

        }

        return &_localReferenceRadius.second;

    }

    void aspBase::setNonLocalReferenceRadius( ){
        /*!
         * Set the non-local reference radius
         */

        _nonLocalReferenceRadius.second = _radius;

        _nonLocalReferenceRadius.first = true;

        return;

    }

    const floatType* aspBase::getNonLocalReferenceRadius( ){
        /*!
         * Get the local reference radius
         */

        if ( !_nonLocalReferenceRadius.first ){

            ERROR_TOOLS_CATCH( setNonLocalReferenceRadius( ) );

        }

        return &_nonLocalReferenceRadius.second;

    }


    void aspBase::setLocalDeformationGradient( ){
        /*!
         * Set the local deformation gradient
         */

        _localDeformationGradient.second = _deformationGradient;

        _localDeformationGradient.first = true;

        return;

    }

    const floatVector* aspBase::getLocalDeformationGradient( ){
        /*!
         * Get the local deformation gradient
         */

        if ( !_localDeformationGradient.first ){

            ERROR_TOOLS_CATCH( setLocalDeformationGradient( ) );

        }

        return &_localDeformationGradient.second;

    }

    void aspBase::setLocalMicroDeformation( ){
        /*!
         * Set the local micro deformation
         */

        _localMicroDeformation.second = _microDeformation;

        _localMicroDeformation.first = true;

        return;

    }

    const floatVector* aspBase::getLocalMicroDeformation( ){
        /*!
         * Get the local micro deformation
         */

        if ( !_localMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setLocalMicroDeformation( ) );

        }

        return &_localMicroDeformation.second;

    }

    void aspBase::setLocalReferenceParticleSpacing( ){
        /*!
         * Set the local particle spacing
         * 
         * \f$dX_I = \Xi_I^{local} + D_I - \Xi_I^{non-local}\f$
         */

        const floatVector* localSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( localSurfaceReferenceRelativePositionVector = getLocalSurfaceReferenceRelativePositionVector( ) );

        const floatVector* nonLocalSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( nonLocalSurfaceReferenceRelativePositionVector = getNonLocalSurfaceReferenceRelativePositionVector( ) );

        const floatVector* referenceDistanceVector;
        ERROR_TOOLS_CATCH( referenceDistanceVector = getReferenceDistanceVector( ) );

        _localReferenceParticleSpacing.second = ( *localSurfaceReferenceRelativePositionVector )
                                              + ( *referenceDistanceVector )
                                              - ( *nonLocalSurfaceReferenceRelativePositionVector );

        _localReferenceParticleSpacing.first = true;

        return;

    }

    const floatVector* aspBase::getLocalReferenceParticleSpacing( ){
        /*!
         * Get the local reference particle spacing
         */

        if ( !_localReferenceParticleSpacing.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceParticleSpacing( ) );

        }

        return &_localReferenceParticleSpacing.second;

    }

    void aspBase::setNonLocalMicroDeformation( ){
        /*!
         * Set the non-local micro deformation
         * 
         * \f$ \chi_{iI}^{non-local} = \chi_{iI} + \chi_{iI,J} dX_J\f$
         */

        const floatVector* localReferenceParticleSpacing;
        ERROR_TOOLS_CATCH( localReferenceParticleSpacing = getLocalReferenceParticleSpacing( ) );

        const floatVector* microDeformationBase;
        ERROR_TOOLS_CATCH( microDeformationBase = getNonLocalMicroDeformationBase( ) );

        const floatVector* gradientMicroDeformation;
        ERROR_TOOLS_CATCH( gradientMicroDeformation = getGradientMicroDeformation( ) );

        _nonLocalMicroDeformation.second = *microDeformationBase;

        for ( unsigned int i = 0; i < _dimension; i++ ){

            for ( unsigned int I = 0; I < _dimension; I++ ){

                for ( unsigned int J = 0; J < _dimension; J++ ){

                    _nonLocalMicroDeformation.second[ _dimension * i + I ]
                        += ( *gradientMicroDeformation )[ _dimension * _dimension * i + _dimension * I + J ]
                         * ( *localReferenceParticleSpacing )[ J ];

                }

            }

        }

        _nonLocalMicroDeformation.first = true;

        return;

    }

    const floatVector* aspBase::getNonLocalMicroDeformation( ){
        /*!
         * Get the non-local micro deformation tensor
         */

        if ( !_nonLocalMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setNonLocalMicroDeformation( ) );

        }

        return &_nonLocalMicroDeformation.second;

    }

    const floatVector* aspBase::getNonLocalMicroDeformationBase( ){
        /*!
         * Get the base of the non-local micro deformation tensor
         */

        if ( !_nonLocalMicroDeformationBase.first ){

            ERROR_TOOLS_CATCH( setNonLocalMicroDeformationBase( ) );

        }

        return &_nonLocalMicroDeformationBase.second;

    }

    void aspBase::setNonLocalMicroDeformationBase( ){
        /*!
         * Set the non-local micro-deformation tensor base
         */

        const floatVector *microDeformationTensor;
        ERROR_TOOLS_CATCH( microDeformationTensor = getMicroDeformation( ) );

        _nonLocalMicroDeformationBase.second = *microDeformationTensor;

        _nonLocalMicroDeformationBase.first = true;

    }

    void aspBase::setCurrentDistanceVector( ){
        /*!
         * Set the current distance vector
         */

        const floatVector* localSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( localSurfaceReferenceRelativePositionVector = getLocalSurfaceReferenceRelativePositionVector( ) );

        const floatVector* nonLocalSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( nonLocalSurfaceReferenceRelativePositionVector = getNonLocalSurfaceReferenceRelativePositionVector( ) );

        const floatVector* referenceDistanceVector;
        ERROR_TOOLS_CATCH( referenceDistanceVector = getReferenceDistanceVector( ) );

        const floatVector* localDeformationGradient;
        ERROR_TOOLS_CATCH( localDeformationGradient = getLocalDeformationGradient( ) );

        const floatVector* localMicroDeformation;
        ERROR_TOOLS_CATCH( localMicroDeformation = getLocalMicroDeformation( ) );

        const floatVector* nonLocalMicroDeformation;
        ERROR_TOOLS_CATCH( nonLocalMicroDeformation = getNonLocalMicroDeformation( ) );

        // Compute the current distance
        ERROR_TOOLS_CATCH( tractionSeparation::computeCurrentDistanceGeneral( *localSurfaceReferenceRelativePositionVector,
                                                                              *nonLocalSurfaceReferenceRelativePositionVector,
                                                                              *referenceDistanceVector,
                                                                              *localDeformationGradient,
                                                                              *localMicroDeformation,
                                                                              *nonLocalMicroDeformation,
                                                                              _currentDistanceVector.second ) );

        _currentDistanceVector.first = true;

        return;

    }

    const floatVector* aspBase::getCurrentDistanceVector( ){
        /*!
         * Get the current distance vector
         */

        if ( !_currentDistanceVector.first ){

            ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

        }

        return &_currentDistanceVector.second;

    }

    void aspBase::setLocalCurrentNormal( ){
        /*!
         * Set the current local normal
         */

        const floatVector* localReferenceNormal;
        ERROR_TOOLS_CATCH( localReferenceNormal = getLocalReferenceNormal( ) );

        const floatVector* localMicroDeformation;
        ERROR_TOOLS_CATCH( localMicroDeformation = getLocalMicroDeformation( ) );

        // Compute the current normal
        ERROR_TOOLS_CATCH( tractionSeparation::computeNansonsRelation( *localMicroDeformation, *localReferenceNormal,
                                                                       _localCurrentNormal.second ) );

        _localCurrentNormal.second /= vectorTools::l2norm( _localCurrentNormal.second );

        _localCurrentNormal.first = true;

        return;

    }

    const floatVector* aspBase::getLocalCurrentNormal( ){
        /*!
         * Get the local current normal vector
         */

        if ( !_localCurrentNormal.first ){

            ERROR_TOOLS_CATCH( setLocalCurrentNormal( ) );

        }

        return &_localCurrentNormal.second;

    }

    void aspBase::initializeSurfaceIntegrandQuantities( ){
        /*!
         * Initialize the surface integrand quantities
         */

        ERROR_TOOLS_CATCH( setLocalReferenceNormal( ) );

        ERROR_TOOLS_CATCH( setLocalSurfaceReferenceRelativePositionVector( ) );

        ERROR_TOOLS_CATCH( setNonLocalSurfaceReferenceRelativePositionVector( ) );

        ERROR_TOOLS_CATCH( setReferenceDistanceVector( ) );

        ERROR_TOOLS_CATCH( setLocalDeformationGradient( ) );

        ERROR_TOOLS_CATCH( setLocalMicroDeformation( ) );

        ERROR_TOOLS_CATCH( setNonLocalMicroDeformation( ) );

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

        ERROR_TOOLS_CATCH( setLocalCurrentNormal( ) );

        ERROR_TOOLS_CATCH( setSurfaceParameters( ) );

        return;

    }

    void aspBase::setSurfaceParameters( ){
        /*!
         * Set the surface parameters for the integrand
         */

        _surfaceParameters.first = true;

        return;
    }

    const floatVector* aspBase::getSurfaceParameters( ){
        /*!
         * Get the surface parameters
         */

        if ( !_surfaceParameters.first ){

            ERROR_TOOLS_CATCH( setSurfaceParameters( ) );

        }

        return &_surfaceParameters.second;

    }

    void aspBase::setReferenceDistanceVector( ){
        /*!
         * Set the initial distance between two particles
         */

        _referenceDistanceVector.second = floatVector( _dimension, 0. );

        _referenceDistanceVector.first = true;

        return;

    }

    const floatVector* aspBase::getReferenceDistanceVector( ){
        /*!
         * Get the initial distance between two particles
         */

        if ( !_referenceDistanceVector.first ){

            ERROR_TOOLS_CATCH( setReferenceDistanceVector( ) );

        }

        return &_referenceDistanceVector.second;

    }

    void aspBase::resetSurface( ){
        /*!
         * Reset the surface to the base state
         */

        _localReferenceRadius = std::make_pair( false, 0. );

        _nonLocalReferenceRadius = std::make_pair( false, 0. );

        _localReferenceNormal = std::make_pair( false, floatVector( 0, 0 ) );

        _localSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( 0, 0 ) );

        _nonLocalSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( 0, 0 ) );

        _referenceDistanceVector = std::make_pair( false, floatVector( 0, 0 ) );

        _localReferenceParticleSpacing = std::make_pair( false, floatVector( 0, 0 ) );

        _localDeformationGradient = std::make_pair( false, floatVector( 0, 0 ) );

        _localMicroDeformation = std::make_pair( false, floatVector( 0, 0 ) );

        _nonLocalMicroDeformation = std::make_pair( false, floatVector( 0, 0 ) );

        _currentDistanceVector = std::make_pair( false, floatVector( 0, 0 ) );

        _localCurrentNormal = std::make_pair( false, floatVector( 0, 0 ) );

        _surfaceParameters = std::make_pair( false, floatVector( 0, 0 ) );

        return;

    }

    void aspBase::computeSurfaceAdhesionEnergyDensity( floatType &surfaceAdhesionEnergyDensity ){
        /*!
         * Compute the surface adhesion energy density in the current configuration ( energy / da )
         * 
         * \param &surfaceAdhesionEnergyDensity: The surface adhesion energy density i.e., the energy bonding the local and non-local particles together
         */

        const floatVector* currentDistanceVector;
        ERROR_TOOLS_CATCH( currentDistanceVector = getCurrentDistanceVector( ) );

        const floatVector* localCurrentNormal;
        ERROR_TOOLS_CATCH( localCurrentNormal = getLocalCurrentNormal( ) );

        const floatVector* surfaceParameters;
        ERROR_TOOLS_CATCH( surfaceParameters = getSurfaceParameters( ) );

        // Decompose the traction into normal and tangential directions
        floatVector dn, dt;
        ERROR_TOOLS_CATCH_NODE_POINTER( tractionSeparation::decomposeVector( *currentDistanceVector, *localCurrentNormal, dn, dt ) );

        floatType energyDensity;
        ERROR_TOOLS_CATCH_NODE_POINTER( tractionSeparation::computeLinearTractionEnergy( dn, dt, *surfaceParameters, energyDensity ) );

        surfaceAdhesionEnergyDensity = 0.5 * energyDensity * vectorTools::l2norm( dn );

        return;

    }

    bool aspBase::pointInBoundingBox( const floatVector &point, const floatMatrix &boundingBox ){
        /*!
         * Determine if the point is inside of the bounding box and return a boolean value
         * 
         * \param &point: The point in [x, y, z, ...] space. Must be of the same length of the bounding box
         * \param &boundingBox: The bouning box in the form (point.size( ), 2 ) where for each dimension
         *    the row is of the form (lower bound, upper bound)
         */

        if ( point.size( ) != boundingBox.size( ) ){

            ERROR_TOOLS_CATCH( throw std::runtime_error( "point and boundingBox must be the same size.\n  point.size( ): " + std::to_string( point.size( ) ) + "\n  boundingBox.size( ): " + std::to_string( boundingBox.size( ) ) ) );

        }

        for ( unsigned int i = 0; i < point.size( ); i++ ){

            if ( boundingBox[ i ].size( ) != 2 ){

                ERROR_TOOLS_CATCH( throw std::runtime_error( "boundingBox row " + std::to_string( i ) + " has a length of " + std::to_string( boundingBox[ i ].size( ) ) + " and it should be of length 2" ) );

            }

            if ( ( point[ i ] < boundingBox[ i ][ 0 ] ) || ( point[ i ] > boundingBox[ i ][ 1 ] ) ){

                return false;

            }

        }

        return true;

    }

    void aspBase::idBoundingBoxContainedPoints( const floatVector &points, const floatMatrix &boundingBox, std::vector< unsigned int > &containedPoints ){
        /*!
         * Determine which of the points are in the bounding box
         * 
         * \param &points: The points in [x1, y1, z1, x2, y2, z2, ... ] format.
         * \param &boundingBox: The bounding box which must be in the form (_dimension, 2 ) where for each dimension
         *    the row is of the form (lower bound, upper bound)
         * \param &containedPoints: The points contained in the vector
         */

        containedPoints = std::vector< unsigned int >( points.size( ) );

        unsigned int numContainedPoints = 0;

        for ( unsigned int i = 0; i < points.size( ); i += _dimension ){

            if ( pointInBoundingBox( floatVector( points.begin( ) + i, points.begin( ) + i + _dimension ), boundingBox ) ){

                containedPoints[ numContainedPoints ] = i / _dimension;

                numContainedPoints++;

            }

        }

        containedPoints.resize( numContainedPoints );

        return;

    }

    void aspBase::computeSurfaceOverlapEnergyDensity( std::unordered_map< unsigned int, floatType > &surfaceOverlapEnergyDensity ){
        /*!
         * Compute the surface overlap energy density for the local particle and a given interaction
         * 
         * \param &surfaceOverlapEnergyDensities: The index is the index of the local points which may be overlapping with a neighboring particle and the values are the overlap energies. It is not a scalar because two particles could be overlapping at multiple points
         */

        surfaceOverlapEnergyDensity.clear( );

        const std::unordered_map< unsigned int, floatVector > *particlePairOverlap;
        ERROR_TOOLS_CATCH( particlePairOverlap = getParticlePairOverlap( ) );

        const floatVector *overlapParameters;
        ERROR_TOOLS_CATCH( overlapParameters = getSurfaceOverlapParameters( ) );

        for ( auto overlap = particlePairOverlap->begin( ); overlap != particlePairOverlap->end( ); overlap++ ){

            surfaceOverlapEnergyDensity.insert( { overlap->first, 0.5 * ( *overlapParameters )[ 0 ] * vectorTools::dot( overlap->second, overlap->second ) } );

        }

        return;

    }

    void aspBase::setParticlePairOverlap( ){
        /*!
         * Set the particle overlap for the current local to non-local pair
         */

        _particlePairOverlap.second.clear( );

        const floatMatrix *nonLocalBoundingBox;
        ERROR_TOOLS_CATCH( nonLocalBoundingBox = getNonLocalParticleCurrentBoundingBox( ) );

        const floatVector *localReferenceSurfacePoints;
        ERROR_TOOLS_CATCH( localReferenceSurfacePoints = getLocalReferenceSurfacePoints( ) );

        const floatVector *localCurrentSurfacePoints;
        ERROR_TOOLS_CATCH( localCurrentSurfacePoints = getLocalCurrentSurfacePoints( ) );

        const floatVector *localDeformationGradient;
        ERROR_TOOLS_CATCH( localDeformationGradient = getLocalDeformationGradient( ) );

        const floatType *nonLocalReferenceRadius;
        ERROR_TOOLS_CATCH( nonLocalReferenceRadius = getNonLocalReferenceRadius( ) );

        const floatVector *localMicroDeformation;
        ERROR_TOOLS_CATCH( localMicroDeformation = getLocalMicroDeformation( ) );

        const floatVector *localReferenceParticleSpacing;
        ERROR_TOOLS_CATCH( localReferenceParticleSpacing = getLocalReferenceParticleSpacing( ) );

        const floatVector *localGradientMicroDeformation;
        ERROR_TOOLS_CATCH( localGradientMicroDeformation = getLocalGradientMicroDeformation( ) );

        // Check which of the local points are contained in the non-local bounding box
        std::vector< unsigned int > possiblePoints;
        ERROR_TOOLS_CATCH( idBoundingBoxContainedPoints( *localCurrentSurfacePoints, *nonLocalBoundingBox, possiblePoints ) );

        const floatVector *nonLocalMicroDeformationBase;
        ERROR_TOOLS_CATCH( nonLocalMicroDeformationBase = getNonLocalMicroDeformationBase( ) );

        for ( auto p = possiblePoints.begin( ); p != possiblePoints.end( ); p++ ){

            // Compute the overlap between the local and non-local particles
            floatVector overlap;

            ERROR_TOOLS_CATCH( tractionSeparation::computeParticleOverlap( floatVector( localReferenceSurfacePoints->begin( ) + _dimension * ( *p ),
                                                                                        localReferenceSurfacePoints->begin( ) + _dimension * ( ( *p ) + 1 ) ),
                                                                           *localReferenceParticleSpacing, *nonLocalReferenceRadius, *localDeformationGradient,
                                                                           *localMicroDeformation, *nonLocalMicroDeformationBase, *localGradientMicroDeformation,
                                                                           overlap ) );

            _particlePairOverlap.second.insert( { *p, overlap } );

        }

        _particlePairOverlap.first = true;

        return;

    }

    const std::unordered_map< unsigned int, floatVector >* aspBase::getParticlePairOverlap( ){
        /*!
         * Get an unordered map of the overlap distance between the current local to non-local particle pair
         */

        if ( !_particlePairOverlap.first ){

            ERROR_TOOLS_CATCH( setParticlePairOverlap( ) );

        }

        return &_particlePairOverlap.second;

    }

    const floatType* aspBase::getSurfaceAdhesionEnergyDensity( ){
        /*!
         * Get the surface adhesion energy density
         */

        if ( !_surfaceAdhesionEnergyDensity.first ){

            ERROR_TOOLS_CATCH( setSurfaceAdhesionEnergyDensity( ) );

        }

        return &_surfaceAdhesionEnergyDensity.second;

    }

    const std::unordered_map< unsigned int, floatType >* aspBase::getSurfaceOverlapEnergyDensity( ){
        /*!
         * Get the surface overlap energy density
         */

        if ( !_surfaceOverlapEnergyDensity.first ){

            ERROR_TOOLS_CATCH( setSurfaceOverlapEnergyDensity( ) );

        }

        return &_surfaceOverlapEnergyDensity.second;

    }

    const floatVector *aspBase::getGradientMicroDeformation( ){
        /*!
         * Get the gradient of the micro-deformation
         */

        return &_gradientMicroDeformation;

    }

    const floatVector *aspBase::getMicroDeformation( ){
        /*!
         * Get the gradient of the micro-deformation
         */

        return &_microDeformation;

    }

    void aspBase::setLocalGradientMicroDeformation( ){
        /*!
         * Set the local gradient of the micro-deformation
         */

        const floatVector *gradientMicroDeformation;
        ERROR_TOOLS_CATCH( gradientMicroDeformation = getGradientMicroDeformation( ) );

        _localGradientMicroDeformation.second = *gradientMicroDeformation;

        _localGradientMicroDeformation.first = true;

    }

    const floatVector *aspBase::getLocalGradientMicroDeformation( ){
        /*!
         * Get the local gradient of the micro-deformation
         */

        if ( !_localGradientMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setLocalGradientMicroDeformation( ) );

        }

        return &_localGradientMicroDeformation.second;

    }

    void aspBase::setSurfaceAdhesionEnergyDensity( ){
        /*!
         * Set the surface adhesion energy density if required.
         */

        ERROR_TOOLS_CATCH( computeSurfaceAdhesionEnergyDensity( _surfaceAdhesionEnergyDensity.second ) );

        _surfaceAdhesionEnergyDensity.first = true;

        return;

    }

    void aspBase::setSurfaceOverlapEnergyDensity( ){
        /*!
         * Set the surface overlap energy density if required.
         */

        ERROR_TOOLS_CATCH( computeSurfaceOverlapEnergyDensity( _surfaceOverlapEnergyDensity.second ) );

        _surfaceOverlapEnergyDensity.first = true;

        return;

    }

    void aspBase::setSurfaceOverlapParameters( ){
        /*!
         * Set the surface overlap parameters
         */

        _surfaceOverlapParameters.second = { 1 };

        _surfaceOverlapParameters.first = true;

        return;

    }

    const floatVector *aspBase::getSurfaceOverlapParameters( ){
        /*!
         * Get the surface overlap parameters
         */

        if ( !_surfaceOverlapParameters.first ){

            ERROR_TOOLS_CATCH( setSurfaceOverlapParameters( ) );

        }

        return &_surfaceOverlapParameters.second;

    }

    void aspBase::setLocalReferenceSurfacePoints( ){
        /*!
         * Set the collection of points which are on the surface of the local particle in the reference configuration
         */

        const floatVector *unitSpherePoints;
        ERROR_TOOLS_CATCH( unitSpherePoints = getUnitSpherePoints( ) );

        const floatType *localReferenceRadius;
        ERROR_TOOLS_CATCH( localReferenceRadius = getLocalReferenceRadius( ) );

        _localReferenceSurfacePoints.second = ( *localReferenceRadius ) * ( *unitSpherePoints );

        _localReferenceSurfacePoints.first = true;

        return;

    }

    const floatVector* aspBase::getLocalReferenceSurfacePoints( ){
        /*!
         * Get the local reference surface points
         */

        if ( !_localReferenceSurfacePoints.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceSurfacePoints( ) );

        }

        return &_localReferenceSurfacePoints.second;

    }

    void aspBase::setNonLocalReferenceSurfacePoints( ){
        /*!
         * Set the collection of points which are on the surface of the non-local particle in the reference configuration
         */

        const floatVector *unitSpherePoints;
        ERROR_TOOLS_CATCH( unitSpherePoints = getUnitSpherePoints( ) );

        const floatType *nonLocalReferenceRadius;
        ERROR_TOOLS_CATCH( nonLocalReferenceRadius = getNonLocalReferenceRadius( ) );

        _nonLocalReferenceSurfacePoints.second = ( *nonLocalReferenceRadius ) * ( *unitSpherePoints );

        _nonLocalReferenceSurfacePoints.first = true;

        return;

    }

    const floatVector* aspBase::getNonLocalReferenceSurfacePoints( ){
        /*!
         * Get the non-local reference surface points
         */

        if ( !_nonLocalReferenceSurfacePoints.first ){

            ERROR_TOOLS_CATCH( setNonLocalReferenceSurfacePoints( ) );

        }

        return &_nonLocalReferenceSurfacePoints.second;

    }

    void aspBase::setLocalCurrentSurfacePoints( ){
        /*!
         * Set the collection of points on the surface of the local particle in the current configuration
         */

        const floatVector *localReferenceSurfacePoints;
        ERROR_TOOLS_CATCH( localReferenceSurfacePoints = getLocalReferenceSurfacePoints( ) );

        const floatVector *localMicroDeformation;
        ERROR_TOOLS_CATCH( localMicroDeformation = getLocalMicroDeformation( ) );

        ERROR_TOOLS_CATCH( _localCurrentSurfacePoints.second = vectorTools::matrixMultiply( *localReferenceSurfacePoints, *localMicroDeformation,
                                                                                            localReferenceSurfacePoints->size( ) / _dimension,
                                                                                            _dimension, _dimension, _dimension, false, true ) );

        _localCurrentSurfacePoints.first = true;

    }

    const floatVector* aspBase::getLocalCurrentSurfacePoints( ){
        /*!
         * Get the collection of points on the surface of the local particle in the current configuration
         */

        if ( !_localCurrentSurfacePoints.first ){

            ERROR_TOOLS_CATCH( setLocalCurrentSurfacePoints( ) );

        }

        return &_localCurrentSurfacePoints.second;

    }

    void aspBase::setNonLocalCurrentSurfacePoints( ){
        /*!
         * Set the collection of points on the surface of the non-local particle in the current configuration
         */

        const floatVector *nonLocalReferenceSurfacePoints;
        ERROR_TOOLS_CATCH( nonLocalReferenceSurfacePoints = getNonLocalReferenceSurfacePoints( ) );

        const floatVector *nonLocalMicroDeformation;
        ERROR_TOOLS_CATCH( nonLocalMicroDeformation = getNonLocalMicroDeformation( ) );

        ERROR_TOOLS_CATCH( _nonLocalCurrentSurfacePoints.second = vectorTools::matrixMultiply( *nonLocalReferenceSurfacePoints, *nonLocalMicroDeformation,
                                                                                               nonLocalReferenceSurfacePoints->size( ) / _dimension,
                                                                                               _dimension, _dimension, _dimension, false, true ) );

        _nonLocalCurrentSurfacePoints.first = true;

    }

    const floatVector* aspBase::getNonLocalCurrentSurfacePoints( ){
        /*!
         * Get the collection of points on the surface of the non-local particle in the current configuration
         */

        if ( !_nonLocalCurrentSurfacePoints.first ){

            ERROR_TOOLS_CATCH( setNonLocalCurrentSurfacePoints( ) );

        }

        return &_nonLocalCurrentSurfacePoints.second;

    }

    void aspBase::setLocalParticleCurrentBoundingBox( ){
        /*!
         * Set the bounding box of the local particle (dimension, (lower bound, upper bound))
         */

        const floatVector *localCurrentSurfacePoints;
        ERROR_TOOLS_CATCH( localCurrentSurfacePoints = getLocalCurrentSurfacePoints( ) ); 

        ERROR_TOOLS_CATCH( formBoundingBox( *localCurrentSurfacePoints, _localParticleCurrentBoundingBox.second ) );

        _localParticleCurrentBoundingBox.first = true;

        return;

    }

    void aspBase::formBoundingBox( const floatVector &points, floatMatrix &boundingBox ){
        /*!
         * Form a bounding box from the provided points
         * 
         * \param &points: The points to form the bounding box in [x1, y1, z1, x2, y2, z2, ...] format
         * \param &boundingBox: The resulting bounding box
         */

        if ( ( points.size( ) % _dimension ) > 0 ){

            ERROR_TOOLS_CATCH( throw std::runtime_error( "The length of points must be an integer multiple of the dimension.\n  points.size( ) = " + std::to_string( points.size( ) ) + "\n  dimension: " + std::to_string( _dimension ) ) );

        }

        if ( points.size( ) < _dimension ){

            ERROR_TOOLS_CATCH( throw std::runtime_error( "The points need at least one point" ) );

        }

        boundingBox = { { points[ 0 ], points[ 0 ] },
                        { points[ 1 ], points[ 1 ] },
                        { points[ 2 ], points[ 2 ] } };

        for ( unsigned int i = _dimension; i < points.size( ); i += _dimension ){

            for ( unsigned int j = 0; j < _dimension; j++ ){

                boundingBox[ j ][ 0 ] = std::fmin( boundingBox[ j ][ 0 ], points[ i + j ] );
                boundingBox[ j ][ 1 ] = std::fmax( boundingBox[ j ][ 1 ], points[ i + j ] );

            }

        }

        return;

    }

    const floatMatrix* aspBase::getLocalParticleCurrentBoundingBox( ){
        /*!
         * Get the local particle's bounding box
         */

        if ( !_localParticleCurrentBoundingBox.first ){

            ERROR_TOOLS_CATCH( setLocalParticleCurrentBoundingBox( ) );

        }

        return &_localParticleCurrentBoundingBox.second;

    }

    void aspBase::setNonLocalParticleCurrentBoundingBox( ){
        /*!
         * Set the bounding box of the non-local particle (dimension, (lower bound, upper bound))
         */

        const floatVector *nonLocalCurrentSurfacePoints;
        ERROR_TOOLS_CATCH( nonLocalCurrentSurfacePoints = getNonLocalCurrentSurfacePoints( ) ); 

        ERROR_TOOLS_CATCH( formBoundingBox( *nonLocalCurrentSurfacePoints, _nonLocalParticleCurrentBoundingBox.second ) );

        _nonLocalParticleCurrentBoundingBox.first = true;

        return;

    }

    void aspBase::setSurfaceAdhesionTraction( ){
        /*!
         * Set the surface adhesion traction
         */

        ERROR_TOOLS_CATCH( computeSurfaceAdhesionTraction( _surfaceAdhesionTraction.second ) );

        _surfaceAdhesionTraction.first = true;

        return;

    }

    const floatVector* aspBase::getSurfaceAdhesionTraction( ){
        /*!
         * Get the surface adhesion traction
         */

        if ( !_surfaceAdhesionTraction.first ){

            ERROR_TOOLS_CATCH( setSurfaceAdhesionTraction( ) );

        }

        return &_surfaceAdhesionTraction.second;

    }

    void aspBase::computeSurfaceAdhesionTraction( floatVector &surfaceAdhesionTraction ){
        /*!
         * Compute the surface adhesion traction
         * 
         * \param &surfaceAdhesionTraction: The current adhesion traction between the local and non-local particles
         */

        const floatVector* currentDistanceVector;
        ERROR_TOOLS_CATCH( currentDistanceVector = getCurrentDistanceVector( ) );

        const floatVector* localCurrentNormal;
        ERROR_TOOLS_CATCH( localCurrentNormal = getLocalCurrentNormal( ) );

        const floatVector* surfaceParameters;
        ERROR_TOOLS_CATCH( surfaceParameters = getSurfaceParameters( ) );

        // Decompose the traction into normal and tangential directions
        floatVector dn, dt;
        ERROR_TOOLS_CATCH_NODE_POINTER( tractionSeparation::decomposeVector( *currentDistanceVector, *localCurrentNormal, dn, dt ) );

        ERROR_TOOLS_CATCH( tractionSeparation::computeLinearTraction( dn, dt, *surfaceParameters, surfaceAdhesionTraction ) );

        return;

    }

    void aspBase::computeSurfaceOverlapTraction( std::unordered_map< unsigned int, floatVector > &surfaceOverlapTraction ){
        /*!
         * Compute the surface overlap traction
         * 
         * \param &surfaceOverlapTraction: An unordered map mapping from points on the local particle to the amount they overlap the non-local particle
         */

        surfaceOverlapTraction.clear( );

        const std::unordered_map< unsigned int, floatVector > *particlePairOverlap;
        ERROR_TOOLS_CATCH( particlePairOverlap = getParticlePairOverlap( ) );

        const floatVector *overlapParameters;
        ERROR_TOOLS_CATCH( overlapParameters = getSurfaceOverlapParameters( ) );

        for ( auto overlap = particlePairOverlap->begin( ); overlap != particlePairOverlap->end( ); overlap++ ){

            surfaceOverlapTraction.insert( { overlap->first, ( *overlapParameters )[ 0 ] * overlap->second } );

        }

        return;

    }

    void aspBase::setSurfaceOverlapTraction( ){
        /*!
         * Set the surface overlap traction
         */

        ERROR_TOOLS_CATCH( computeSurfaceOverlapTraction( _surfaceOverlapTraction.second ) );

        _surfaceOverlapTraction.first = true;

        return;

    }

    const std::unordered_map< unsigned int, floatVector >* aspBase::getSurfaceOverlapTraction( ){
        /*!
         * Get the surface overlap traction
         */

        if ( !_surfaceOverlapTraction.first ){

            ERROR_TOOLS_CATCH( setSurfaceOverlapTraction( ) );

        }

        return &_surfaceOverlapTraction.second;

    }

    const floatMatrix* aspBase::getNonLocalParticleCurrentBoundingBox( ){
        /*!
         * Get the local particle's bounding box
         */

        if ( !_nonLocalParticleCurrentBoundingBox.first ){

            ERROR_TOOLS_CATCH( setNonLocalParticleCurrentBoundingBox( ) );

        }

        return &_nonLocalParticleCurrentBoundingBox.second;

    }

    /** \brief Define required number of Abaqus material constants for the Abaqus interface. */
    const int nStateVariables = 2;

    /** \brief Define required number of Abaqus material constants for the Abaqus interface. */
    const int nMaterialParameters = 2;

    /// Say hello
    /// @param message The message to print
    errorOut sayHello( std::string message ) {
        if ( message.compare( "George" ) == 0 ){
            errorOut result = new errorNode( __func__, "ERROR: George is a wolf in sheep's clothing!");
            return result;
        }
        std::cout << "Hello " << message << std::endl;
        return NULL;
    }

    errorOut dummyMaterialModel( floatVector &stress,             floatVector &statev,        floatMatrix &ddsdde,       floatType &SSE,            floatType &SPD,
                                 floatType &SCD,                  floatType &RPL,             floatVector &ddsddt,       floatVector &drplde,       floatType &DRPLDT,
                                 const floatVector &strain,       const floatVector &dstrain, const floatVector &time,   const floatType &DTIME,    const floatType &TEMP,
                                 const floatType &DTEMP,          const floatVector &predef,  const floatVector &dpred,  const std::string &cmname, const int &NDI,
                                 const int &NSHR,                 const int &NTENS,           const int &NSTATV,         const floatVector &props,  const int &NPROPS,
                                 const floatVector &coords,       const floatMatrix &drot,    floatType &PNEWDT,         const floatType &CELENT,   const floatMatrix &dfgrd0,
                                 const floatMatrix &dfgrd1,       const int &NOEL,            const int &NPT,            const int &LAYER,          const int &KSPT,
                                 const std::vector< int > &jstep, const int &KINC ){
        /*!
         * A template Abaqus c++ UMAT using c++ STL types. Variables in all caps reference ABAQUS FORTRAN
         * memory directly. Variables in lower case are native c++ type conversions stored separately from the original
         * ABAQUS FORTRAN memory.
         */

        //Call functions of constitutive model to do things
        errorOut error = sayHello( "Abaqus" );

        //Error handling
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when calling sayHello" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    void abaqusInterface( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                          double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                          const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                          const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                          const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                          const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                          const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                          const int *JSTEP,     const int &KINC ){
        /*!
         * A template Abaqus UMAT c++ interface that performs Fortran to C++ type conversions, calculates the material
         * model's expected input, handles tensor shape changes, and calls a c++ material model.
         */

        //Initialize error return codes
        errorOut error = NULL;

        //Provide a variable string message for error nodes
        std::ostringstream message;

        //Map FORTRAN UMAT variables to C++ types as necessary. Use case sensitivity to distinguish.
        //TODO: Decide if case sensitive variable names is a terrible idea or not
        //Vectors can be created directly with pointer arithmetic
        std::vector< double > stress( STRESS, STRESS + NTENS );
        std::vector< double > statev( STATEV, STATEV + NSTATV );
        std::vector< double > ddsddt( DDSDDT, DDSDDT + NTENS );
        std::vector< double > drplde( DRPLDE, DRPLDE + NTENS );
        const std::vector< double > strain( STRAN, STRAN + NTENS );
        const std::vector< double > dstrain( DSTRAN, DSTRAN + NTENS );
        const std::vector< double > time( TIME, TIME + 2 );
        const std::vector< double > predef( PREDEF, PREDEF + 1 );
        const std::vector< double > dpred( DPRED, DPRED + 1 );
        const std::string cmname( abaqusTools::FtoCString( 80, CMNAME ) );
        const std::vector< double > props( PROPS, PROPS + NPROPS );
        const std::vector< double > coords( COORDS, COORDS + spatialDimensions );
        const std::vector< int > jstep( JSTEP, JSTEP + 4 );
        //Fortran two-dimensional arrays require careful column to row major conversions to c++ types
        std::vector< std::vector< double > > ddsdde = abaqusTools::columnToRowMajor( DDSDDE, NTENS, NTENS );
        const std::vector< std::vector< double > > drot = abaqusTools::columnToRowMajor( DROT, spatialDimensions, spatialDimensions );
        const std::vector< std::vector< double > > dfgrd0 = abaqusTools::columnToRowMajor( DFGRD0, spatialDimensions, spatialDimensions );
        const std::vector< std::vector< double > > dfgrd1 = abaqusTools::columnToRowMajor( DFGRD1, spatialDimensions, spatialDimensions );

        //Verify number of state variables against asp expectations
        if ( statev.size( ) != nStateVariables ){
            message.clear();
            message << "ERROR:" << __FILENAME__ << "." << __func__ << ": The asp Abaqus interface requires exactly "
                << nStateVariables << " state variables. Found " << statev.size( ) << ".";
            throw std::runtime_error( message.str( ) );
        }

        //Verify number of material parameters against asp expectations
        if ( props.size( ) != nMaterialParameters ){
            message.clear();
            message << "ERROR:" << __FILENAME__ << "." << __func__ << ": The asp Abaqus interface requires exactly "
                << nMaterialParameters << " material constants. Found " << props.size( ) << ".";
            throw std::runtime_error( message.str( ) );
        }

        //Call the constitutive model c++ interface
        if ( KINC == 1 && NOEL == 1 && NPT == 1 ){
            error = dummyMaterialModel( stress, statev,  ddsdde, SSE,    SPD,
                                        SCD,    RPL,     ddsddt, drplde, DRPLDT,
                                        strain, dstrain, time,   DTIME,  TEMP,
                                        DTEMP,  predef,  dpred,  cmname, NDI,
                                        NSHR,   NTENS,   NSTATV, props,  NPROPS,
                                        coords, drot,    PNEWDT, CELENT, dfgrd0,
                                        dfgrd1, NOEL,    NPT,    LAYER,  KSPT,
                                        jstep,  KINC );
        }

        //Error handling
        if ( error ){
            message.clear();
            message << "ERROR:" << __FILENAME__ << "." << __func__ << ": Error when calling dummyMaterialModel.";
            errorOut result = new errorNode( __func__, message.str( ) );
            result->addNext( error );
            error->print( true );
            //If an error was thrown, but the ratio of new/current time increment is not updated, it was a fatal error.
            if ( vectorTools::fuzzyEquals( PNEWDT, 1. ) ){
                throw std::runtime_error( message.str( ) );
            }
        }

        //Re-pack C++ objects into FORTRAN memory to return values to Abaqus
        //Scalars were passed by reference and will update correctly
        //Vectors don't require row/column major considerations, but do require re-packing to the Fortran pointer
        abaqusTools::rowToColumnMajor( STRESS, stress, 1, NTENS );
        abaqusTools::rowToColumnMajor( DDSDDT, ddsddt, 1, NTENS );
        abaqusTools::rowToColumnMajor( DRPLDE, drplde, 1, NTENS );
        abaqusTools::rowToColumnMajor( STATEV, statev, 1, NSTATV );
        //Arrays require vector of vector to column major conversion
        abaqusTools::rowToColumnMajor( DDSDDE, ddsdde, NTENS, NTENS );

    }

}
