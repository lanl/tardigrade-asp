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
        const unsigned int *dim = getDimension( );

        _localReferenceRadius = dataStorage( false, 0. );

        _nonLocalReferenceRadius = dataStorage( false, 0. );

        _unitSpherePoints = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _localReferenceNormal = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _localSurfaceReferenceRelativePositionVector = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _nonLocalSurfaceReferenceRelativePositionVector = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _referenceDistanceVector = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _localReferenceParticleSpacing = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _localDeformationGradient = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _localMicroDeformation = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _nonLocalMicroDeformation = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _currentDistanceVector = dataStorage( false, floatVector( ( *dim ), 0 ) );

        _localCurrentNormal = dataStorage( false, floatVector( ( *dim ), 0 ) );

    }

    void aspBase::computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
                                                     const floatVector &currentMicroDeformation, const floatVector &previousMicroDeformation,
                                                     const floatType &currentTemperature, const floatType &previousTemperature,
                                                     const floatVector &previousStateVariables,
                                                     const floatVector &parameters,
                                                     floatType &energyDensity, floatVector &cauchyStress, floatVector &stateVariables ){
        /*!
         * Default function for computing the local particle's energy in the current configuration
         * Defaults to a linear-elastic solid with parameters defined in the reference configuration.
         * 
         * \param &previousTime: The previous value of the time
         * \param &deltaTime: The change in time
         * \param &currentMicroDeformation: The current value of the micro-deformation tensor \f$\chi\f$
         * \param &previousMicroDeformation: The previous value of the micro-deformation tensor \f$\chi\f$
         * \param &currentTemperature: The current value of the temperature
         * \param &previousTemperature: The previous value of the temperature
         * \param &parameters: The parameters for the model
         * \param &energyDensity: The current value of the Helmholtz free energy density in the current configuration
         * \param &cauchyStress: The current value of the Cauchy stress
         * \param &stateVariables: The updated values of the state variables
         */

        ERROR_TOOLS_CATCH_NODE_POINTER( stressTools::linearElasticity::evaluateEnergy( currentMicroDeformation, parameters, energyDensity, cauchyStress ) );

        return;

    }

    void aspBase::computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
                                                     const floatVector &currentMicroDeformation, const floatVector &previousMicroDeformation,
                                                     const floatType &currentTemperature, const floatType &previousTemperature,
                                                     const floatVector &previousStateVariables,
                                                     const floatVector &parameters,
                                                     floatType &energyDensity, floatVector &cauchyStress, floatVector &stateVariables,
                                                     floatType &logProbabilityRatio ){
        /*!
         * Default function for computing the local particle's energy density in the current configuration
         * Defaults to a linear-elastic solid with parameters defined in the reference configuration.
         * 
         * \param &previousTime: The previous value of the time
         * \param &deltaTime: The change in time
         * \param &currentMicroDeformation: The current value of the micro-deformation tensor \f$\chi\f$
         * \param &previousMicroDeformation: The previous value of the micro-deformation tensor \f$\chi\f$
         * \param &currentTemperature: The current value of the temperature
         * \param &previousTemperature: The previous value of the temperature
         * \param &parameters: The parameters for the model
         * \param &energyDensity: The current value of the Helmholtz free energy density in the current configuration
         * \param &cauchyStress: The current value of the Cauchy stress
         * \param &stateVariables: The updated values of the state variables
         * \param &logProbabilityRatio: The natural logarithm of the probability ratio between the new
         *     probability density and the previous value.
         */

        ERROR_TOOLS_CATCH_NODE_POINTER( stressTools::linearElasticity::evaluateEnergy( currentMicroDeformation, parameters, energyDensity, cauchyStress ) );

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

    const unsigned int* aspBase::getDimension( ){
        /*!
         * Get the spacial dimension
         */

        return &_dimension;

    }

    const floatType* aspBase::getRadius( ){
        /*!
         * Get the radius
         */

        return &_radius;

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

        floatVector localReferenceNormal;

        ERROR_TOOLS_CATCH( getLocalReferenceNormal( *getLocalSurfaceNodeIndex( ), localReferenceNormal ) );

        setLocalReferenceNormal( localReferenceNormal );

        return;

    }

    void aspBase::setLocalReferenceNormal( const floatVector &value ){
        /*!
         * Set the local reference normal vector
         * 
         * \param &value: The local reference normal vector
         * 
         */

        _localReferenceNormal.second = value;

        _localReferenceNormal.first = true;

        addSurfacePointData( &_localReferenceNormal );

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

        setLocalSurfaceReferenceRelativePositionVector( ( *localReferenceRadius ) * ( *localReferenceNormal ) );

        return;

    }

    void aspBase::setLocalSurfaceReferenceRelativePositionVector( const floatVector &value ){
        /*!
         * Set the local surface relative position vector
         * 
         * \param &value: The value of the local surface relative position vector
         */

        _localSurfaceReferenceRelativePositionVector.second = value;

        _localSurfaceReferenceRelativePositionVector.first = true;

        addSurfacePointData( &_localSurfaceReferenceRelativePositionVector );

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

        setNonLocalSurfaceReferenceRelativePositionVector( -( *nonLocalReferenceRadius ) * ( *localReferenceNormal ) );

        return;

    }

    void aspBase::setNonLocalSurfaceReferenceRelativePositionVector( const floatVector &value ){
        /*!
         * Set the non-local relative position vector
         * 
         * \param &value: The value of the non-local surface reference relative position vector
         */

        _nonLocalSurfaceReferenceRelativePositionVector.second = value;

        _nonLocalSurfaceReferenceRelativePositionVector.first = true;

        addInteractionPairData( &_nonLocalSurfaceReferenceRelativePositionVector );

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

        setLocalReferenceRadius( *getRadius( ) );

        return;

    }

    void aspBase::setLocalReferenceRadius( const floatType &value ){
        /*!
         * Set the local reference radius
         * 
         * \param &value: The value of the local reference radius
         */

        _localReferenceRadius.second = value;

        _localReferenceRadius.first = true;

        addLocalParticleData( &_localReferenceRadius );

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

        setNonLocalReferenceRadius( *getRadius( ) );

        return;

    }

    void aspBase::setNonLocalReferenceRadius( const floatType &value ){
        /*!
         * Set the non-local reference radius
         * 
         * \param &value: The value of the non-local reference radius
         */

        _nonLocalReferenceRadius.second = value;

        _nonLocalReferenceRadius.first = true;

        addInteractionPairData( &_nonLocalReferenceRadius );

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

    const floatVector* aspBase::getDeformationGradient( ){
        /*!
         * Get the deformation gradient
         */

        return &_deformationGradient;

    }

    const floatVector* aspBase::getPreviousDeformationGradient( ){
        /*!
         * Get the previous deformation gradient
         */

        return &_previousDeformationGradient;

    }

    const floatVector* aspBase::getPreviousMicroDeformation( ){
        /*!
         * Get the previous value of the micro deformation
         */

        return &_previousMicroDeformation;

    }

    void aspBase::setLocalDeformationGradient( ){
        /*!
         * Set the local deformation gradient
         */

        const floatVector* deformationGradient;
        ERROR_TOOLS_CATCH( deformationGradient = getDeformationGradient( ) );

        setLocalDeformationGradient( *deformationGradient );

        return;

    }

    void aspBase::setLocalDeformationGradient( const floatVector &value ){
        /*!
         * Set the local deformation gradient
         * 
         * \param &value: The value of the local deformation gradient
         */

        _localDeformationGradient.second = value;

        _localDeformationGradient.first = true;

        addLocalParticleData( &_localDeformationGradient );

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

    void aspBase::setPreviousLocalDeformationGradient( ){
        /*!
         * Set the local previous deformation gradient
         */

        const floatVector* previousDeformationGradient;
        ERROR_TOOLS_CATCH( previousDeformationGradient = getPreviousDeformationGradient( ) );

        setPreviousLocalDeformationGradient( *previousDeformationGradient );

        return;

    }

    void aspBase::setPreviousLocalDeformationGradient( const floatVector &value ){
        /*!
         * Set the local previous deformation gradient
         * 
         * \param &value: The value of the local previous deformation gradient
         */

        _previousLocalDeformationGradient.second = value;

        _previousLocalDeformationGradient.first = true;

        addLocalParticleData( &_previousLocalDeformationGradient );

        return;

    }

    const floatVector* aspBase::getPreviousLocalDeformationGradient( ){
        /*!
         * Get the local previous deformation gradient
         */

        if ( !_previousLocalDeformationGradient.first ){

            ERROR_TOOLS_CATCH( setPreviousLocalDeformationGradient( ) );

        }

        return &_previousLocalDeformationGradient.second;

    }


    void aspBase::setLocalMicroDeformation( ){
        /*!
         * Set the local micro deformation
         */

        const floatVector* microDeformation;
        ERROR_TOOLS_CATCH( microDeformation = getMicroDeformation( ) );

        setLocalMicroDeformation( *microDeformation );

        return;

    }

    void aspBase::setLocalMicroDeformation( const floatVector &value ){
        /*!
         * Set the local micro deformation
         * 
         * \param &value: The value of the local micro deformation
         */

        _localMicroDeformation.second = value;

        _localMicroDeformation.first = true;

        addLocalParticleData( &_localMicroDeformation );

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

    void aspBase::setPreviousLocalMicroDeformation( ){
        /*!
         * Set the local previous micro deformation
         */

        const floatVector* previousMicroDeformation;
        ERROR_TOOLS_CATCH( previousMicroDeformation = getPreviousMicroDeformation( ) );

        setPreviousLocalMicroDeformation( *previousMicroDeformation );

        return;

    }

    void aspBase::setPreviousLocalMicroDeformation( const floatVector &value ){
        /*!
         * Set the local previous micro deformation
         * 
         * \param &value: The previous local micro deformation
         */

        _previousLocalMicroDeformation.second = value;

        _previousLocalMicroDeformation.first = true;

        addLocalParticleData( &_previousLocalMicroDeformation );

        return;

    }

    const floatVector* aspBase::getPreviousLocalMicroDeformation( ){
        /*!
         * Get the local previous micro deformation
         */

        if ( !_previousLocalMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setPreviousLocalMicroDeformation( ) );

        }

        return &_previousLocalMicroDeformation.second;

    }

    void aspBase::setLocalReferenceParticleSpacingVector( ){
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

        floatVector localReferenceParticleSpacing;

        localReferenceParticleSpacing = ( *localSurfaceReferenceRelativePositionVector )
                                      + ( *referenceDistanceVector )
                                      - ( *nonLocalSurfaceReferenceRelativePositionVector );

        setLocalReferenceParticleSpacingVector( localReferenceParticleSpacing );

        return;

    }

    void aspBase::setLocalReferenceParticleSpacingVector( const floatVector &value ){
        /*!
         * Set the local particle spacing
         * 
         * \param &value: The local particle spacing
         */

        _localReferenceParticleSpacing.second = value;

        _localReferenceParticleSpacing.first = true;

        addInteractionPairData( &_localReferenceParticleSpacing );

        return;

    }

    const floatVector* aspBase::getLocalReferenceParticleSpacingVector( ){
        /*!
         * Get the local reference particle spacing
         */

        if ( !_localReferenceParticleSpacing.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceParticleSpacingVector( ) );

        }

        return &_localReferenceParticleSpacing.second;

    }

    void aspBase::setNonLocalMicroDeformation( ){
        /*!
         * Set the non-local micro deformation given
         * 
         * \f$ \chi_{iI}^{\text{NL}} = \chi_{iI}^{\text{NL base}} + \chi_{iI,J} dX_J\f$
         * 
         * where \f$ \chi_{iI}^{\text{NL}} \f$ is the non-local micro-deformation, \f$\chi_{iI}^{\text{NL base}}\f$
         * is the non-local micro-deformation base, \f$\chi_{iI,J}\f$ is the gradient of the micro-deformation w.r.t.
         * the reference spacial coordinate, and \f$dX_J\f$ is the center-to-center spacing of the two particles defined as
         * 
         * \f$dX_I = \Xi_I + D_I - \Xi_I^{\text{NL}} \f$
         * 
         * where \f$\Xi_I\f$ is the local relative position vector in the reference configuration,
         * \f$D_I\f$ is the reference distance, and \f$\Xi_I^{\text{NL}}\f$ is the
         * non-local relative position vector in the reference configuration
         */

        const unsigned int *dim = getDimension( );

        const floatVector* localReferenceParticleSpacing;
        ERROR_TOOLS_CATCH( localReferenceParticleSpacing = getLocalReferenceParticleSpacingVector( ) );

        const floatVector* microDeformationBase;
        ERROR_TOOLS_CATCH( microDeformationBase = getNonLocalMicroDeformationBase( ) );

        const floatVector* gradientMicroDeformation;
        ERROR_TOOLS_CATCH( gradientMicroDeformation = getGradientMicroDeformation( ) );

        floatVector nonLocalMicroDeformation = *microDeformationBase;

        for ( unsigned int i = 0; i < ( *dim ); i++ ){

            for ( unsigned int I = 0; I < ( *dim ); I++ ){

                for ( unsigned int J = 0; J < ( *dim ); J++ ){

                    nonLocalMicroDeformation[ ( *dim ) * i + I ]
                        += ( *gradientMicroDeformation )[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ]
                         * ( *localReferenceParticleSpacing )[ J ];

                }

            }

        }

        setNonLocalMicroDeformation( nonLocalMicroDeformation );

        return;

    }

    void aspBase::setdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the local relative position
         * vector in the reference configuration given
         * 
         * \f$ \chi_{iI}^{\text{NL}} = \chi_{iI}^{\text{NL base}} + \chi_{iI,J} dX_J\f$
         * 
         * where \f$ \chi_{iI}^{\text{NL}} \f$ is the non-local micro-deformation, \f$\chi_{iI}^{\text{NL base}}\f$
         * is the non-local micro-deformation base, \f$\chi_{iI,J}\f$ is the gradient of the micro-deformation w.r.t.
         * the reference spacial coordinate, and \f$dX_J\f$ is the center-to-center spacing of the two particles defined as
         * 
         * \f$dX_I = \Xi_I + D_I - \Xi_I^{\text{NL}} \f$
         * 
         * where \f$\Xi_I\f$ is the local relative position vector in the reference configuration,
         * \f$D_I\f$ is the reference distance, and \f$\Xi_I^{\text{NL}}\f$ is the
         * non-local relative position vector in the reference configuration
         */

        const unsigned int *dim = getDimension( );

        setdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( vectorTools::inflate( *getGradientMicroDeformation( ), ( *dim ) * ( *dim ), ( *dim ) ) );

    }

    void aspBase::setdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( const floatMatrix &value ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the local relative position
         * vector in the reference configuration
         * 
         * \param &value: The value of the gradient of the non-local micro-deformation w.r.t. the local reference relative
         *     position vector
         */

        _dNonLocalMicroDeformationdLocalReferenceRelativePositionVector.second = value;

        _dNonLocalMicroDeformationdLocalReferenceRelativePositionVector.first = true;

        addInteractionPairData( &_dNonLocalMicroDeformationdLocalReferenceRelativePositionVector );

    }

    const floatMatrix* aspBase::getdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ){
        /*!
         * Get the gradient of the non-local micro-deformation w.r.t. the local reference relative position vector
         */

        if ( !_dNonLocalMicroDeformationdLocalReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ) );

        }

        return &_dNonLocalMicroDeformationdLocalReferenceRelativePositionVector.second;

    }

    void aspBase::setdNonLocalMicroDeformationdLocalReferenceDistanceVector( ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the local relative position
         * vector in the reference configuration given
         * 
         * \f$ \chi_{iI}^{\text{NL}} = \chi_{iI}^{\text{NL base}} + \chi_{iI,J} dX_J\f$
         * 
         * where \f$ \chi_{iI}^{\text{NL}} \f$ is the non-local micro-deformation, \f$\chi_{iI}^{\text{NL base}}\f$
         * is the non-local micro-deformation base, \f$\chi_{iI,J}\f$ is the gradient of the micro-deformation w.r.t.
         * the reference spacial coordinate, and \f$dX_J\f$ is the center-to-center spacing of the two particles defined as
         * 
         * \f$dX_I = \Xi_I + D_I - \Xi_I^{\text{NL}} \f$
         * 
         * where \f$\Xi_I\f$ is the local relative position vector in the reference configuration,
         * \f$D_I\f$ is the reference distance, and \f$\Xi_I^{\text{NL}}\f$ is the
         * non-local relative position vector in the reference configuration
         */

        setdNonLocalMicroDeformationdLocalReferenceDistanceVector( *getdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ) );

    }

    void aspBase::setdNonLocalMicroDeformationdLocalReferenceDistanceVector( const floatMatrix &value ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the local distance
         * vector in the reference configuration
         * 
         * \param &value: The value of the gradient of the non-local micro-deformation w.r.t. the relative
         *     spacing vector in the reference configuration
         */

        _dNonLocalMicroDeformationdLocalReferenceDistanceVector.second = value;

        _dNonLocalMicroDeformationdLocalReferenceDistanceVector.first = true;

        addInteractionPairData( &_dNonLocalMicroDeformationdLocalReferenceDistanceVector );

    }

    const floatMatrix* aspBase::getdNonLocalMicroDeformationdLocalReferenceDistanceVector( ){
        /*!
         * Get the gradient of the non-local micro-deformation w.r.t. the local reference distance vector
         */

        if ( !_dNonLocalMicroDeformationdLocalReferenceDistanceVector.first ){

            ERROR_TOOLS_CATCH( setdNonLocalMicroDeformationdLocalReferenceDistanceVector( ) );

        }

        return &_dNonLocalMicroDeformationdLocalReferenceDistanceVector.second;

    }

    void aspBase::setdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the non-local relative position
         * vector in the reference configuration given
         * 
         * \f$ \chi_{iI}^{\text{NL}} = \chi_{iI}^{\text{NL base}} + \chi_{iI,J} dX_J\f$
         * 
         * where \f$ \chi_{iI}^{\text{NL}} \f$ is the non-local micro-deformation, \f$\chi_{iI}^{\text{NL base}}\f$
         * is the non-local micro-deformation base, \f$\chi_{iI,J}\f$ is the gradient of the micro-deformation w.r.t.
         * the reference spacial coordinate, and \f$dX_J\f$ is the center-to-center spacing of the two particles defined as
         * 
         * \f$dX_I = \Xi_I + D_I - \Xi_I^{\text{NL}} \f$
         * 
         * where \f$\Xi_I\f$ is the local relative position vector in the reference configuration,
         * \f$D_I\f$ is the reference distance, and \f$\Xi_I^{\text{NL}}\f$ is the
         * non-local relative position vector in the reference configuration
         */

        setdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( -( *getdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ) ) );

    }

    void aspBase::setdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( const floatMatrix &value ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the non-local relative position
         * vector in the reference configuration
         * 
         * \param &value: The value of the gradient of the non-local micro-deformation w.r.t. the non-local reference relative
         *     position vector
         */

        _dNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector.second = value;

        _dNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector.first = true;

        addInteractionPairData( &_dNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector );

    }

    const floatMatrix* aspBase::getdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( ){
        /*!
         * Get the gradient of the non-local micro-deformation w.r.t. the non-local reference relative position vector
         */

        if ( !_dNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( ) );

        }

        return &_dNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector.second;

    }

    void aspBase::setdNonLocalMicroDeformationdNonLocalMicroDeformationBase( ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the non-local micro-deformation base given
         * 
         * \f$ \chi_{iI}^{\text{NL}} = \chi_{iI}^{\text{NL base}} + \chi_{iI,J} dX_J\f$
         * 
         * where \f$ \chi_{iI}^{\text{NL}} \f$ is the non-local micro-deformation, \f$\chi_{iI}^{\text{NL base}}\f$
         * is the non-local micro-deformation base, \f$\chi_{iI,J}\f$ is the gradient of the micro-deformation w.r.t.
         * the reference spacial coordinate, and \f$dX_J\f$ is the center-to-center spacing of the two particles defined as
         * 
         * \f$dX_I = \Xi_I + D_I - \Xi_I^{\text{NL}} \f$
         * 
         * where \f$\Xi_I\f$ is the local relative position vector in the reference configuration,
         * \f$D_I\f$ is the reference inter-particle spacing, and \f$\Xi_I^{\text{NL}}\f$ is the
         * non-local relative position vector in the reference configuration
         */

        setdNonLocalMicroDeformationdNonLocalMicroDeformationBase( vectorTools::eye< floatType >( getNonLocalMicroDeformation( )->size( ) ) );

    }

    void aspBase::setdNonLocalMicroDeformationdNonLocalMicroDeformationBase( const floatMatrix &value ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the non-local micro-deformation base
         * 
         * \param &value: The value of the gradient of the non-local micro-deformation w.r.t. the non-local micro-deformation base
         *     position vector
         */

        _dNonLocalMicroDeformationdNonLocalMicroDeformationBase.second = value;

        _dNonLocalMicroDeformationdNonLocalMicroDeformationBase.first = true;

        addInteractionPairData( &_dNonLocalMicroDeformationdNonLocalMicroDeformationBase );

    }

    const floatMatrix* aspBase::getdNonLocalMicroDeformationdNonLocalMicroDeformationBase( ){
        /*!
         * Get the gradient of the non-local micro-deformation w.r.t. the non-local micro deformation base
         */

        if ( !_dNonLocalMicroDeformationdNonLocalMicroDeformationBase.first ){

            ERROR_TOOLS_CATCH( setdNonLocalMicroDeformationdNonLocalMicroDeformationBase( ) );

        }

        return &_dNonLocalMicroDeformationdNonLocalMicroDeformationBase.second;

    }

    void aspBase::setdNonLocalMicroDeformationdGradientMicroDeformation( ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the reference coordinate gradient of the micro-deformation
         * 
         * \f$ \chi_{iI}^{\text{NL}} = \chi_{iI}^{\text{NL base}} + \chi_{iI,J} dX_J\f$
         * 
         * where \f$ \chi_{iI}^{\text{NL}} \f$ is the non-local micro-deformation, \f$\chi_{iI}^{\text{NL base}}\f$
         * is the non-local micro-deformation base, \f$\chi_{iI,J}\f$ is the gradient of the micro-deformation w.r.t.
         * the reference spacial coordinate, and \f$dX_J\f$ is the center-to-center spacing of the two particles defined as
         * 
         * \f$dX_I = \Xi_I + D_I - \Xi_I^{\text{NL}} \f$
         * 
         * where \f$\Xi_I\f$ is the local relative position vector in the reference configuration,
         * \f$D_I\f$ is the reference inter-particle spacing, and \f$\Xi_I^{\text{NL}}\f$ is the
         * non-local relative position vector in the reference configuration
         */

        const unsigned int *dim = getDimension( );

        const floatVector *dX;
        ERROR_TOOLS_CATCH( dX = getLocalReferenceParticleSpacingVector( ) );

        floatVector eye( ( *dim ) * ( *dim ), 0 );
        vectorTools::eye( eye );

        floatMatrix value( ( *dim ) * ( *dim ), floatVector( ( *dim ) * ( *dim ) * ( *dim ), 0 ) );

        for ( unsigned int i = 0; i < *dim; i++ ){

            for ( unsigned int I = 0; I < *dim; I++ ){

                for ( unsigned int a = 0; a < *dim; a++ ){

                    for ( unsigned int A = 0; A < *dim; A++ ){

                        for ( unsigned int B = 0; B < *dim; B++ ){

                            value[ ( *dim ) * i + I ][ ( *dim ) * ( *dim ) * a + ( *dim ) * A + B ]
                                += eye[ ( *dim ) * i + a ] * eye[ ( *dim ) * I + A ] * ( *dX )[ B ];

                        }

                    }

                }

            }

        }

        setdNonLocalMicroDeformationdGradientMicroDeformation( value );

    }

    void aspBase::setdNonLocalMicroDeformationdGradientMicroDeformation( const floatMatrix &value ){
        /*!
         * Set the derivative of the non-local micro deformation w.r.t. the reference coordinate gradient of the micro-deformation
         * 
         * \param &value: The value of the gradient of the non-local micro-deformation w.r.t. the reference coordinate gradient
         *     of the micro-deformation.
         */

        _dNonLocalMicroDeformationdGradientMicroDeformation.second = value;

        _dNonLocalMicroDeformationdGradientMicroDeformation.first = true;

        addInteractionPairData( &_dNonLocalMicroDeformationdGradientMicroDeformation );

    }

    const floatMatrix* aspBase::getdNonLocalMicroDeformationdGradientMicroDeformation( ){
        /*!
         * Get the gradient of the non-local micro-deformation w.r.t. the reference coordinate gradient of the micro-deformation
         */

        if ( !_dNonLocalMicroDeformationdGradientMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setdNonLocalMicroDeformationdGradientMicroDeformation( ) );

        }

        return &_dNonLocalMicroDeformationdGradientMicroDeformation.second;

    }

    void aspBase::setNonLocalMicroDeformation( const floatVector &value ){
        /*!
         * Set the non-local micro deformation
         * 
         * \param &value: The value of the non-local micro deformation
         */

        _nonLocalMicroDeformation.second = value;

        _nonLocalMicroDeformation.first = true;

        addInteractionPairData( &_nonLocalMicroDeformation );

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

        setNonLocalMicroDeformationBase( *microDeformationTensor );

    }

    void aspBase::setNonLocalMicroDeformationBase( const floatVector &value ){
        /*!
         * Set the non-local micro-deformation tensor base
         * 
         * \param &value: The value of the non-local micro-deformation tensor base
         */

        _nonLocalMicroDeformationBase.second = value;

        _nonLocalMicroDeformationBase.first = true;

        addInteractionPairData( &_nonLocalMicroDeformationBase );

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

        floatVector currentDistanceVector;

        floatMatrix dddXi, dddXiNL, dddD, dddF, dddchi, dddchiNL;

        ERROR_TOOLS_CATCH( tractionSeparation::computeCurrentDistanceGeneral( *localSurfaceReferenceRelativePositionVector,
                                                                              *nonLocalSurfaceReferenceRelativePositionVector,
                                                                              *referenceDistanceVector,
                                                                              *localDeformationGradient,
                                                                              *localMicroDeformation,
                                                                              *nonLocalMicroDeformation,
                                                                              currentDistanceVector,
                                                                              dddXi, dddXiNL, dddD, dddF, dddchi, dddchiNL ) );

        setCurrentDistanceVector( currentDistanceVector );

        // Set first order derivatives
        dddXi   += vectorTools::dot( dddchiNL, *getdNonLocalMicroDeformationdLocalReferenceRelativePositionVector( ) );

        dddXiNL += vectorTools::dot( dddchiNL, *getdNonLocalMicroDeformationdNonLocalReferenceRelativePositionVector( ) );

        dddD    += vectorTools::dot( dddchiNL, *getdNonLocalMicroDeformationdLocalReferenceDistanceVector( ) );

        floatMatrix dddChiNLBase = vectorTools::dot( dddchiNL, *getdNonLocalMicroDeformationdNonLocalMicroDeformationBase( ) );

        floatMatrix dddGradChi = vectorTools::dot( dddchiNL, *getdNonLocalMicroDeformationdGradientMicroDeformation( ) );

        setdCurrentDistanceVectordLocalReferenceRelativePositionVector( dddXi );

        setdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( dddXiNL );

        setdCurrentDistanceVectordReferenceDistanceVector( dddD );

        setdCurrentDistanceVectordLocalMicroDeformation( dddchi );

        setdCurrentDistanceVectordNonLocalMicroDeformationBase( dddChiNLBase );

        setdCurrentDistanceVectordGradientMicroDeformation( dddGradChi );

        return;

    }

    void aspBase::setCurrentDistanceVector( const floatVector &value ){
        /*!
         * Set the current distance vector
         * 
         * \param &value: The value of the current distance vector
         */

        _currentDistanceVector.second = value;

        _currentDistanceVector.first = true;

        addInteractionPairData( &_currentDistanceVector );

        return;

    }

    void aspBase::setdCurrentDistanceVectordLocalReferenceRelativePositionVector( ){
        /*!
         * set the derivative of the current distance vector w.r.t. the local reference
         * relative position vector.
         */

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

    }

    void aspBase::setdCurrentDistanceVectordLocalReferenceRelativePositionVector( const floatMatrix &value ){
        /*!
         * set the derivative of the current distance vector w.r.t. the local reference
         * relative position vector.
         * 
         * \param &value: The derivative of the current distance vector w.r.t. the local
         *     reference relative position vector
         */

        _dCurrentDistanceVectordLocalReferenceRelativePositionVector.second = value;

        _dCurrentDistanceVectordLocalReferenceRelativePositionVector.first = true;

        addInteractionPairData( &_dCurrentDistanceVectordLocalReferenceRelativePositionVector );

    }

    const floatMatrix *aspBase::getdCurrentDistanceVectordLocalReferenceRelativePositionVector( ){
        /*!
         * Get the derivative of the current distance vector w.r.t. the local reference
         * relative position vector.
         */

        if ( !_dCurrentDistanceVectordLocalReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setdCurrentDistanceVectordLocalReferenceRelativePositionVector( ) )

        }

        return &_dCurrentDistanceVectordLocalReferenceRelativePositionVector.second;

    }

    void aspBase::setdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ){
        /*!
         * set the derivative of the current distance vector w.r.t. the non-local reference
         * relative position vector.
         */

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

    }

    void aspBase::setdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( const floatMatrix &value ){
        /*!
         * set the derivative of the current distance vector w.r.t. the non-local reference
         * relative position vector.
         * 
         * \param &value: The value of the derivative of the current distance vector w.r.t. the non-local relative position vector.
         */

        _dCurrentDistanceVectordNonLocalReferenceRelativePositionVector.second = value;

        _dCurrentDistanceVectordNonLocalReferenceRelativePositionVector.first = true;

        addInteractionPairData( &_dCurrentDistanceVectordNonLocalReferenceRelativePositionVector );

    }

    const floatMatrix *aspBase::getdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ){
        /*!
         * Get the derivative of the current distance vector w.r.t. the non-local reference
         * relative position vector.
         */

        if ( !_dCurrentDistanceVectordNonLocalReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setdCurrentDistanceVectordNonLocalReferenceRelativePositionVector( ) )

        }

        return &_dCurrentDistanceVectordNonLocalReferenceRelativePositionVector.second;

    }

    void aspBase::setdCurrentDistanceVectordReferenceDistanceVector( ){
        /*!
         * set the derivative of the current distance vector w.r.t. the reference
         * distance vector
         */

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

    }

    void aspBase::setdCurrentDistanceVectordReferenceDistanceVector( const floatMatrix &value ){
        /*!
         * set the derivative of the current distance vector w.r.t. the reference
         * distance vector
         * 
         * \param &value: The derivative of the current distance vector w.r.t. the reference distance vector
         */

        _dCurrentDistanceVectordReferenceDistanceVector.second = value;

        _dCurrentDistanceVectordReferenceDistanceVector.first = true;

        addInteractionPairData( &_dCurrentDistanceVectordReferenceDistanceVector );

    }

    const floatMatrix *aspBase::getdCurrentDistanceVectordReferenceDistanceVector( ){
        /*!
         * Get the derivative of the current distance vector w.r.t. the reference distance vector.
         */

        if ( !_dCurrentDistanceVectordReferenceDistanceVector.first ){

            ERROR_TOOLS_CATCH( setdCurrentDistanceVectordReferenceDistanceVector( ) )

        }

        return &_dCurrentDistanceVectordReferenceDistanceVector.second;

    }

    void aspBase::setdCurrentDistanceVectordLocalDeformationGradient( ){
        /*!
         * set the derivative of the current distance vector w.r.t. the local
         * deformation gradient.
         */

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

    }

    void aspBase::setdCurrentDistanceVectordLocalDeformationGradient( const floatMatrix &value ){
        /*!
         * set the derivative of the current distance vector w.r.t. the local
         * deformation gradient.
         * 
         * \param &value: The derivative of the current distance vector w.r.t. the local
         *     deformation gradient.
         */

        _dCurrentDistanceVectordLocalDeformationGradient.second = value;

        _dCurrentDistanceVectordLocalDeformationGradient.first = true;

        addInteractionPairData( &_dCurrentDistanceVectordLocalDeformationGradient );

    }

    const floatMatrix *aspBase::getdCurrentDistanceVectordLocalDeformationGradient( ){
        /*!
         * Get the derivative of the current distance vector w.r.t. the local deformation gradient.
         */

        if ( !_dCurrentDistanceVectordLocalDeformationGradient.first ){

            ERROR_TOOLS_CATCH( setdCurrentDistanceVectordLocalDeformationGradient( ) )

        }

        return &_dCurrentDistanceVectordLocalDeformationGradient.second;

    }

    void aspBase::setdCurrentDistanceVectordLocalMicroDeformation( ){
        /*!
         * set the derivative of the current distance vector w.r.t. the local
         * micro-deformation
         */

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

    }

    void aspBase::setdCurrentDistanceVectordLocalMicroDeformation( const floatMatrix &value ){
        /*!
         * set the derivative of the current distance vector w.r.t. the local
         * micro-deformation
         * 
         * \param &value: The derivative of the current distance vector w.r.t. the local
         *     micro-deformation
         */

        _dCurrentDistanceVectordLocalMicroDeformation.second = value;

        _dCurrentDistanceVectordLocalMicroDeformation.first = true;

        addInteractionPairData( &_dCurrentDistanceVectordLocalMicroDeformation );

    }

    const floatMatrix *aspBase::getdCurrentDistanceVectordLocalMicroDeformation( ){
        /*!
         * Get the derivative of the current distance vector w.r.t. the local micro-deformation.
         */

        if ( !_dCurrentDistanceVectordLocalMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setdCurrentDistanceVectordLocalMicroDeformation( ) )

        }

        return &_dCurrentDistanceVectordLocalMicroDeformation.second;

    }

    void aspBase::setdCurrentDistanceVectordNonLocalMicroDeformationBase( ){
        /*!
         * set the derivative of the current distance vector w.r.t. the non-local
         * micro-deformation
         */

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

    }

    void aspBase::setdCurrentDistanceVectordNonLocalMicroDeformationBase( const floatMatrix &value ){
        /*!
         * set the derivative of the current distance vector w.r.t. the non-local
         * micro-deformation
         * 
         * \param &value: The derivative of the current distance vector w.r.t. the non-local micro-deformation
         */

        _dCurrentDistanceVectordNonLocalMicroDeformationBase.second = value;

        _dCurrentDistanceVectordNonLocalMicroDeformationBase.first = true;

        addInteractionPairData( &_dCurrentDistanceVectordNonLocalMicroDeformationBase );

    }

    const floatMatrix *aspBase::getdCurrentDistanceVectordNonLocalMicroDeformationBase( ){
        /*!
         * Get the derivative of the current distance vector w.r.t. the non-local micro-deformation.
         */

        if ( !_dCurrentDistanceVectordNonLocalMicroDeformationBase.first ){

            ERROR_TOOLS_CATCH( setdCurrentDistanceVectordNonLocalMicroDeformationBase( ) )

        }

        return &_dCurrentDistanceVectordNonLocalMicroDeformationBase.second;

    }

    void aspBase::setdCurrentDistanceVectordGradientMicroDeformation( ){
        /*!
         * set the derivative of the current distance vector w.r.t. the reference
         * spatial gradient of the micro-deformation.
         */

        ERROR_TOOLS_CATCH( setCurrentDistanceVector( ) );

    }

    void aspBase::setdCurrentDistanceVectordGradientMicroDeformation( const floatMatrix &value ){
        /*!
         * set the derivative of the current distance vector w.r.t. the reference
         * spatial gradient of the micro-deformation.
         * 
         * \param &value:
         */

        _dCurrentDistanceVectordGradientMicroDeformation.second = value;

        _dCurrentDistanceVectordGradientMicroDeformation.first = true;

        addInteractionPairData( &_dCurrentDistanceVectordGradientMicroDeformation );

    }

    const floatMatrix *aspBase::getdCurrentDistanceVectordGradientMicroDeformation( ){
        /*!
         * Get the derivative of the current distance vector w.r.t. the gradient of the micro-deformation.
         */

        if ( !_dCurrentDistanceVectordGradientMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setdCurrentDistanceVectordGradientMicroDeformation( ) )

        }

        return &_dCurrentDistanceVectordGradientMicroDeformation.second;

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

        floatVector dan, localCurrentNormal;

        floatMatrix ddandChi, ddandN;

        ERROR_TOOLS_CATCH( tractionSeparation::computeNansonsRelation( *localMicroDeformation, *localReferenceNormal,
                                                                       dan, ddandChi, ddandN ) );

        floatType da = vectorTools::l2norm( dan );

        localCurrentNormal = dan / da;

        floatMatrix dndChi = ddandChi / da - vectorTools::dyadic( dan, vectorTools::Tdot( ddandChi, localCurrentNormal ) ) / ( da * da );

        floatMatrix dndN = ddandN / da - vectorTools::dyadic( dan, vectorTools::Tdot( ddandN, localCurrentNormal ) ) / ( da * da );

        setLocalCurrentNormal( localCurrentNormal );

        setdLocalCurrentNormaldLocalMicroDeformation( dndChi );

        setdLocalCurrentNormaldLocalReferenceNormal( dndN );

        return;

    }

    void aspBase::setLocalCurrentNormal( const floatVector &value ){
        /*!
         * Set the current local normal
         * 
         * \param &value: The value of the local current normal
         */

        _localCurrentNormal.second = value;

        _localCurrentNormal.first = true;

        addSurfacePointData( &_localCurrentNormal );

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

    void aspBase::setdLocalCurrentNormaldLocalReferenceNormal( ){
        /*!
         * Set the derivative of the local current normal vector w.r.t. the local reference normal
         */

        ERROR_TOOLS_CATCH( setLocalCurrentNormal( ) );

    }

    void aspBase::setdLocalCurrentNormaldLocalReferenceNormal( const floatMatrix &value ){
        /*!
         * Set the derivative of the local current normal vector w.r.t. the local reference normal
         * 
         * \param &value: The gradient of the local current normal w.r.t. the reference current normal
         */

        _dLocalCurrentNormaldLocalReferenceNormal.second = value;

        _dLocalCurrentNormaldLocalReferenceNormal.first = true;

        addSurfacePointData( &_dLocalCurrentNormaldLocalReferenceNormal );

    }

    const floatMatrix* aspBase::getdLocalCurrentNormaldLocalReferenceNormal( ){
        /*!
         * Get the derivative of the local current normal w.r.t. the local reference normal
         */

        if ( !_dLocalCurrentNormaldLocalReferenceNormal.first ){

            ERROR_TOOLS_CATCH( setdLocalCurrentNormaldLocalReferenceNormal( ) );

        }

        return &_dLocalCurrentNormaldLocalReferenceNormal.second;

    }

    void aspBase::setdLocalCurrentNormaldLocalMicroDeformation( ){
        /*!
         * Set the derivative of the local current normal vector w.r.t. the local micro-deformation
         */

        ERROR_TOOLS_CATCH( setLocalCurrentNormal( ) );

    }

    void aspBase::setdLocalCurrentNormaldLocalMicroDeformation( const floatMatrix &value ){
        /*!
         * Set the derivative of the local current normal vector w.r.t. the local micro-deformation
         * 
         * \param &value: The gradient of the local current normal w.r.t. the local micro-deformation
         */

        _dLocalCurrentNormaldLocalMicroDeformation.second = value;

        _dLocalCurrentNormaldLocalMicroDeformation.first = true;

        addSurfacePointData( &_dLocalCurrentNormaldLocalMicroDeformation );

    }

    const floatMatrix* aspBase::getdLocalCurrentNormaldLocalMicroDeformation( ){
        /*!
         * Get the derivative of the local current normal w.r.t. the local micro-deformation
         */

        if ( !_dLocalCurrentNormaldLocalMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setdLocalCurrentNormaldLocalMicroDeformation( ) );

        }

        return &_dLocalCurrentNormaldLocalMicroDeformation.second;

    }

    void aspBase::getLocalReferenceNormal( const unsigned int &index, floatVector &localReferenceNormal ){
        /*!
         * compute the local reference normal at the referenced surface point index
         * 
         * \param &index: The index of the local surface point
         * \param &localReferenceNormal: The resulting local reference normal vector
         */

        const unsigned int *dim = getDimension( );

        const floatVector* unitSpherePoints;
        ERROR_TOOLS_CATCH( unitSpherePoints = getUnitSpherePoints( ) );

        if ( ( *dim ) * ( index + 1 ) > unitSpherePoints->size( ) ){

            std::string message = "The requested index is greater than the number of points available on the unit sphere.\n";
            message            += "  index:            " + std::to_string( index ) + "\n";
            message            += "  number of points: " + std::to_string( unitSpherePoints->size( ) / ( *dim ) );

            ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        localReferenceNormal = floatVector( unitSpherePoints->begin( ) + ( *dim ) * index,
                                            unitSpherePoints->begin( ) + ( *dim ) * ( index + 1 ) );

        localReferenceNormal /= vectorTools::l2norm( localReferenceNormal );

        return;

    }

    void aspBase::getLocalCurrentNormal( const unsigned int &index, floatVector &normal ){
        /*!
         * compute the local current normal at the referenced surface point index
         * 
         * \param &index: The index of the local surface point
         * \param &normal: The resulting normal vector
         */

        floatVector referenceNormal;
        ERROR_TOOLS_CATCH( getLocalReferenceNormal( index, referenceNormal ) );

        const floatVector* localMicroDeformation;
        ERROR_TOOLS_CATCH( localMicroDeformation = getLocalMicroDeformation( ) );

        ERROR_TOOLS_CATCH( tractionSeparation::computeNansonsRelation( *localMicroDeformation, referenceNormal,
                                                                       normal ) );

        normal /= vectorTools::l2norm( normal );

        return;

    }

    void aspBase::setSurfaceParameters( ){
        /*!
         * Set the surface parameters for the integrand
         */

        setSurfaceParameters( floatVector( 0, 0 ) );        

        return;
    }

    void aspBase::setSurfaceParameters( const floatVector &value ){
        /*!
         * Set the surface parameters for the integrand
         * 
         * \param &value: The values of the surface parameters
         */

        _surfaceParameters.second = value;

        _surfaceParameters.first = true;

        addInteractionPairData( &_surfaceParameters );

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


        setReferenceDistanceVector( floatVector( *getDimension( ), 0. ) );

        return;

    }

    void aspBase::setReferenceDistanceVector( const floatVector &value ){
        /*!
         * Set the initial distance between two particles
         * 
         * \param &value: The value of the reference distance vector
         */

        _referenceDistanceVector.second = value;

        _referenceDistanceVector.first = true;

        addInteractionPairData( &_referenceDistanceVector );

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

    void aspBase::resetInteractionPairData( ){
        /*!
         * Reset the interaction pair data to the base state
         */

        for ( auto d = _interactionPairData.begin( ); d != _interactionPairData.end( ); d++ ){

            ( *d )->clear( );

        }

        _interactionPairData.clear( );

        return;

    }

    void aspBase::resetSurfacePointData( ){
        /*!
         * Reset the surface point data to the base state
         */

        // Reset the interaction pair
        resetInteractionPairData( );

        for ( auto d = _surfacePointData.begin( ); d != _surfacePointData.end( ); d++ ){

            ( *d )->clear( );

        }

        _surfacePointData.clear( );

        return;

    }

    void aspBase::resetLocalParticleData( ){
        /*!
         * Reset the local particle data to the base state
         */

        resetSurfacePointData( );

        for ( auto d = _localParticleData.begin( ); d != _localParticleData.end( ); d++ ){

            ( *d )->clear( );

        }

        _localParticleData.clear( );

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

        const unsigned int *dim = getDimension( );

        containedPoints = std::vector< unsigned int >( points.size( ) );

        unsigned int numContainedPoints = 0;

        for ( unsigned int i = 0; i < points.size( ); i += ( *dim ) ){

            if ( pointInBoundingBox( floatVector( points.begin( ) + i, points.begin( ) + i + ( *dim ) ), boundingBox ) ){

                containedPoints[ numContainedPoints ] = i / ( *dim );

                numContainedPoints++;

            }

        }

        containedPoints.resize( numContainedPoints );

        return;

    }

    void aspBase::computeSurfaceOverlapEnergyDensity( mapFloatType &surfaceOverlapEnergyDensity ){
        /*!
         * Compute the surface overlap energy density for the local particle and a given interaction
         * 
         * \param &surfaceOverlapEnergyDensities: The index is the index of the local points which may be overlapping with a neighboring particle and the values are the overlap energies. It is not a scalar because two particles could be overlapping at multiple points
         */

        surfaceOverlapEnergyDensity.clear( );

        const mapFloatVector *particlePairOverlap;
        ERROR_TOOLS_CATCH( particlePairOverlap = getParticlePairOverlap( ) );

        const floatVector *overlapParameters;
        ERROR_TOOLS_CATCH( overlapParameters = getSurfaceOverlapParameters( ) );

        for ( auto overlap = particlePairOverlap->begin( ); overlap != particlePairOverlap->end( ); overlap++ ){

            floatVector normal;

            getLocalCurrentNormal( overlap->first, normal );

            surfaceOverlapEnergyDensity.insert( { overlap->first, 0.5 * ( *overlapParameters )[ 0 ] * vectorTools::dot( overlap->second, overlap->second ) * std::fabs( vectorTools::dot( overlap->second, normal ) ) } );

        }

        return;

    }

    void aspBase::setParticlePairOverlap( ){
        /*!
         * Set the particle overlap for the current local to non-local pair
         */

        mapFloatVector particlePairOverlap;

        const unsigned int *dim = getDimension( );

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
        ERROR_TOOLS_CATCH( localReferenceParticleSpacing = getLocalReferenceParticleSpacingVector( ) );

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

            ERROR_TOOLS_CATCH( tractionSeparation::computeParticleOverlap( floatVector( localReferenceSurfacePoints->begin( ) + ( *dim ) * ( *p ),
                                                                                        localReferenceSurfacePoints->begin( ) + ( *dim ) * ( ( *p ) + 1 ) ),
                                                                           *localReferenceParticleSpacing, *nonLocalReferenceRadius, *localDeformationGradient,
                                                                           *localMicroDeformation, *nonLocalMicroDeformationBase, *localGradientMicroDeformation,
                                                                           overlap ) );

            particlePairOverlap.insert( { *p, overlap } );

        }

        setParticlePairOverlap( particlePairOverlap );

        return;

    }

    void aspBase::setParticlePairOverlap( const mapFloatVector &value ){
        /*!
         * Set the vector of the overlap of the current pariticle pair
         * 
         * \param &value: The vector defining the overlap of the two particles.
         */

        _particlePairOverlap.second = value;

        _particlePairOverlap.first = true;

        addInteractionPairData( &_particlePairOverlap );

    }

    const mapFloatVector* aspBase::getParticlePairOverlap( ){
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

    const mapFloatType* aspBase::getSurfaceOverlapEnergyDensity( ){
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

        setLocalGradientMicroDeformation( *gradientMicroDeformation );

    }

    void aspBase::setLocalGradientMicroDeformation( const floatVector &value ){
        /*!
         * Set the local gradient of the micro-deformation
         * 
         * \param &value: The current value of the local gradient of the micro deformation
         */

        _localGradientMicroDeformation.second = value;

        _localGradientMicroDeformation.first = true;

        addLocalParticleData( &_localGradientMicroDeformation );

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

        floatType surfaceAdhesionEnergyDensity;

        ERROR_TOOLS_CATCH( computeSurfaceAdhesionEnergyDensity( surfaceAdhesionEnergyDensity ) );

        setSurfaceAdhesionEnergyDensity( surfaceAdhesionEnergyDensity );

        return;

    }

    void aspBase::setSurfaceAdhesionEnergyDensity( const floatType& value ){
        /*!
         * Set the surface adhesion energy density if required.
         * 
         * \param &value: the current value of the surface adhesion density
         */

        _surfaceAdhesionEnergyDensity.second = value;

        _surfaceAdhesionEnergyDensity.first = true;

        addInteractionPairData( &_surfaceAdhesionEnergyDensity );

        return;

    }

    void aspBase::setSurfaceOverlapEnergyDensity( ){
        /*!
         * Set the surface overlap energy density if required.
         */

        mapFloatType surfaceOverlapEnergyDensity;

        ERROR_TOOLS_CATCH( computeSurfaceOverlapEnergyDensity( surfaceOverlapEnergyDensity ) );

        setSurfaceOverlapEnergyDensity( surfaceOverlapEnergyDensity );

        return;

    }

    void aspBase::setSurfaceOverlapEnergyDensity( const mapFloatType &value ){
        /*!
         * Set the surface overlap energy density if required.
         * 
         * \param &value: The map to the surface overlap energy density
         */

        _surfaceOverlapEnergyDensity.second = value;

        _surfaceOverlapEnergyDensity.first = true;

        addInteractionPairData( &_surfaceOverlapEnergyDensity );

        return;

    }

    void aspBase::setSurfaceOverlapParameters( ){
        /*!
         * Set the surface overlap parameters
         */

        setSurfaceOverlapParameters( { 1 } );

        return;

    }

    void aspBase::setSurfaceOverlapParameters( const floatVector &value ){
        /*!
         * Set the surface overlap parameters
         * 
         * \param &value: Set the value of the surface overlap parameters
         */

        _surfaceOverlapParameters.second = value;

        _surfaceOverlapParameters.first = true;

        addInteractionPairData( &_surfaceOverlapParameters );

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

        setLocalReferenceSurfacePoints( ( *localReferenceRadius ) * ( *unitSpherePoints ) );

        return;

    }

    void aspBase::setLocalReferenceSurfacePoints( const floatVector &value ){
        /*!
         * Set the collection of points which are on the surface of the local particle in the reference configuration
         * \param &value: The values of the local reference surface points
         */

        _localReferenceSurfacePoints.second = value;

        _localReferenceSurfacePoints.first = true;

        addLocalParticleData( &_localReferenceSurfacePoints );

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

        setNonLocalReferenceSurfacePoints( ( *nonLocalReferenceRadius ) * ( *unitSpherePoints ) );

        return;

    }

    void aspBase::setNonLocalReferenceSurfacePoints( const floatVector &value ){
        /*!
         * Set the collection of points which are on the surface of the non-local particle in the reference configuration
         * 
         * \param &value: The value of the non-local reference surface points
         */

        _nonLocalReferenceSurfacePoints.second = value;

        _nonLocalReferenceSurfacePoints.first = true;

        addInteractionPairData( &_nonLocalReferenceSurfacePoints );

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

        const unsigned int *dim = getDimension( );

        const floatVector *localReferenceSurfacePoints;
        ERROR_TOOLS_CATCH( localReferenceSurfacePoints = getLocalReferenceSurfacePoints( ) );

        const floatVector *localMicroDeformation;
        ERROR_TOOLS_CATCH( localMicroDeformation = getLocalMicroDeformation( ) );

        floatVector localCurrentSurfacePoints;

        ERROR_TOOLS_CATCH( localCurrentSurfacePoints = vectorTools::matrixMultiply( *localReferenceSurfacePoints, *localMicroDeformation,
                                                                                    localReferenceSurfacePoints->size( ) / ( *dim ), *dim, *dim, *dim, false, true ) );

        setLocalCurrentSurfacePoints( localCurrentSurfacePoints );

    }

    void aspBase::setLocalCurrentSurfacePoints( const floatVector &value ){
        /*!
         * Set the collection of points on the surface of the local particle in the current configuration
         *
         * \param &value: The current points on the surface of the local particle
         */

        _localCurrentSurfacePoints.second = value;

        _localCurrentSurfacePoints.first = true;

        addLocalParticleData( &_localCurrentSurfacePoints );

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

        const unsigned int *dim = getDimension( );

        const floatVector *nonLocalReferenceSurfacePoints;
        ERROR_TOOLS_CATCH( nonLocalReferenceSurfacePoints = getNonLocalReferenceSurfacePoints( ) );

        const floatVector *nonLocalMicroDeformation;
        ERROR_TOOLS_CATCH( nonLocalMicroDeformation = getNonLocalMicroDeformation( ) );

        floatVector nonLocalCurrentSurfacePoints;

        ERROR_TOOLS_CATCH( nonLocalCurrentSurfacePoints = vectorTools::matrixMultiply( *nonLocalReferenceSurfacePoints, *nonLocalMicroDeformation,
                                                                                       nonLocalReferenceSurfacePoints->size( ) / ( *dim ),
                                                                                               *dim, *dim, *dim, false, true ) );

        setNonLocalCurrentSurfacePoints( nonLocalCurrentSurfacePoints );

    }

    void aspBase::setNonLocalCurrentSurfacePoints( const floatVector &value ){
        /*!
         * Set the collection of points on the surface of the non-local particle in the current configuration
         * 
         * \param &value: The current non-local surface points
         */

        _nonLocalCurrentSurfacePoints.second = value;

        _nonLocalCurrentSurfacePoints.first = true;

        addInteractionPairData( &_nonLocalCurrentSurfacePoints );

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

        floatMatrix localParticleCurrentBoundingBox;

        ERROR_TOOLS_CATCH( formBoundingBox( *localCurrentSurfacePoints, localParticleCurrentBoundingBox ) );

        setLocalParticleCurrentBoundingBox( localParticleCurrentBoundingBox );

        return;

    }

    void aspBase::setLocalParticleCurrentBoundingBox( const floatMatrix &value ){
        /*!
         * Set the bounding box of the local particle (dimension, (lower bound, upper bound))
         * 
         * \param &value: The local particle's current bounding box
         */

        _localParticleCurrentBoundingBox.second = value;

        _localParticleCurrentBoundingBox.first = true;

        addLocalParticleData( &_localParticleCurrentBoundingBox );

        return;

    }

    void aspBase::formBoundingBox( const floatVector &points, floatMatrix &boundingBox ){
        /*!
         * Form a bounding box from the provided points
         * 
         * \param &points: The points to form the bounding box in [x1, y1, z1, x2, y2, z2, ...] format
         * \param &boundingBox: The resulting bounding box
         */

        const unsigned int *dim = getDimension( );

        if ( ( points.size( ) % ( *dim ) ) > 0 ){

            ERROR_TOOLS_CATCH( throw std::runtime_error( "The length of points must be an integer multiple of the dimension.\n  points.size( ) = " + std::to_string( points.size( ) ) + "\n  dimension: " + std::to_string( ( *dim ) ) ) );

        }

        if ( points.size( ) < ( *dim ) ){

            ERROR_TOOLS_CATCH( throw std::runtime_error( "The points need at least one point" ) );

        }

        boundingBox = { { points[ 0 ], points[ 0 ] },
                        { points[ 1 ], points[ 1 ] },
                        { points[ 2 ], points[ 2 ] } };

        for ( unsigned int i = ( *dim ); i < points.size( ); i += ( *dim ) ){

            for ( unsigned int j = 0; j < ( *dim ); j++ ){

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

        floatMatrix nonLocalParticleCurrentBoundingBox;

        ERROR_TOOLS_CATCH( formBoundingBox( *nonLocalCurrentSurfacePoints, nonLocalParticleCurrentBoundingBox ) );

        setNonLocalParticleCurrentBoundingBox( nonLocalParticleCurrentBoundingBox );

        return;

    }

    void aspBase::setNonLocalParticleCurrentBoundingBox( const floatMatrix &value ){
        /*!
         * Set the bounding box of the non-local particle (dimension, (lower bound, upper bound))
         * 
         * \param &value: The non-local particle's current bounding box
         */

        _nonLocalParticleCurrentBoundingBox.second = value;

        _nonLocalParticleCurrentBoundingBox.first = true;

        addInteractionPairData( &_nonLocalParticleCurrentBoundingBox );

        return;

    }

    void aspBase::setSurfaceAdhesionTraction( ){
        /*!
         * Set the surface adhesion traction
         */

        floatVector surfaceAdhesionTraction;

        ERROR_TOOLS_CATCH( computeSurfaceAdhesionTraction( surfaceAdhesionTraction ) );

        setSurfaceAdhesionTraction( surfaceAdhesionTraction );

        return;

    }

    void aspBase::setSurfaceAdhesionTraction( const floatVector &value ){
        /*!
         * Set the surface adhesion traction
         *
         * \param &value: The value of the surface adhesion traction.
         */

        _surfaceAdhesionTraction.second = value;

        _surfaceAdhesionTraction.first = true;

        addInteractionPairData( &_surfaceAdhesionTraction );

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

    void aspBase::setSurfaceAdhesionThickness( ){
        /*!
         * Set the current thickness of the surface adhesion
         */

        const floatVector* currentDistanceVector;
        ERROR_TOOLS_CATCH( currentDistanceVector = getCurrentDistanceVector( ) );

        const floatVector* localCurrentNormal;
        ERROR_TOOLS_CATCH( localCurrentNormal = getLocalCurrentNormal( ) );
        
        // Decompose the distance into normal and tangential directions
        floatVector dn, dt;
        ERROR_TOOLS_CATCH_NODE_POINTER( tractionSeparation::decomposeVector( *currentDistanceVector, *localCurrentNormal, dn, dt ) );

        setSurfaceAdhesionThickness( vectorTools::l2norm( dn ) );

    }

    void aspBase::setSurfaceAdhesionThickness( const floatType &value ){
        /*!
         * Set the current thickness of the surface adhesion
         * 
         * \param &value: The thickness associated with the surface adhesion
         */

        _surfaceAdhesionThickness.second = value;

        _surfaceAdhesionThickness.first = true;

        addInteractionPairData( &_surfaceAdhesionThickness );

    }

    const floatType* aspBase::getSurfaceAdhesionThickness( ){
        /*!
         * Get the current value of the surface adhesion thickness
         */

        if ( !_surfaceAdhesionThickness.first ){

            ERROR_TOOLS_CATCH( setSurfaceAdhesionThickness( ) );

        }

        return &_surfaceAdhesionThickness.second;

    }

    void aspBase::setSurfaceOverlapThickness( ){
        /*!
         * Set the current thickness of the surface overlap
         */

        const mapFloatVector* particlePairOverlap;
        ERROR_TOOLS_CATCH( particlePairOverlap = getParticlePairOverlap( ) );

        mapFloatType surfaceOverlapThickness;

        for ( auto overlap = particlePairOverlap->begin( ); overlap != particlePairOverlap->end( ); overlap++ ){

            floatVector normal;

            ERROR_TOOLS_CATCH( getLocalCurrentNormal( overlap->first, normal ) );
        
            surfaceOverlapThickness.insert( { overlap->first, std::fabs( vectorTools::dot( overlap->second, normal ) ) } );

        }

        setSurfaceOverlapThickness( surfaceOverlapThickness );

    }

    void aspBase::setSurfaceOverlapThickness( const mapFloatType &value ){
        /*!
         * Set the current thickness of the surface overlap
         * 
         * \param &value: The map from local surface node indices to the overlap thickness
         */

        _surfaceOverlapThickness.second = value;

        _surfaceOverlapThickness.first = true;

        addInteractionPairData( &_surfaceOverlapThickness );

    }

    const mapFloatType* aspBase::getSurfaceOverlapThickness( ){
        /*!
         * Get the current value of the surface overlap thickness
         */

        if ( !_surfaceOverlapThickness.first ){

            ERROR_TOOLS_CATCH( setSurfaceOverlapThickness( ) );

        }

        return &_surfaceOverlapThickness.second;

    }

    const std::vector< std::vector< floatVector > >* aspBase::getAssembledSurfaceAdhesionThicknesses( ){
        /*!
         * Get the current value of the surface adhesion thicknesses
         */

        if ( !_assembledSurfaceAdhesionThicknesses.first ){

            ERROR_TOOLS_CATCH( assembleSurfaceResponses( ) );

        }

        return &_assembledSurfaceAdhesionThicknesses.second;

    }

    const std::vector< std::vector< floatVector > >* aspBase::getAssembledSurfaceAdhesionEnergyDensities( ){
        /*!
         * Get the current value of the surface adhesion energy densities
         */

        if ( !_assembledSurfaceAdhesionEnergyDensities.first ){

            ERROR_TOOLS_CATCH( assembleSurfaceResponses( ) );

        }

        return &_assembledSurfaceAdhesionEnergyDensities.second;

    }

    const std::vector< std::vector< floatMatrix > >* aspBase::getAssembledSurfaceAdhesionTractions( ){
        /*!
         * Get the current value of the surface adhesion thicknesses
         */

        if ( !_assembledSurfaceAdhesionTractions.first ){

            ERROR_TOOLS_CATCH( assembleSurfaceResponses( ) );

        }

        return &_assembledSurfaceAdhesionTractions.second;

    }

    const std::vector< std::vector< std::vector< mapFloatType > > >* aspBase::getAssembledSurfaceOverlapThicknesses( ){
        /*!
         * Get the current value of the surface overlap thicknesses
         */

        if ( !_assembledSurfaceOverlapThicknesses.first ){

            ERROR_TOOLS_CATCH( assembleSurfaceResponses( ) );

        }

        return &_assembledSurfaceOverlapThicknesses.second;

    }

    const std::vector< std::vector< std::vector< mapFloatType > > >* aspBase::getAssembledSurfaceOverlapEnergyDensities( ){
        /*!
         * Get the current value of the surface overlap energy densities
         */

        if ( !_assembledSurfaceOverlapEnergyDensities.first ){

            ERROR_TOOLS_CATCH( assembleSurfaceResponses( ) );

        }

        return &_assembledSurfaceOverlapEnergyDensities.second;

    }

    const std::vector< std::vector< std::vector< mapFloatVector > > >* aspBase::getAssembledSurfaceOverlapTractions( ){
        /*!
         * Get the current value of the surface overlap thicknesses
         */

        if ( !_assembledSurfaceOverlapTractions.first ){

            ERROR_TOOLS_CATCH( assembleSurfaceResponses( ) );

        }

        return &_assembledSurfaceOverlapTractions.second;

    }

    void aspBase::computeSurfaceOverlapTraction( mapFloatVector &surfaceOverlapTraction ){
        /*!
         * Compute the surface overlap traction
         * 
         * \param &surfaceOverlapTraction: An unordered map mapping from points on the local particle to the amount they overlap the non-local particle
         */

        surfaceOverlapTraction.clear( );

        const mapFloatVector *particlePairOverlap;
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

        mapFloatVector surfaceOverlapTraction;

        ERROR_TOOLS_CATCH( computeSurfaceOverlapTraction( surfaceOverlapTraction ) );

        setSurfaceOverlapTraction( surfaceOverlapTraction );

        return;

    }

    void aspBase::setSurfaceOverlapTraction( const mapFloatVector &value ){
        /*!
         * Set the surface overlap traction
         * 
         * \param &value: The map from the local node indices to the overlap traction vectors
         */

        _surfaceOverlapTraction.second = value;

        _surfaceOverlapTraction.first = true;

        addInteractionPairData( &_surfaceOverlapTraction );

        return;

    }

    const mapFloatVector* aspBase::getSurfaceOverlapTraction( ){
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

    const unsigned int* aspBase::getNumLocalParticles( ){
        /*!
         * Get the number of local particles
         */

        return &_numLocalParticles;
    }

    void aspBase::setLocalParticleQuantities( ){
        /*!
         * Set the energy, microstress, and log probability ratio for the local particle
         */

        const floatType *previousTime;
        ERROR_TOOLS_CATCH( previousTime = getPreviousTime( ) );

        const floatType *deltaTime;
        ERROR_TOOLS_CATCH( deltaTime = getDeltaTime( ) );

        const floatVector *currentLocalMicroDeformation;
        ERROR_TOOLS_CATCH( currentLocalMicroDeformation = getLocalMicroDeformation( ) );

        const floatVector *previousLocalMicroDeformation;
        ERROR_TOOLS_CATCH( previousLocalMicroDeformation = getPreviousLocalMicroDeformation( ) );

        const floatType *currentTemperature;
        ERROR_TOOLS_CATCH( currentTemperature = getTemperature( ) );

        const floatType *previousTemperature;
        ERROR_TOOLS_CATCH( previousTemperature = getPreviousTemperature( ) );

        const floatVector *previousStateVariables;
        ERROR_TOOLS_CATCH( previousStateVariables = getPreviousLocalStateVariables( ) );

        const floatVector *parameters;
        ERROR_TOOLS_CATCH( parameters = getLocalParticleParameters( ) );

        floatType localParticleEnergyDensity;

        floatVector localParticleStateVariables;

        floatVector localParticleMicroCauchyStress;

        floatType localParticleLogProbabilityRatio;

        ERROR_TOOLS_CATCH( computeLocalParticleEnergyDensity( *previousTime, *deltaTime, *currentLocalMicroDeformation, *previousLocalMicroDeformation,
                                                              *currentTemperature, *previousTemperature, *previousStateVariables, *parameters,
                                                              localParticleEnergyDensity, localParticleMicroCauchyStress,
                                                              localParticleStateVariables, localParticleLogProbabilityRatio ) );

        setLocalParticleEnergyDensity( localParticleEnergyDensity );

        setLocalParticleMicroCauchyStress( localParticleMicroCauchyStress );

        setLocalParticleStateVariables( localParticleStateVariables );

        setLocalParticleLogProbabilityRatio( localParticleLogProbabilityRatio );

    }

    void aspBase::setLocalParticleEnergy( ){
        /*!
         * Set the local particle's energy
         */

        floatType localParticleEnergy;

        localParticleEnergy = ( *getLocalParticleEnergyDensity( ) ) * ( *getLocalParticleCurrentVolume( ) );

        setLocalParticleEnergy( localParticleEnergy );

    }

    void aspBase::setLocalParticleEnergy( const floatType &value ){
        /*!
         * Set the local particle's energy
         * 
         * \param &value: The local particle's energy
         */

        _localParticleEnergy.second = value;

        _localParticleEnergy.first = true;

        addLocalParticleData( &_localParticleEnergy );

    }

    const floatType* aspBase::getLocalParticleEnergy( ){
        /*!
         * Get the local particle's energy
         */

        if ( !_localParticleEnergy.first ){

            ERROR_TOOLS_CATCH( setLocalParticleEnergy( ) );

        }

        return &_localParticleEnergy.second;

    }

    void aspBase::setLocalParticleEnergyDensity( const floatType &value ){
        /*!
         * Set the local particle's energy density
         *
         * \param &value: The value of the local particle's energy density
         */

        _localParticleEnergyDensity.second = value;

        _localParticleEnergyDensity.first = true;

        addLocalParticleData( &_localParticleEnergyDensity );

    }

    const floatType* aspBase::getLocalParticleEnergyDensity( ){
        /*!
         * Get the local particle's energy density
         */

        if ( !_localParticleEnergyDensity.first ){

            ERROR_TOOLS_CATCH( setLocalParticleQuantities( ) );

        }

        return &_localParticleEnergyDensity.second;

    }

    void aspBase::setLocalParticleLogProbabilityRatio( const floatType &value ){
        /*!
         * Set the local particle's log-probability ratio (i.e. the change in
         * the local particle's energy log probability (log(p_current) - log(p_previous))
         * 
         * \param &value: The value of the change in log probability.
         */

        _localParticleLogProbabilityRatio.second = value;

        _localParticleLogProbabilityRatio.first = true;

        addLocalParticleData( &_localParticleLogProbabilityRatio );

    }

    const floatType* aspBase::getLocalParticleLogProbabilityRatio( ){
        /*!
         * Get the local particle's log probability ratio
         */

        if ( !_localParticleLogProbabilityRatio.first ){

            ERROR_TOOLS_CATCH( setLocalParticleQuantities( ) );

        }

        return &_localParticleLogProbabilityRatio.second;

    }

    void aspBase::setLocalParticleMicroCauchyStress( const floatVector &value ){
        /*!
         * Set the local particle's value of the micro Cauchy stress
         * 
         * \param &value: The value of the micro-cauchy stress
         */

        _localParticleMicroCauchyStress.second = value;

        _localParticleMicroCauchyStress.first = true;

        addLocalParticleData( &_localParticleMicroCauchyStress );

    }

    const floatVector* aspBase::getLocalParticleMicroCauchyStress( ){
        /*!
         * Get the local particle's micro cauchy stress
         */

        if ( !_localParticleMicroCauchyStress.first ){

            ERROR_TOOLS_CATCH( setLocalParticleQuantities( ) );

        }

        return &_localParticleMicroCauchyStress.second;

    }

    void aspBase::setLocalParticleStateVariables( const floatVector &value ){
        /*!
         * Set the local particle's state variables
         * 
         * \param &value: The value of the local particle's state variables
         */

        _localParticleStateVariables.second = value;

        _localParticleStateVariables.first = true;

        addLocalParticleData( &_localParticleStateVariables );

    }

    const floatVector* aspBase::getLocalParticleStateVariables( ){
        /*!
         * Get the local particle's state variables
         */

        if ( !_localParticleStateVariables.first ){

            ERROR_TOOLS_CATCH( setLocalParticleQuantities( ) );

        }

        return &_localParticleStateVariables.second;

    }

    const floatType* aspBase::getPreviousTime( ){
        /*!
         * Get the previous value of the time
         */

        return &_previousTime;

    }
    
    const floatType* aspBase::getDeltaTime( ){
        /*!
         * Get the change in time
         */

        return &_deltaTime;

    }

    const floatType* aspBase::getPreviousTemperature( ){
        /*!
         * Get the previous value of the temperature
         */

        return &_previousTemperature;

    }
    
    const floatType* aspBase::getTemperature( ){
        /*!
         * Get the current temperature
         */

        return &_temperature;

    }

    const floatType* aspBase::getLocalParticleReferenceVolume( ){
        /*!
         * Get the reference configuration volume for the local particle
         */ 

        if ( !_localParticleReferenceVolume.first ){

            ERROR_TOOLS_CATCH( setLocalParticleReferenceVolume( ) );

        }

        return &_localParticleReferenceVolume.second;

    }

    void aspBase::setLocalParticleReferenceVolume( ){
        /*!
         * Set the reference configuration volume for the local particle
         */ 

        const floatType* localReferenceRadius;
        ERROR_TOOLS_CATCH( localReferenceRadius = getLocalReferenceRadius( ) );

        floatType localParticleReferenceVolume;

        localParticleReferenceVolume = 4./3 * _pi * std::pow( *localReferenceRadius, 3 );

        setLocalParticleReferenceVolume( localParticleReferenceVolume );

        return;

    }

    void aspBase::setLocalParticleReferenceVolume( const floatType &value ){
        /*!
         * Set the reference configuration volume for the local particle
         * 
         * \param &value: The value of the reference configuration volume of the
         *     local particle.
         */

        _localParticleReferenceVolume.second = value;

        _localParticleReferenceVolume.first = true;

        addLocalParticleData( &_localParticleReferenceVolume );

        return;

    }

    const floatType* aspBase::getLocalParticleCurrentVolume( ){
        /*!
         * Set the current configuration volume for the local particle
         */ 

        if ( !_localParticleCurrentVolume.first ){

            ERROR_TOOLS_CATCH( setLocalParticleCurrentVolume( ) );

        }

        return &_localParticleCurrentVolume.second;

    }

    void aspBase::setLocalParticleCurrentVolume( ){
        /*!
         * Set the current configuration local particle volume
         */ 

        const unsigned int *dim = getDimension( );

        const floatVector* microDeformation;
        ERROR_TOOLS_CATCH( microDeformation = getLocalMicroDeformation( ) );

        const floatType* referenceVolume;
        ERROR_TOOLS_CATCH( referenceVolume = getLocalParticleReferenceVolume( ) );

        floatType J;
        ERROR_TOOLS_CATCH( J = vectorTools::determinant( *microDeformation, ( *dim ), ( *dim ) ) );

        floatType currentVolume = J * ( *referenceVolume );

        setLocalParticleCurrentVolume( currentVolume );

    }

    void aspBase::setLocalParticleCurrentVolume( const floatType &value ){
        /*!
         * Set the current configuration local particle volume
         * 
         * \param &value: The current configuration local particle volume
         */

        _localParticleCurrentVolume.second = value;

        _localParticleCurrentVolume.first = true;

        addLocalParticleData( &_localParticleCurrentVolume );

    }

    const floatVector* aspBase::getPreviousStateVariables( ){
        /*!
         * Get the previous state variables
         */

        return &_previousStateVariables;

    }

    const floatVector* aspBase::getPreviousLocalStateVariables( ){
        /*!
         * Get the previous state variables for the local particle
         */

        const floatVector* previousStateVariables;
        ERROR_TOOLS_CATCH( previousStateVariables = getPreviousStateVariables( ) );

        return previousStateVariables;

    }

    const floatVector* aspBase::getParticleParameters( ){
        /*!
         * Get the parameters for the internal particle computation
         */

        return &_particleParameters;

    }

    void aspBase::setLocalParticleParameters( ){
        /*!
         * Set the parameters for the local particle
         */

        setLocalParticleParameters( *getParticleParameters( ) );

    }

    void aspBase::setLocalParticleParameters( const floatVector &value ){
        /*!
         * Set the parameters for the local particle
         * 
         * \param &value: The current parameters for the local particle
         */

        _localParticleParameters.second = value;

        _localParticleParameters.first = true;

        addLocalParticleData( &_localParticleParameters );

    }


    const floatVector* aspBase::getLocalParticleParameters( ){

        if ( !_localParticleParameters.first ){

            ERROR_TOOLS_CATCH( setLocalParticleParameters( ) );

        }

        return &_localParticleParameters.second;

    }

    void aspBase::assembleLocalParticles( ){
        /*!
         * Assemble all the required quantities for the local particles
         */

        _assembledLocalParticleEnergies.second = floatVector( *getNumLocalParticles( ) );

        _assembledLocalParticleMicroCauchyStress.second = floatMatrix( *getNumLocalParticles( ) );

        _assembledLocalParticleVolumes.second = floatVector( *getNumLocalParticles( ) );

        _assembledLocalParticleLogProbabilityRatios.second = floatVector( *getNumLocalParticles( ) );

        for ( unsigned int i = 0; i < *getNumLocalParticles( ); i++ ){

            _localIndex = i; // Set the current local index

            // Quantities required for the energy calculation
            _assembledLocalParticleEnergies.second[ i ] = *getLocalParticleEnergy( );

            _assembledLocalParticleMicroCauchyStress.second[ i ] = *getLocalParticleMicroCauchyStress( );

            _assembledLocalParticleVolumes.second[ i ] = *getLocalParticleCurrentVolume( );

            _assembledLocalParticleLogProbabilityRatios.second[ i ] = *getLocalParticleLogProbabilityRatio( );

            // Quantities required for the gradient calculation

            // Quantities required for the Hessian calculation

            resetLocalParticleData( );

        }

        _assembledLocalParticleEnergies.first = true;

        _assembledLocalParticleMicroCauchyStress.first = true;

        _assembledLocalParticleVolumes.first = true;

        _assembledLocalParticleLogProbabilityRatios.first = true;

    }

    void aspBase::assembleSurfaceResponses( ){
        /*!
         * Assemble the surface responses of the particles
         */

        const unsigned int *dim = getDimension( );

        unsigned int numLocalParticles = *getNumLocalParticles( );

        unsigned int numSurfacePoints = getUnitSpherePoints( )->size( ) / ( *dim );

        _assembledSurfaceAdhesionEnergyDensities.second = std::vector< std::vector< floatVector > >( numLocalParticles, std::vector< floatVector >( numSurfacePoints, floatVector( numLocalParticles ) ) );

        _assembledSurfaceAdhesionThicknesses.second = std::vector< std::vector< floatVector > >( numLocalParticles, std::vector< floatVector >( numSurfacePoints, floatVector( numLocalParticles ) ) );

        _assembledSurfaceAdhesionTractions.second = std::vector< std::vector< floatMatrix > >( numLocalParticles, std::vector< floatMatrix >( numSurfacePoints, floatMatrix( numLocalParticles ) ) );

        _assembledSurfaceOverlapEnergyDensities.second = std::vector< std::vector< std::vector< mapFloatType > > >( numLocalParticles, std::vector< std::vector< mapFloatType > >( numSurfacePoints, std::vector< mapFloatType >( numLocalParticles ) ) );

        _assembledSurfaceOverlapThicknesses.second = std::vector< std::vector< std::vector< mapFloatType > > >( numLocalParticles, std::vector< std::vector< mapFloatType > >( numSurfacePoints, std::vector< mapFloatType >( numLocalParticles ) ) );

        _assembledSurfaceOverlapTractions.second = std::vector< std::vector< std::vector< mapFloatVector > > >( numLocalParticles, std::vector< std::vector< mapFloatVector > >( numSurfacePoints, std::vector< mapFloatVector >( numLocalParticles ) ) );

        for ( unsigned int i = 0; i < *getNumLocalParticles( ); i++ ){

            _localIndex = i; // Set the current local index

            for ( unsigned int j = 0; j < numSurfacePoints; j++ ){

                _localSurfaceNodeIndex = j; // Set the local surface node index

                for ( unsigned int k = 0; k < *getNumLocalParticles( ); k++ ){

                    _nonLocalIndex = k; // Set the interaction index

                    // Quantities required for the energy calculation
                    _assembledSurfaceAdhesionEnergyDensities.second[ i ][ j ][ k ] = *getSurfaceAdhesionEnergyDensity( );

                    _assembledSurfaceAdhesionTractions.second[ i ][ j ][ k ] = *getSurfaceAdhesionTraction( );

                    _assembledSurfaceAdhesionThicknesses.second[ i ][ j ][ k ] = *getSurfaceAdhesionThickness( );

                    _assembledSurfaceOverlapEnergyDensities.second[ i ][ j ][ k ] = *getSurfaceOverlapEnergyDensity( );

                    _assembledSurfaceOverlapTractions.second[ i ][ j ][ k ] = *getSurfaceOverlapTraction( );

                    _assembledSurfaceOverlapThicknesses.second[ i ][ j ][ k ] = *getSurfaceOverlapThickness( );

                    // Quantities required for the gradient calculation

                    // Quantities required for the Hessian calculation

                    resetInteractionPairData( );

                }

                resetSurfacePointData( );

            }

            resetLocalParticleData( );

        }

        _assembledSurfaceAdhesionEnergyDensities.first = true;

        _assembledSurfaceAdhesionTractions.first = true;

        _assembledSurfaceAdhesionThicknesses.first = true;

        _assembledSurfaceOverlapEnergyDensities.first = true;

        _assembledSurfaceOverlapTractions.first = true;

        _assembledSurfaceOverlapThicknesses.first = true;

    }

    const floatVector* aspBase::getAssembledLocalParticleEnergies( ){
        /*!
         * Get the assembled local particle energies
         */

        if ( !_assembledLocalParticleEnergies.first ){

            ERROR_TOOLS_CATCH( assembleLocalParticles( ) );

        }

        return &_assembledLocalParticleEnergies.second;

    }

    const floatMatrix* aspBase::getAssembledLocalParticleMicroCauchyStresses( ){
        /*!
         * Get the assembled local particle micro Cauchy stress
         */

        if ( !_assembledLocalParticleMicroCauchyStress.first ){

            ERROR_TOOLS_CATCH( assembleLocalParticles( ) );

        }

        return &_assembledLocalParticleMicroCauchyStress.second;

    }

    const floatVector* aspBase::getAssembledLocalParticleVolumes( ){
        /*!
         * Get the assembled local particle volumes
         */

        if ( !_assembledLocalParticleVolumes.first ){

            ERROR_TOOLS_CATCH( assembleLocalParticles( ) );

        }

        return &_assembledLocalParticleVolumes.second;

    }

    const floatVector* aspBase::getAssembledLocalParticleLogProbabilityRatios( ){
        /*!
         * Get the assembled local particle log probability ratios
         */

        if ( !_assembledLocalParticleLogProbabilityRatios.first ){

            ERROR_TOOLS_CATCH( assembleLocalParticles( ) );

        }

        return &_assembledLocalParticleLogProbabilityRatios.second;

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
