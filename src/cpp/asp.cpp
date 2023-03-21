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

        _nonlocalReferenceRadius = std::make_pair( false, 0 );

        _unitSpherePoints = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localReferenceNormal = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _nonlocalSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _referenceDistanceVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localReferenceParticleSpacing = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localDeformationGradient = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localMicroDeformation = std::make_pair( false, floatVector( _dimension, 0 ) );

        _nonlocalMicroDeformation = std::make_pair( false, floatVector( _dimension, 0 ) );

        _currentDistanceVector = std::make_pair( false, floatVector( _dimension, 0 ) );

        _localCurrentNormal = std::make_pair( false, floatVector( _dimension, 0 ) );

    }

    errorOut aspBase::computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
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

        errorOut error = stressTools::linearElasticity::evaluateEnergy( currentMicroDeformation, parameters, energy, cauchyStress );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in linear elasticity" );

            result->addNext( error );

        }

        return NULL;

    }

    errorOut aspBase::computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
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

        errorOut error = stressTools::linearElasticity::evaluateEnergy( currentMicroDeformation, parameters, energy, cauchyStress );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in linear elasticity" );

            result->addNext( error );

        }

        logProbabilityRatio = 0.;

        return NULL;

    }

    errorOut aspBase::initializeUnitSphere( ){
        /*!
         * Initialize the unit sphere (i.e. a sphere of radius 1) to integrate over
         */

        surfaceIntegration::decomposeSphere( 1.0, _surfaceElementCount, _unitSpherePoints.second, _unitSphereConnectivity.second );

        _unitSpherePoints.first = true;

        _unitSphereConnectivity.first = true;

        return NULL;

    }

    errorOut aspBase::setLocalReferenceNormal( ){
        /*!
         * Set the local reference normal vector
         * 
         * Because we have a unit sphere centered on zero as a basis the normal vector is the same
         * as any position vector on the surface of the sphere.
         */

        errorOut error;

        if ( ! _unitSpherePoints.first ){

            error = initializeUnitSphere( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when initializing the unit sphere" );

                result->addNext( error );

                return result;

            }

        }

        if ( _dimension * ( _localIndex + 1 ) > _unitSpherePoints.second.size( ) ){

            std::string message = "The requested index is greater than the number of points available on the unit sphere.\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";
            message            += "  number of points:      " + std::to_string( _unitSpherePoints.second.size( ) / _dimension );

            return new errorNode( __func__, message );

        }

        _localReferenceNormal.second = floatVector( _unitSpherePoints.second.begin( ) + _dimension * _localSurfaceNodeIndex,
                                                    _unitSpherePoints.second.begin( ) + _dimension * ( _localSurfaceNodeIndex + 1 ) );

        _localReferenceNormal.second /= vectorTools::l2norm( _localReferenceNormal.second );

        _localReferenceNormal.first = true;

        return NULL;

    }

    errorOut aspBase::setLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Set the local surface relative position vector
         * 
         * \f$\Xi_I = R^{local} N_I\f$
         * 
         * where \f$R^{local}\f$ is the local radius and \f$N_I\f$ is the local reference normal.
         */

        errorOut error;

        if ( ! _localReferenceNormal.first ){

            error = setLocalReferenceNormal( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when initializing the local reference normal vector" );

                result->addNext( error );

                return result;

            }

        }

        if ( ! _localReferenceRadius.first ){

            error = setLocalReferenceRadius ( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when initializing the local reference radius" );

                result->addNext( error );

                return result;

            }

        }

        _localSurfaceReferenceRelativePositionVector.second = _localReferenceRadius.second * _localReferenceNormal.second;

        _localSurfaceReferenceRelativePositionVector.first = true;

        return NULL;

    }

    errorOut aspBase::setNonLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Set the non-local relative position vector
         * 
         * \f$\Xi_I^{non-local} = -R^{non-local} N_I\f$
         * 
         * where \f$R^{non-local}\f$ is the non-local radius and \f$N_I\f$ is the local reference normal.
         */

        errorOut error;

        if ( ! _localReferenceNormal.first ){

            error = setLocalReferenceNormal( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when initializing the local reference normal vector" );

                result->addNext( error );

                return result;

            }

        }

        if ( ! _nonlocalReferenceRadius.first ){

            error = setNonLocalReferenceRadius( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when initializing the non-local reference radius" );

                result->addNext( error );

                return result;

            }

        }

        _nonlocalSurfaceReferenceRelativePositionVector.second = -_nonlocalReferenceRadius.second * _localReferenceNormal.second;

        _nonlocalSurfaceReferenceRelativePositionVector.first = true;

        return NULL;

    }

    errorOut aspBase::setLocalReferenceRadius( ){
        /*!
         * Set the local reference radius
         */

        _localReferenceRadius.second = _radius;

        _localReferenceRadius.first = true;

        return NULL;

    }

    errorOut aspBase::setNonLocalReferenceRadius( ){
        /*!
         * Set the non-local reference radius
         */

        _nonlocalReferenceRadius.second = _radius;

        _nonlocalReferenceRadius.first = true;

        return NULL;

    }

    errorOut aspBase::setLocalDeformationGradient( ){
        /*!
         * Set the local deformation gradient
         */

        _localDeformationGradient.second = _deformationGradient;

        _localDeformationGradient.first = true;

        return NULL;

    }

    errorOut aspBase::setLocalMicroDeformation( ){
        /*!
         * Set the local micro deformation
         */

        _localMicroDeformation.second = _microDeformation;

        _localMicroDeformation.first = true;

        return NULL;

    }

    errorOut aspBase::setLocalReferenceParticleSpacing( ){
        /*!
         * Set the local particle spacing
         * 
         * \f$dX_I = \Xi_I^{local} + D_I - \Xi_I^{non-local}\f$
         */

        errorOut error;

        if ( ! _localSurfaceReferenceRelativePositionVector.first ){

            error = setLocalSurfaceReferenceRelativePositionVector( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error while building the local surface relative position vector" );

                result->addNext( error );

                return result;

            }

        }

        if ( ! _nonlocalSurfaceReferenceRelativePositionVector.first ){

            error = setNonLocalSurfaceReferenceRelativePositionVector( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error while building the non-local surface relative position vector" );

                result->addNext( error );

                return result;

            }

        }

        if ( ! _referenceDistanceVector.first ){

            error = setReferenceDistanceVector( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error while building the reference distance vector" );

                result->addNext( error );

                return result;

            }

        }

        _localReferenceParticleSpacing.second = _localSurfaceReferenceRelativePositionVector.second
                                              + _referenceDistanceVector.second
                                              - _nonlocalSurfaceReferenceRelativePositionVector.second;

        _localReferenceParticleSpacing.first = true;

        return NULL;

    }

    errorOut aspBase::setNonLocalMicroDeformation( ){
        /*!
         * Set the non-local micro deformation
         * 
         * \f$ \chi_{iI}^{non-local} = \chi_{iI} + \chi_{iI,J} dX_J\f$
         */

        errorOut error;

        if ( ! _localReferenceParticleSpacing.first ){

            error = setLocalReferenceParticleSpacing( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the refernece relative particle spacing" );

                result->addNext( error );

                return result; 

            }

        }

        _nonlocalMicroDeformation.second = _microDeformation;

        for ( unsigned int i = 0; i < _dimension; i++ ){

            for ( unsigned int I = 0; I < _dimension; I++ ){

                for ( unsigned int J = 0; J < _dimension; J++ ){

                    _nonlocalMicroDeformation.second[ _dimension * i + I ]
                        += _gradientMicroDeformation[ _dimension * _dimension * i + _dimension * I + J ]
                         * _localReferenceParticleSpacing.second[ J ];

                }

            }

        }

        _nonlocalMicroDeformation.first = true;

        return NULL;

    }

    errorOut aspBase::setCurrentDistance( ){
        /*!
         * Set the current distance vector
         */

        errorOut error;

        if ( ! _localSurfaceReferenceRelativePositionVector.first ){

            error = setLocalSurfaceReferenceRelativePositionVector( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the local reference surface relative position vector" );
    
                result->addNext( error );
    
                return result;

            }

        }

        if ( ! _nonlocalSurfaceReferenceRelativePositionVector.first ){

            error = setNonLocalSurfaceReferenceRelativePositionVector( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the non-local reference surface relative position vector" );
    
                result->addNext( error );
    
                return result;

            }

        }

        if ( ! _referenceDistanceVector.first ){

            error = setReferenceDistanceVector( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing distance vector in the reference configuration" );
    
                result->addNext( error );
    
                return result;

            }

        }

        if ( ! _localDeformationGradient.first ){

            error = setLocalDeformationGradient( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the local deformation gradient" );
    
                result->addNext( error );
    
                return result;

            }

        }

        if ( ! _localMicroDeformation.first ){

            error = setLocalMicroDeformation( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the local micro deformation" );
    
                result->addNext( error );
    
                return result;

            }

        }

        if ( ! _nonlocalMicroDeformation.first ){

            error = setNonLocalMicroDeformation( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the non-local micro deformation" );
    
                result->addNext( error );
    
                return result;

            }

        }

        // Compute the current distance
        error = tractionSeparation::computeCurrentDistanceGeneral( _localSurfaceReferenceRelativePositionVector.second,
                                                                   _nonlocalSurfaceReferenceRelativePositionVector.second,
                                                                   _referenceDistanceVector.second,
                                                                   _localDeformationGradient.second,
                                                                   _localMicroDeformation.second,
                                                                   _nonlocalMicroDeformation.second,
                                                                   _currentDistanceVector.second );

        if ( error ){

            std::string message = "Error when computing the current distance";
            
            errorOut result = new errorNode( __func__, message );

            return result;

        }

        _currentDistanceVector.first = true;

        return NULL;

    }

    errorOut aspBase::setLocalCurrentNormal( ){
        /*!
         * Set the current local normal
         */

        errorOut error;

        if ( ! _localReferenceNormal.first ){

            error = setLocalReferenceNormal( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the local reference normal" );

                result->addNext( error );

                return result;

            }

        }

        if ( ! _localMicroDeformation.first ){

            error = setLocalMicroDeformation( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when computing the local micro deformation" );

                result->addNext( error );

                return result;

            }

        }

        // Compute the current normal
        error = tractionSeparation::computeNansonsRelation( _localMicroDeformation.second, _localReferenceNormal.second,
                                                            _localCurrentNormal.second );

        if ( error ){

            std::string message = "Error when computing the current normal";
            
            errorOut result = new errorNode( __func__, message );

            return result;

        }

        _localCurrentNormal.second /= vectorTools::l2norm( _localCurrentNormal.second );

        _localCurrentNormal.first = true;

        return NULL;

    }

    errorOut aspBase::initializeSurfaceIntegrandQuantities( ){
        /*!
         * Initialize the surface integrand quantities
         */

        errorOut error;

        error = setLocalReferenceNormal( );

        if ( error ){

            std::string message = "Error in setting the local reference normal\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setLocalSurfaceReferenceRelativePositionVector( );

        if ( error ){

            std::string message = "Error in setting the local surface relative position vector\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setNonLocalSurfaceReferenceRelativePositionVector( );

        if ( error ){

            std::string message = "Error in setting the non-local surface relative position vector\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  nonlocalIndex:         " + std::to_string( _nonlocalIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setReferenceDistanceVector( );

        if ( error ){

            std::string message = "Error in setting the reference distance vector\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  nonlocalIndex:         " + std::to_string( _nonlocalIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setLocalDeformationGradient( );

        if ( error ){

            std::string message = "Error in setting the local deformation gradient\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setLocalMicroDeformation( );

        if ( error ){

            std::string message = "Error in setting the local micro-deformation tensor\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setNonLocalMicroDeformation( );

        if ( error ){

            std::string message = "Error in setting the non-local micro-deformation tensor\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  nonlocalIndex:         " + std::to_string( _nonlocalIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setCurrentDistance( );

        if ( error ){

            std::string message = "Error in setting the current distance\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  nonlocalIndex:         " + std::to_string( _nonlocalIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setLocalCurrentNormal( );

        if ( error ){

            std::string message = "Error in setting the local normal vector\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        error = setSurfaceParameters( );

        if ( error ){

            std::string message = "Error in setting the surface parameters vector\n";
            message            += "  localIndex:            " + std::to_string( _localIndex ) + "\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";

            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        return NULL;

    }

    errorOut aspBase::setSurfaceParameters( ){
        /*!
         * Set the surface parameters for the integrand
         */

        _surfaceParameters.first = true;

        return NULL;
    }

    errorOut aspBase::setReferenceDistanceVector( ){
        /*!
         * Set the initial distance between two particles
         */

        _referenceDistanceVector.second = floatVector( _dimension, 0. );

        _referenceDistanceVector.first = true;

        return NULL;

    }

    errorOut aspBase::resetSurface( ){
        /*!
         * Reset the surface to the base state
         */

        _localReferenceRadius = std::make_pair( false, 0. );

        _nonlocalReferenceRadius = std::make_pair( false, 0. );

        _localReferenceNormal = std::make_pair( false, floatVector( 0, 0 ) );

        _localSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( 0, 0 ) );

        _nonlocalSurfaceReferenceRelativePositionVector = std::make_pair( false, floatVector( 0, 0 ) );

        _referenceDistanceVector = std::make_pair( false, floatVector( 0, 0 ) );

        _localReferenceParticleSpacing = std::make_pair( false, floatVector( 0, 0 ) );

        _localDeformationGradient = std::make_pair( false, floatVector( 0, 0 ) );

        _localMicroDeformation = std::make_pair( false, floatVector( 0, 0 ) );

        _nonlocalMicroDeformation = std::make_pair( false, floatVector( 0, 0 ) );

        _currentDistanceVector = std::make_pair( false, floatVector( 0, 0 ) );

        _localCurrentNormal = std::make_pair( false, floatVector( 0, 0 ) );

        _surfaceParameters = std::make_pair( false, floatVector( 0, 0 ) );

        return NULL;

    }

    errorOut aspBase::computeSurfaceEnergyDensity( ){
        /*!
         * Compute the surface energy density in the current configuration ( energy / da )
         * 
         * It is expected that the user will use the function `setSurfaceEnergyDensity` to
         * set the current value of the surface energy density.
         */

        errorOut error;

        // Decompose the traction into normal and tangential directions
        floatVector dn, dt;
        error = tractionSeparation::decomposeVector( _currentDistanceVector.second, _localCurrentNormal.second, dn, dt );

        if ( error ){

            std::string message = "Error when decomposing the distance";
            
            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        floatType energyDensity;
        error = tractionSeparation::computeLinearTractionEnergy( dn, dt, _surfaceParameters.second, energyDensity );

        if ( error ){

            std::string message = "Error when computing the surface traction energy density";
            
            errorOut result = new errorNode( __func__, message );

            result->addNext( error );

            return result;

        }

        setSurfaceEnergyDensity( energyDensity * vectorTools::l2norm( dn ) );

        return NULL;

    }

    floatType aspBase::getSurfaceEnergyDensity( ){
        /*!
         * Get the surface energy density
         */

        errorOut error = NULL;

        if ( !_surfaceEnergyDensity.first ){

            error = setSurfaceEnergyDensity( );

        }

        //TODO: Catch this error using boost::stacktrace

        return _surfaceEnergyDensity.second;

    }

    errorOut aspBase::setSurfaceEnergyDensity( ){
        /*!
         * Set the surface energy density if required.
         */

        errorOut error;

        if ( !_surfaceEnergyDensity.first ){

            error = computeSurfaceEnergyDensity( );

        }

        return error;

    }

    errorOut aspBase::setSurfaceEnergyDensity( const floatType &surfaceEnergyDensity ){
        /*!
         * Set the surface energy density with a new value
         * 
         * \param &surfaceEnergyDensity: The value of the surface energy density in the current configuration ( e / da )
         */

        _surfaceEnergyDensity.first = true;

        _surfaceEnergyDensity.second = surfaceEnergyDensity;

        return NULL;

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
