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

    void aspBase::setLocalReferenceNormal( ){
        /*!
         * Set the local reference normal vector
         * 
         * Because we have a unit sphere centered on zero as a basis the normal vector is the same
         * as any position vector on the surface of the sphere.
         */

        if ( ! _unitSpherePoints.first ){

            ERROR_TOOLS_CATCH( initializeUnitSphere( ) );

        }

        if ( _dimension * ( _localIndex + 1 ) > _unitSpherePoints.second.size( ) ){

            std::string message = "The requested index is greater than the number of points available on the unit sphere.\n";
            message            += "  localSurfaceNodeIndex: " + std::to_string( _localSurfaceNodeIndex ) + "\n";
            message            += "  number of points:      " + std::to_string( _unitSpherePoints.second.size( ) / _dimension );

            ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        _localReferenceNormal.second = floatVector( _unitSpherePoints.second.begin( ) + _dimension * _localSurfaceNodeIndex,
                                                    _unitSpherePoints.second.begin( ) + _dimension * ( _localSurfaceNodeIndex + 1 ) );

        _localReferenceNormal.second /= vectorTools::l2norm( _localReferenceNormal.second );

        _localReferenceNormal.first = true;

        return;

    }

    floatVector aspBase::getLocalReferenceNormal( ){
        /*!
         * Get the local reference normal vector
         * 
         * Get the local reference normal vector or set it and return it if necessary
         */

        if ( !_localReferenceNormal.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceNormal( ) );

        }

        return _localReferenceNormal.second;

    }

    void aspBase::setLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Set the local surface relative position vector
         * 
         * \f$\Xi_I = R^{local} N_I\f$
         * 
         * where \f$R^{local}\f$ is the local radius and \f$N_I\f$ is the local reference normal.
         */

        floatVector localReferenceNormal;
        ERROR_TOOLS_CATCH( localReferenceNormal = getLocalReferenceNormal( ) );

        if ( ! _localReferenceRadius.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceRadius ( ) );

        }

        _localSurfaceReferenceRelativePositionVector.second = _localReferenceRadius.second * localReferenceNormal;

        _localSurfaceReferenceRelativePositionVector.first = true;

        return;

    }

    floatVector aspBase::getLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Get the value of the local surface reference relative position vector
         */

        if ( !_localSurfaceReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setLocalSurfaceReferenceRelativePositionVector( ) );

        }

        return _localSurfaceReferenceRelativePositionVector.second;

    }

    void aspBase::setNonLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Set the non-local relative position vector
         * 
         * \f$\Xi_I^{non-local} = -R^{non-local} N_I\f$
         * 
         * where \f$R^{non-local}\f$ is the non-local radius and \f$N_I\f$ is the local reference normal.
         */

        floatVector localReferenceNormal;
        ERROR_TOOLS_CATCH( localReferenceNormal = getLocalReferenceNormal( ) );

        if ( ! _nonlocalReferenceRadius.first ){

            ERROR_TOOLS_CATCH( setNonLocalReferenceRadius( ) );

        }

        _nonlocalSurfaceReferenceRelativePositionVector.second = -_nonlocalReferenceRadius.second * localReferenceNormal;

        _nonlocalSurfaceReferenceRelativePositionVector.first = true;

        return;

    }

    floatVector aspBase::getNonLocalSurfaceReferenceRelativePositionVector( ){
        /*!
         * Get the non-local surface reference relative position vector
         */

        if ( !_nonlocalSurfaceReferenceRelativePositionVector.first ){

            ERROR_TOOLS_CATCH( setNonLocalSurfaceReferenceRelativePositionVector( ) );

        }

        return _nonlocalSurfaceReferenceRelativePositionVector.second;

    }

    void aspBase::setLocalReferenceRadius( ){
        /*!
         * Set the local reference radius
         */

        _localReferenceRadius.second = _radius;

        _localReferenceRadius.first = true;

        return;

    }

    void aspBase::setNonLocalReferenceRadius( ){
        /*!
         * Set the non-local reference radius
         */

        _nonlocalReferenceRadius.second = _radius;

        _nonlocalReferenceRadius.first = true;

        return;

    }

    void aspBase::setLocalDeformationGradient( ){
        /*!
         * Set the local deformation gradient
         */

        _localDeformationGradient.second = _deformationGradient;

        _localDeformationGradient.first = true;

        return;

    }

    void aspBase::setLocalMicroDeformation( ){
        /*!
         * Set the local micro deformation
         */

        _localMicroDeformation.second = _microDeformation;

        _localMicroDeformation.first = true;

        return;

    }

    void aspBase::setLocalReferenceParticleSpacing( ){
        /*!
         * Set the local particle spacing
         * 
         * \f$dX_I = \Xi_I^{local} + D_I - \Xi_I^{non-local}\f$
         */

        floatVector localSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( localSurfaceReferenceRelativePositionVector = getLocalSurfaceReferenceRelativePositionVector( ) );

        floatVector nonlocalSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( nonlocalSurfaceReferenceRelativePositionVector = getNonLocalSurfaceReferenceRelativePositionVector( ) );

        if ( ! _referenceDistanceVector.first ){

            ERROR_TOOLS_CATCH( setReferenceDistanceVector( ) );

        }

        _localReferenceParticleSpacing.second = localSurfaceReferenceRelativePositionVector
                                              + _referenceDistanceVector.second
                                              - nonlocalSurfaceReferenceRelativePositionVector;

        _localReferenceParticleSpacing.first = true;

        return;

    }

    void aspBase::setNonLocalMicroDeformation( ){
        /*!
         * Set the non-local micro deformation
         * 
         * \f$ \chi_{iI}^{non-local} = \chi_{iI} + \chi_{iI,J} dX_J\f$
         */

        if ( ! _localReferenceParticleSpacing.first ){

            ERROR_TOOLS_CATCH( setLocalReferenceParticleSpacing( ) );

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

        return;

    }

    void aspBase::setCurrentDistance( ){
        /*!
         * Set the current distance vector
         */

        floatVector localSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( localSurfaceReferenceRelativePositionVector = getLocalSurfaceReferenceRelativePositionVector( ) );

        floatVector nonlocalSurfaceReferenceRelativePositionVector;
        ERROR_TOOLS_CATCH( nonlocalSurfaceReferenceRelativePositionVector = getNonLocalSurfaceReferenceRelativePositionVector( ) );

        if ( ! _referenceDistanceVector.first ){

            ERROR_TOOLS_CATCH( setReferenceDistanceVector( ) );

        }

        if ( ! _localDeformationGradient.first ){

            ERROR_TOOLS_CATCH( setLocalDeformationGradient( ) );

        }

        if ( ! _localMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setLocalMicroDeformation( ) );

        }

        if ( ! _nonlocalMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setNonLocalMicroDeformation( ) );

        }

        // Compute the current distance
        ERROR_TOOLS_CATCH( tractionSeparation::computeCurrentDistanceGeneral( localSurfaceReferenceRelativePositionVector,
                                                                              nonlocalSurfaceReferenceRelativePositionVector,
                                                                              _referenceDistanceVector.second,
                                                                              _localDeformationGradient.second,
                                                                              _localMicroDeformation.second,
                                                                              _nonlocalMicroDeformation.second,
                                                                              _currentDistanceVector.second ) );

        _currentDistanceVector.first = true;

        return;

    }

    void aspBase::setLocalCurrentNormal( ){
        /*!
         * Set the current local normal
         */

        floatVector localReferenceNormal;
        ERROR_TOOLS_CATCH( localReferenceNormal = getLocalReferenceNormal( ) );

        if ( ! _localMicroDeformation.first ){

            ERROR_TOOLS_CATCH( setLocalMicroDeformation( ) );

        }

        // Compute the current normal
        ERROR_TOOLS_CATCH( tractionSeparation::computeNansonsRelation( _localMicroDeformation.second, localReferenceNormal,
                                                                       _localCurrentNormal.second ) );

        _localCurrentNormal.second /= vectorTools::l2norm( _localCurrentNormal.second );

        _localCurrentNormal.first = true;

        return;

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

        ERROR_TOOLS_CATCH( setCurrentDistance( ) );

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

    void aspBase::setReferenceDistanceVector( ){
        /*!
         * Set the initial distance between two particles
         */

        _referenceDistanceVector.second = floatVector( _dimension, 0. );

        _referenceDistanceVector.first = true;

        return;

    }

    void aspBase::resetSurface( ){
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

        return;

    }

    void aspBase::computeSurfaceEnergyDensity( floatType &surfaceEnergyDensity ){
        /*!
         * Compute the surface energy density in the current configuration ( energy / da )
         * 
         * It is expected that the user will use the function `setSurfaceEnergyDensity` to
         * set the current value of the surface energy density.
         */

        // Decompose the traction into normal and tangential directions
        floatVector dn, dt;
        ERROR_TOOLS_CATCH_NODE_POINTER( tractionSeparation::decomposeVector( _currentDistanceVector.second, _localCurrentNormal.second, dn, dt ) );

        floatType energyDensity;
        ERROR_TOOLS_CATCH_NODE_POINTER( tractionSeparation::computeLinearTractionEnergy( dn, dt, _surfaceParameters.second, energyDensity ) );

        surfaceEnergyDensity = 0.5 * energyDensity * vectorTools::l2norm( dn );

        return;

    }

    floatType aspBase::getSurfaceEnergyDensity( ){
        /*!
         * Get the surface energy density
         */

        if ( !_surfaceEnergyDensity.first ){

            ERROR_TOOLS_CATCH( setSurfaceEnergyDensity( ) );

        }

        return _surfaceEnergyDensity.second;

    }

    void aspBase::setSurfaceEnergyDensity( ){
        /*!
         * Set the surface energy density if required.
         */

        ERROR_TOOLS_CATCH( computeSurfaceEnergyDensity( _surfaceEnergyDensity.second ) );

        _surfaceEnergyDensity.first = true;

        return;

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
