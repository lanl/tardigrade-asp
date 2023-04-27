/**
  ******************************************************************************
  * \file asp.h
  ******************************************************************************
  * A C++ library for printing messages to stdout. Used as a stub repo example.
  ******************************************************************************
  */

#ifndef ASP_H
#define ASP_H

#include<sstream>

#include<error_tools.h>
#define USE_EIGEN
#include<vector_tools.h>
#include<abaqus_tools.h>

namespace asp{

    namespace unit_test{
        class aspBaseTester;
    }

    constexpr const char* str_end(const char *str) {
        /*! Recursively search string for last character
         * \param *str: pointer to string START of UNIX path like string
         * \return *str: pointer to last character in string
         */
        return *str ? str_end(str + 1) : str;
    }
    constexpr bool str_slant(const char *str) {
        /*! Recursively search string for leftmost UNIX path separator from the left
         * \param *str: pointer to string START of UNIX path like string
         * \return bool: True if string contains UNIX path separator. Else false.
         */
        return *str == '/' ? true : (*str ? str_slant(str + 1) : false);
    }
    constexpr const char* r_slant(const char* str) {
        /*! Recursively search string for rightmost UNIX path separator from the right
         * \param *str: pointer to string END of UNIX path like string
         * \return *str: pointer to start of base name
         */
        return *str == '/' ? (str + 1) : r_slant(str - 1);
    }
    constexpr const char* file_name(const char* str) {
        /*! Return the current file name with extension at compile time
         * \param *str: pointer to string START of UNIX path like string
         * \return str: file base name
         */
        return str_slant(str) ? r_slant(str_end(str)) : str;
    }
    //Return filename for constructing debugging messages
    //https://stackoverflow.com/questions/31050113/how-to-extract-the-source-filename-without-path-and-suffix-at-compile-time
    const std::string __BASENAME__ = file_name(__FILE__);
    const std::string __FILENAME__ = __BASENAME__.substr(0, __BASENAME__.find_last_of("."));

    typedef errorTools::Node errorNode; //!< Redefinition for the error node
    typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    /// Say hello
    /// @param message The message to print
    errorOut sayHello(std::string message);

    void abaqusInterface( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                          double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                          const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                          const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                          const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                          const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                          const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                          const int *JSTEP,     const int &KINC );

    errorOut dummyMaterialModel( floatVector &stress,             floatVector &statev,        floatMatrix &ddsdde,       floatType &SSE,            floatType &SPD,
                             floatType &SCD,                  floatType &RPL,             floatVector &ddsddt,       floatVector &drplde,       floatType &DRPLDT,
                             const floatVector &strain,       const floatVector &dstrain, const floatVector &time,   const floatType &DTIME,    const floatType &TEMP,
                             const floatType &DTEMP,          const floatVector &predef,  const floatVector &dpred,  const std::string &cmname, const int &NDI,
                             const int &NSHR,                 const int &NTENS,           const int &NSTATV,         const floatVector &props,  const int &NPROPS,
                             const floatVector &coords,       const floatMatrix &drot,    floatType &PNEWDT,         const floatType &CELENT,   const floatMatrix &dfgrd0,
                             const floatMatrix &dfgrd1,       const int &NOEL,            const int &NPT,            const int &LAYER,          const int &KSPT,
                             const std::vector< int > &jstep, const int &KINC );

    class aspBase{
        /*!
         * The base class for all Anisotropic Stochastic Particle (ASP) models.
         */

        public:

            // Public parameters

            // Constructors
            aspBase( );

            // Public member functions
            virtual void computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
                                                            const floatVector &currentMicroDeformation, const floatVector &previousMicroDeformation,
                                                            const floatType &currentTemperature, const floatType &previousTemperature,
                                                            const floatVector &previousStateVariables,
                                                            const floatVector &parameters,
                                                            floatType &energy, floatVector &cauchyStress );

            virtual void computeLocalParticleEnergyDensity( const floatType &previousTime, const floatType &deltaTime,
                                                            const floatVector &currentMicroDeformation, const floatVector &previousMicroDeformation,
                                                            const floatType &currentTemperature, const floatType &previousTemperature,
                                                            const floatVector &previousStateVariables,
                                                            const floatVector &parameters,
                                                            floatType &energy, floatVector &cauchyStress, floatType &logProbabilityRatio );

            virtual void computeSurfaceAdhesionTraction( floatVector &surfaceAdhesionTraction );

            virtual void computeSurfaceAdhesionEnergyDensity( floatType &surfaceAdhesionEnergyDensity );

            virtual void computeSurfaceOverlapEnergyDensity( std::unordered_map< unsigned int, floatType > &surfaceOverlapEnergyDensity );

            // Getter functions
            const floatType* getSurfaceAdhesionEnergyDensity( );

            const floatType* getLocalReferenceRadius( );

            const floatType* getNonLocalReferenceRadius( );

            const floatVector* getLocalReferenceNormal( );

            const floatVector* getLocalSurfaceReferenceRelativePositionVector( );

            const floatVector* getNonLocalSurfaceReferenceRelativePositionVector( );

            const floatVector* getLocalDeformationGradient( );

            const floatVector* getLocalMicroDeformation( );

            const floatVector* getLocalReferenceParticleSpacing( );

            const floatVector* getMicroDeformation( );

            const floatVector* getGradientMicroDeformation( );

            const floatVector* getNonLocalMicroDeformationBase( );

            const floatVector* getNonLocalMicroDeformation( );

            const floatVector* getLocalGradientMicroDeformation( );

            const floatVector* getCurrentDistanceVector( );

            const floatVector* getLocalCurrentNormal( );

            const floatVector* getSurfaceParameters( );

            const floatVector* getSurfaceOverlapParameters( );

            const floatVector* getReferenceDistanceVector( );

            const floatVector* getUnitSpherePoints( );

            const floatVector* getLocalReferenceSurfacePoints( );

            const floatVector* getNonLocalReferenceSurfacePoints( );

            const floatVector* getLocalCurrentSurfacePoints( );

            const floatVector* getNonLocalCurrentSurfacePoints( );

            const floatVector* getSurfaceAdhesionTraction( );

            const floatMatrix* getLocalParticleCurrentBoundingBox( );

            const floatMatrix* getNonLocalParticleCurrentBoundingBox( );

            const std::unordered_map< unsigned int, floatType >* getSurfaceOverlapEnergyDensity( );

            const std::unordered_map< unsigned int, floatVector >* getParticlePairOverlap( );

            const std::vector< unsigned int >* getUnitSphereConnectivity( );

        protected:

            // Protected parameters
            unsigned int _dimension = 3;

            unsigned int _surfaceElementCount = 1;

            bool pointInBoundingBox( const floatVector &point, const floatMatrix &boundingBox );

            void formBoundingBox( const floatVector &points, floatMatrix &boundingBox );

            void idBoundingBoxContainedPoints( const floatVector &points, const floatMatrix &boundingBox, std::vector< unsigned int > &containedPoints );

        private:
            // Friend classes
            friend class unit_test::aspBaseTester;

            // Private parameters
            unsigned int _localIndex = 0;

            unsigned int _nonLocalIndex = 0;

            unsigned int _localSurfaceNodeIndex = 0;

            floatType   _radius;

            floatVector _deformationGradient;

            floatVector _microDeformation;

            floatVector _gradientMicroDeformation;

            std::pair< bool, floatMatrix > _localParticleCurrentBoundingBox;

            std::pair< bool, floatVector > _localReferenceSurfacePoints;

            std::pair< bool, floatVector > _localCurrentSurfacePoints;

            // ALL OF THESE MUST BE CLEARED AFTER EACH SURFACE INTEGRAND CALCULATION
            std::pair< bool, floatType > _localReferenceRadius;

            std::pair< bool, floatType > _nonLocalReferenceRadius;

            std::pair< bool, std::vector< unsigned int > > _unitSphereConnectivity;

            std::pair< bool, floatVector > _unitSpherePoints;

            std::pair< bool, floatVector > _localReferenceNormal;

            std::pair< bool, floatVector > _localSurfaceReferenceRelativePositionVector;

            std::pair< bool, floatVector > _nonLocalSurfaceReferenceRelativePositionVector;

            std::pair< bool, floatVector > _referenceDistanceVector;

            std::pair< bool, floatVector > _localReferenceParticleSpacing;

            std::pair< bool, floatVector > _localDeformationGradient;

            std::pair< bool, floatVector > _localMicroDeformation;

            std::pair< bool, floatVector > _nonLocalMicroDeformation;

            std::pair< bool, floatVector > _nonLocalMicroDeformationBase;

            std::pair< bool, floatVector > _localGradientMicroDeformation;

            std::pair< bool, floatVector > _currentDistanceVector;

            std::pair< bool, floatVector > _localCurrentNormal;

            std::pair< bool, floatVector > _surfaceParameters;

            std::pair< bool, floatVector > _surfaceOverlapParameters;

            std::pair< bool, floatType > _surfaceAdhesionEnergyDensity;

            std::pair< bool, std::unordered_map< unsigned int, floatType > > _surfaceOverlapEnergyDensity;

            std::pair< bool, floatVector > _nonLocalReferenceSurfacePoints;

            std::pair< bool, floatVector > _nonLocalCurrentSurfacePoints;

            std::pair< bool, floatMatrix > _nonLocalParticleCurrentBoundingBox;

            std::pair< bool, std::unordered_map< unsigned int, floatVector > > _particlePairOverlap;

            std::pair< bool, floatVector > _surfaceAdhesionTraction;
            // END OF MEMBERS WHICH MUST BE CLEARED AFTER EACH SURFACE INTEGRAND CALCULATION

            // Private member functions
            virtual void initializeUnitSphere( );

            virtual void setLocalReferenceRadius( );

            virtual void setNonLocalReferenceRadius( );

            virtual void setLocalReferenceNormal( );

            virtual void setLocalSurfaceReferenceRelativePositionVector( );

            virtual void setNonLocalSurfaceReferenceRelativePositionVector( );

            virtual void setReferenceDistanceVector( );

            virtual void setLocalDeformationGradient( );

            virtual void setLocalMicroDeformation( );

            virtual void setNonLocalMicroDeformation( );

            virtual void setNonLocalMicroDeformationBase( );

            virtual void setLocalGradientMicroDeformation( );

            virtual void setLocalCurrentNormal( );

            virtual void setLocalReferenceParticleSpacing( );

            virtual void setCurrentDistanceVector( );

            virtual void setSurfaceParameters( );

            virtual void setSurfaceOverlapParameters( );

            virtual void initializeSurfaceIntegrandQuantities( );

            virtual void setSurfaceAdhesionEnergyDensity( );

            virtual void setSurfaceOverlapEnergyDensity( );

            virtual void setSurfaceAdhesionTraction( );

            virtual void setLocalParticleCurrentBoundingBox( );

            virtual void setNonLocalParticleCurrentBoundingBox( );

            virtual void setLocalReferenceSurfacePoints( );

            virtual void setNonLocalReferenceSurfacePoints( );

            virtual void setLocalCurrentSurfacePoints( );

            virtual void setNonLocalCurrentSurfacePoints( );

            virtual void setParticlePairOverlap( );

            virtual void resetSurface( );

    };

}

#endif
