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

    class dataBase{

        public:

            virtual void clear( ){

                ERROR_TOOLS_CATCH( throw std::runtime_error( "clear not implemented!" ) );

            }

    };

    template < typename T >
    class dataStorage : public dataBase{

        public:

            bool first = false;

            T second;

            dataStorage( ){ };

            dataStorage( const bool &_first, const T &_second ) : first( _first ), second( _second ) { }

            virtual void clear( ){

                first = false;

                second.clear( );

            }

    };

    template <>
    inline void dataStorage< int >::clear( ){

        first = false;

        second = 0;

    }

    template <>
    inline void dataStorage< unsigned int >::clear( ){

        first = false;

        second = 0;

    }

    template <>
    inline void dataStorage< floatType >::clear( ){

        first = false;

        second = 0;

    }

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

            virtual void computeSurfaceOverlapTraction( std::unordered_map< unsigned int, floatVector > &surfaceOverlapTraction );

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

            const std::unordered_map< unsigned int, floatVector>* getSurfaceOverlapTraction( );

            const floatMatrix* getLocalParticleCurrentBoundingBox( );

            const floatMatrix* getNonLocalParticleCurrentBoundingBox( );

            const std::unordered_map< unsigned int, floatType >* getSurfaceOverlapEnergyDensity( );

            const std::unordered_map< unsigned int, floatVector >* getParticlePairOverlap( );

            const std::vector< unsigned int >* getUnitSphereConnectivity( );

            // Add functions
            void addLocalParticleData( dataBase *data ){ _localParticleData.push_back( data ); }

            void addSurfacePointData( dataBase *data ){ _surfacePointData.push_back( data ); }

            void addInteractionPairData( dataBase *data ){ _interactionPairData.push_back( data ); }

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

            dataStorage< floatMatrix > _localParticleCurrentBoundingBox;

            dataStorage< floatVector > _localReferenceSurfacePoints;

            dataStorage< floatVector > _localCurrentSurfacePoints;

            // ALL OF THESE MUST BE CLEARED AFTER EACH SURFACE INTEGRAND CALCULATION
            dataStorage< floatType > _localReferenceRadius;

            dataStorage< floatType > _nonLocalReferenceRadius;

            dataStorage< std::vector< unsigned int > > _unitSphereConnectivity;

            dataStorage< floatVector > _unitSpherePoints;

            dataStorage< floatVector > _localReferenceNormal;

            dataStorage< floatVector > _localSurfaceReferenceRelativePositionVector;

            dataStorage< floatVector > _nonLocalSurfaceReferenceRelativePositionVector;

            dataStorage< floatVector > _referenceDistanceVector;

            dataStorage< floatVector > _localReferenceParticleSpacing;

            dataStorage< floatVector > _localDeformationGradient;

            dataStorage< floatVector > _localMicroDeformation;

            dataStorage< floatVector > _nonLocalMicroDeformation;

            dataStorage< floatVector > _nonLocalMicroDeformationBase;

            dataStorage< floatVector > _localGradientMicroDeformation;

            dataStorage< floatVector > _currentDistanceVector;

            dataStorage< floatVector > _localCurrentNormal;

            dataStorage< floatVector > _surfaceParameters;

            dataStorage< floatVector > _surfaceOverlapParameters;

            dataStorage< floatType > _surfaceAdhesionEnergyDensity;

            dataStorage< std::unordered_map< unsigned int, floatType > > _surfaceOverlapEnergyDensity;

            dataStorage< floatVector > _nonLocalReferenceSurfacePoints;

            dataStorage< floatVector > _nonLocalCurrentSurfacePoints;

            dataStorage< floatMatrix > _nonLocalParticleCurrentBoundingBox;

            dataStorage< std::unordered_map< unsigned int, floatVector > > _particlePairOverlap;

            dataStorage< floatVector > _surfaceAdhesionTraction;

            dataStorage< std::unordered_map< unsigned int, floatVector > > _surfaceOverlapTraction;

            dataStorage< floatVector > _allParticleSurfaceAdhesionEnergy;

            dataStorage< floatMatrix > _allParticleSurfaceAdhesionTraction;

            dataStorage< std::vector< std::unordered_map< unsigned int, floatType > > > _allParticleSurfaceOverlapEnergy;

            dataStorage< std::vector< std::unordered_map< unsigned int, floatVector > > > _allParticleSurfaceOverlapTraction;

            dataStorage< floatVector > _allParticleSurfaceConstraintEnergy;

            std::vector< dataBase* > _localParticleData; //! A vector of pointers to quantities required for a local particle

            std::vector< dataBase* > _surfacePointData; //! A vector of pointers to quantities required for a local surface point

            std::vector< dataBase* > _interactionPairData; //! A vector of pointers to quantities required for a particle interaction

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

            virtual void setSurfaceAdhesionEnergyDensity( );

            virtual void setSurfaceOverlapEnergyDensity( );

            virtual void setSurfaceAdhesionTraction( );

            virtual void setSurfaceOverlapTraction( );

            virtual void setLocalParticleCurrentBoundingBox( );

            virtual void setNonLocalParticleCurrentBoundingBox( );

            virtual void setLocalReferenceSurfacePoints( );

            virtual void setNonLocalReferenceSurfacePoints( );

            virtual void setLocalCurrentSurfacePoints( );

            virtual void setNonLocalCurrentSurfacePoints( );

            virtual void setParticlePairOverlap( );

            virtual void resetInteractionPairData( );

            virtual void resetSurfacePointData( );

            virtual void resetLocalParticleData( );

    };

}

#endif
