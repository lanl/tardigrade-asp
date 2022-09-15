#include<linear_elasticity.h>
#include<sstream>

namespace linearElasticity{

    /** Define the expected number of spatial dimensions */
    unsigned int spatialDimensions = 3;

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy ){
        /*!
         * Compute the value of the linear elastic energy which we define via
         * 
         * \f$\rho \psi = \frac{1}{J} \left[ \frac{\lambda}{2} \left( E_{II} \right)^2 + \mu E_{IJ} E_{JI} \right] $\f
         * 
         * \param &chi: The micro-deformation
         * \param &parameters: The parameters used in the calculation. The two Lame parameters are expected lambda and mu.
         * \param &energy: The resulting free energy in the current configuration
         */

        if ( chi.size( ) != spatialDimensions * spatialDimensions ){
            return new errorNode( __func__, "The spatial dimension of " + std::to_string( spatialDimensions ) + " is not reflected by chi which has a size of " + std::to_string( chi.size( ) ) + " rather than " + std::to_string( spatialDimensions * spatialDimensions ) );
        }

        floatVector E;
        errorOut error = constitutiveTools::computeGreenLagrangeStrain( chi, E );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        floatType detChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );

        }
        catch ( std::exception &e ) {
            std::ostringstream message;
            message << "Error in calculation of det( chi ). Error follows:\n";
            message << e.what( );
            return new errorNode( __func__, message.str( ) );

        }

        if ( parameters.size( ) != 2 ){
            return new errorNode( __func__, "The expected parameters are the Lame parameters. " + std::to_string( parameters.size( ) ) + " were provided rather than 2." );
        }

        floatType lambda = parameters[ 0 ];
        floatType mu     = parameters[ 1 ];

        floatType trE = 0;

        energy = 0;
        for ( unsigned int I = 0; I < spatialDimensions; I++ ){
            trE += E[ spatialDimensions * I + I ];
        }

        energy += ( 1. / detChi ) * 0.5 * lambda * trE * trE;

        for ( unsigned int I = 0; I < spatialDimensions; I++ ){

            for ( unsigned int J = 0; J < spatialDimensions; J++ ){

                energy += ( 1. / detChi ) * mu * E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ];

            }

        }

        return NULL;
    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress ){
        /*!
         * Compute the value of the linear elastic energy which we define via
         * 
         * \f$\rho \psi = \frac{1}{J} \left[ \frac{\lambda}{2} \left( E_{II} \right)^2 + \mu E_{IJ} E_{JI} \right] $\f
         * 
         * and the value of the Cauchy stress which is defined via
         * 
         * \f$\sigma_{ij} = \frac{1}{J} \frac{ \partial \left( \rho \psi \right ) }{\partial F_{iI}} F_{jI} $\f
         * 
         * \param &chi: The micro-deformation
         * \param &parameters: The parameters used in the calculation. The two Lame parameters are expected lambda and mu.
         * \param &energy: The resulting free energy in the current configuration
         * \param &cauchyStress; The expected cauchy stress
         */

        if ( chi.size( ) != spatialDimensions * spatialDimensions ){
            return new errorNode( __func__, "The spatial dimension of " + std::to_string( spatialDimensions ) + " is not reflected by chi which has a size of " + std::to_string( chi.size( ) ) + " rather than " + std::to_string( spatialDimensions * spatialDimensions ) );
        }

        floatVector E;
        floatMatrix dEdChi;
        errorOut error = constitutiveTools::computeGreenLagrangeStrain( chi, E, dEdChi );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        floatType detChi;
        floatVector dDetChidChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            dDetChidChi = vectorTools::computeDDetAdJ( chi, spatialDimensions, spatialDimensions );

        }
        catch ( std::exception &e ) {

            std::ostringstream message;
            message << "Error in calculation of det( chi ). Error follows:\n";
            message << e.what( );
            return new errorNode( __func__, message.str( ) );

        }

        if ( parameters.size( ) != 2 ){
            return new errorNode( __func__, "The expected parameters are the Lame parameters. " + std::to_string( parameters.size( ) ) + " were provided rather than 2." );
        }

        floatType lambda = parameters[ 0 ];
        floatType mu     = parameters[ 1 ];

        // Compute the energy and the gradient of the energy w.r.t. chi.
        floatType trE = 0;
        floatVector dtrEdchi( spatialDimensions * spatialDimensions, 0. );
        energy = 0;
        floatVector dEnergydChi( spatialDimensions * spatialDimensions, 0. );
        
        for ( unsigned int I = 0; I < spatialDimensions; I++ ){
            trE += E[ spatialDimensions * I + I ];

            for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                    dtrEdchi[ spatialDimensions * k + K ] += dEdChi[ spatialDimensions * I + I ][ spatialDimensions * k + K ];

                }

            }

        }

        energy += ( 1. / detChi ) * 0.5 * lambda * trE * trE;

        dEnergydChi = -( 1. / ( detChi * detChi ) ) * dDetChidChi * 0.5 * lambda * trE * trE +  ( 1. / detChi ) * lambda * trE * dtrEdchi;

        for ( unsigned int I = 0; I < spatialDimensions; I++ ){

            for ( unsigned int J = 0; J < spatialDimensions; J++ ){

                energy += ( 1. / detChi ) * mu * E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ];

                for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                    for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                        dEnergydChi[ spatialDimensions * k + K ] +=  -( 1. / ( detChi * detChi ) ) * dDetChidChi[ spatialDimensions * k + K ]* mu * E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ]
                                                                  + ( 1. / detChi ) * mu * dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * E[ spatialDimensions * J + I ]
                                                                  + ( 1. / detChi ) * mu * E[ spatialDimensions * I + J ] * dEdChi[ spatialDimensions * J + I ][ spatialDimensions * k + K ];

                    }

                }

            }

        }

        // Use the energy gradient to compute the Cauchy stress
        cauchyStress = floatVector( spatialDimensions * spatialDimensions, 0 );
        for ( unsigned int i = 0; i < spatialDimensions; i++ ){
            for ( unsigned int j = 0; j < spatialDimensions; j++ ){
                for ( unsigned int I = 0; I < spatialDimensions; I++ ){
                    cauchyStress[ spatialDimensions * i + j ] += dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] / detChi;
                }
            }
        }


        return NULL;
    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatVector &dCauchyStressdChi ){
        /*!
         * Compute the value of the linear elastic energy which we define via
         * 
         * \f$\rho \psi = \frac{1}{J} \left[ \frac{\lambda}{2} \left( E_{II} \right)^2 + \mu E_{IJ} E_{JI} \right] $\f
         * 
         * and the value of the Cauchy stress which is defined via
         * 
         * \f$\sigma_{ij} = \frac{1}{J} \frac{ \partial \left( \rho \psi \right ) }{\partial F_{iI}} F_{jI} $\f
         * 
         * \param &chi: The micro-deformation
         * \param &parameters: The parameters used in the calculation. The two Lame parameters are expected lambda and mu.
         * \param &energy: The resulting free energy in the current configuration
         * \param &cauchyStress; The expected cauchy stress
         * \param &dEnergydChi: The gradient of the energy w.r.t. the micro deformation
         * \param &dCauchyStressdChi: The gradient of the Cauchy stress w.r.t. the micro deformation
         * 
         */

        if ( chi.size( ) != spatialDimensions * spatialDimensions ){
            return new errorNode( __func__, "The spatial dimension of " + std::to_string( spatialDimensions ) + " is not reflected by chi which has a size of " + std::to_string( chi.size( ) ) + " rather than " + std::to_string( spatialDimensions * spatialDimensions ) );
        }

        floatVector E;
        floatMatrix dEdChi;
        errorOut error = constitutiveTools::computeGreenLagrangeStrain( chi, E, dEdChi );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        floatType detChi;
        floatVector dDetChidChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            dDetChidChi = vectorTools::computeDDetAdJ( chi, spatialDimensions, spatialDimensions );

        }
        catch ( std::exception &e ) {

            std::ostringstream message;
            message <<  "Error in calculation of det( chi ). Error follows:\n";
            message << e.what( );
            return new errorNode( __func__, message.str( ) );

        }

        if ( parameters.size( ) != 2 ){
            return new errorNode( __func__, "The expected parameters are the Lame parameters. " + std::to_string( parameters.size( ) ) + " were provided rather than 2." );
        }

        floatType lambda = parameters[ 0 ];
        floatType mu     = parameters[ 1 ];

        // Compute the energy and the gradient of the energy w.r.t. chi.
        energy = 0;
        dEnergydChi = floatVector( spatialDimensions * spatialDimensions, 0. );
        
        for ( unsigned int I = 0; I < spatialDimensions; I++ ){

            energy += ( 1. / detChi ) * 0.5 * lambda * E[ spatialDimensions * I + I ];

            for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                    dEnergydChi[ spatialDimensions * k + K ] += -( 1. / ( detChi * detChi ) ) * dDetChidChi[ spatialDimensions * k + K ] * 0.5 * lambda * E[ spatialDimensions * I + I ]
                                                              +   ( 1. / detChi ) * 0.5 * lambda * dEdChi[ spatialDimensions * I + I ][ spatialDimensions * k + K ];

                }

            }

            for ( unsigned int J = 0; J < spatialDimensions; J++ ){

                energy += ( 1. / detChi ) * mu * E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ];

                for ( unsigned int k = 0; k < spatialDimensions; k++ ){
    
                    for ( unsigned int K = 0; K < spatialDimensions; K++ ){
    
                        dEnergydChi[ spatialDimensions * k + K ] += -( 1. / ( detChi * detChi ) ) * dDetChidChi[ spatialDimensions * k + K ] * mu * E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ]
                                                                  +  ( 1. / detChi ) * mu * dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * E[ spatialDimensions * J + I ]
                                                                  +  ( 1. / detChi ) * mu * E[ spatialDimensions * I + J ] * dEdChi[ spatialDimensions * J + I ][ spatialDimensions * k + K ];

    
                    }
    
                }

            }

        }

        // Use the energy gradient to compute the Cauchy stress
        cauchyStress = floatVector( spatialDimensions * spatialDimensions, 0 );
        for ( unsigned int i = 0; i < spatialDimensions; i++ ){
            for ( unsigned int j = 0; j < spatialDimensions; j++ ){
                for ( unsigned int I = 0; I < spatialDimensions; I++ ){
                    cauchyStress[ spatialDimensions * i + j ] += dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] / detChi;
                }
            }
        }


        return NULL;
    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatVector &dCauchyStressdChi,
                             floatVector &d2EnergydChi2, floatVector &d2CauchyStressdChi2 );

}

