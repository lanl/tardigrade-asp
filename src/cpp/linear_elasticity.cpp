#include<linear_elasticity.h>
#include<sstream>

namespace linearElasticity{

    /** Define the expected number of spatial dimensions */
    unsigned int spatialDimensions = 3;

    errorOut formReferenceStiffnessTensor( const floatVector &parameters, floatMatrix &C ){
        /*!
         * Form the stiffness tensor in the reference configuration.
         * 
         * \param &parameters: The parameters. The first index is the type of stiffness tensor and the later values are the coefficients.
         *     -type 0: Isotropic stiffness parameterized by lambda and mu
         * \param &C: The resulting stiffness tensor.
         */

        if ( parameters.size( ) < 1 ){

            return new errorNode( __func__, "The parameters must at least define the type" );

        }

        // Form isotropic stiffness tensor
        if ( int( parameters[ 0 ] + 0.5 ) == 0 ){

            if ( parameters.size( ) != 3 ){

                return new errorNode( __func__, "Isotropic stiffeness requires two parameters. Parameters only defines " + std::to_string( parameters.size( ) - 1 ) );

            }

            floatType lambda = parameters[ 1 ];
            floatType mu     = parameters[ 2 ];

            C = {
                    { lambda + 2 * mu,      0,      0,      0,          lambda,      0,      0,      0,          lambda },
                    {               0,      0,      0, 2 * mu,               0,      0,      0,      0,               0 },
                    {               0,      0,      0,      0,               0,      0, 2 * mu,      0,               0 },
                    {               0, 2 * mu,      0,      0,               0,      0,      0,      0,               0 },
                    {          lambda,      0,      0,      0, lambda + 2 * mu,      0,      0,      0,          lambda },
                    {               0,      0,      0,      0,               0,      0,      0, 2 * mu,               0 },
                    {               0,      0, 2 * mu,      0,               0,      0,      0,      0,               0 },
                    {               0,      0,      0,      0,               0, 2 * mu,      0,      0,               0 },
                    {          lambda,      0,      0,      0,          lambda,      0,      0,      0, lambda + 2 * mu }
                }; 

        }
        else{

            return new errorNode( __func__, "Type " + std::to_string( parameters[ 0 ] ) + " is not recognized" );

        }

        return NULL;

    }

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
        floatType trEsq = 0;

        energy = 0;
        for ( unsigned int I = 0; I < spatialDimensions; I++ ){
            trE += E[ spatialDimensions * I + I ];
            for ( unsigned int J = 0; J < spatialDimensions; J++ ){
                trEsq += E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ];
            }
        }

        energy += ( 1. / detChi ) * ( 0.5 * lambda * trE * trE + mu * trEsq );

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
        floatVector invChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            invChi = vectorTools::inverse( chi, spatialDimensions, spatialDimensions );

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
        floatType trEsq = 0;
        floatVector dtrEdChi( spatialDimensions * spatialDimensions, 0. );
        floatVector dtrEsqdChi( spatialDimensions * spatialDimensions, 0. );
        energy = 0;
        floatVector dEnergydChi( spatialDimensions * spatialDimensions, 0. );

        energy = 0;
        for ( unsigned int I = 0; I < spatialDimensions; I++ ){
            trE += E[ spatialDimensions * I + I ];
            for ( unsigned int J = 0; J < spatialDimensions; J++ ){
                trEsq += E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ];

                for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                    dtrEdChi[ spatialDimensions * J + k ] += dEdChi[ spatialDimensions * I + I ][ spatialDimensions * J + k ];

                    for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                        dtrEsqdChi[ spatialDimensions * k + K ] += 2 * dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * E[ spatialDimensions * I + J ];

                    }

                }

            }

        }

        energy += ( 1. / detChi ) * ( 0.5 * lambda * trE * trE + mu * trEsq );

        for ( unsigned int k = 0; k < spatialDimensions; k++ ){

            for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                dEnergydChi[ spatialDimensions * k + K ] += ( 1 / detChi ) * (   lambda * dtrEdChi[ spatialDimensions * k + K ] * trE + mu * dtrEsqdChi[ spatialDimensions * k + K ]
                                                                               - invChi[ spatialDimensions * K + k ] * ( 0.5 * lambda * trE * trE + mu * trEsq ) );

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
                             floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi ){
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
        floatVector invChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            dDetChidChi = vectorTools::computeDDetAdJ( chi, spatialDimensions, spatialDimensions );

            invChi = vectorTools::inverse( chi, spatialDimensions, spatialDimensions );

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
        floatType trEsq = 0;
        floatVector dtrEdChi( spatialDimensions * spatialDimensions, 0. );
        floatVector dtrEsqdChi( spatialDimensions * spatialDimensions, 0. );
        floatVector d2trEdChi2( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0. );
        vectorTools::eye( d2trEdChi2 );
        floatVector d2trEsqdChi2( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0. );
        energy = 0;
        dEnergydChi = floatVector( spatialDimensions * spatialDimensions, 0. );

        floatVector eye( spatialDimensions * spatialDimensions, 0 );
        vectorTools::eye( eye );

        energy = 0;
        for ( unsigned int I = 0; I < spatialDimensions; I++ ){
            trE += E[ spatialDimensions * I + I ];
            for ( unsigned int J = 0; J < spatialDimensions; J++ ){
                trEsq += E[ spatialDimensions * I + J ] * E[ spatialDimensions * J + I ];

                for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                    dtrEdChi[ spatialDimensions * J + k ] += dEdChi[ spatialDimensions * I + I ][ spatialDimensions * J + k ];

                    for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                        dtrEsqdChi[ spatialDimensions * k + K ] += 2 * dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * E[ spatialDimensions * I + J ];
                        d2trEsqdChi2[ spatialDimensions * spatialDimensions * spatialDimensions * I + spatialDimensions * spatialDimensions * J + spatialDimensions * k + K ]
                            += 2 * eye[ spatialDimensions * I + k ] * E[ spatialDimensions * J + K ];

                        for ( unsigned int l = 0; l < spatialDimensions; l++ ){

                            for ( unsigned int L = 0; L < spatialDimensions; L++ ){

                                d2trEsqdChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                                    += 2 * dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] *  dEdChi[ spatialDimensions * J + I ][ spatialDimensions * l + L ];

                            }

                        }

                    }

                }

            }

        }

        energy += ( 1. / detChi ) * ( 0.5 * lambda * trE * trE + mu * trEsq );

        floatVector d2EnergydChi2( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0 );

        for ( unsigned int k = 0; k < spatialDimensions; k++ ){

            for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                dEnergydChi[ spatialDimensions * k + K ] += ( 1 / detChi ) * (   lambda * dtrEdChi[ spatialDimensions * k + K ] * trE + mu * dtrEsqdChi[ spatialDimensions * k + K ]
                                                                               - invChi[ spatialDimensions * K + k ] * ( 0.5 * lambda * trE * trE + mu * trEsq ) );

                for ( unsigned int l = 0; l < spatialDimensions; l++ ){

                    for ( unsigned int L = 0; L < spatialDimensions; L++ ){

                        d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                            += ( 1 / detChi ) * ( - invChi[ spatialDimensions * L + l ] *  (   lambda * dtrEdChi[ spatialDimensions * k + K ] * trE + mu * dtrEsqdChi[ spatialDimensions * k + K ]
                                                                                             - invChi[ spatialDimensions * K + k ] * ( 0.5 * lambda * trE * trE + mu * trEsq ) )
                                                  + lambda * d2trEdChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ] * trE + lambda * dtrEdChi[ spatialDimensions * k + K ] * dtrEdChi[ spatialDimensions * l + L ] + mu * d2trEsqdChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                                                  + invChi[ spatialDimensions * K + l ] * invChi[ spatialDimensions * L + k ] * ( 0.5 * lambda * trE * trE + mu * trEsq )
                                                  - invChi[ spatialDimensions * K + k ] * ( lambda * dtrEdChi[ spatialDimensions * l + L ] * trE + mu * dtrEsqdChi[ spatialDimensions * l + L ] ) );

                    }

                }

            }

        }

        // Use the energy gradient to compute the Cauchy stress
        cauchyStress = floatVector( spatialDimensions * spatialDimensions, 0 );
        dCauchyStressdChi = floatMatrix( spatialDimensions * spatialDimensions, floatVector( spatialDimensions * spatialDimensions, 0 ) );
        for ( unsigned int i = 0; i < spatialDimensions; i++ ){
            for ( unsigned int j = 0; j < spatialDimensions; j++ ){
                for ( unsigned int I = 0; I < spatialDimensions; I++ ){
                    cauchyStress[ spatialDimensions * i + j ] += dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] / detChi;

                    for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                        dCauchyStressdChi[ spatialDimensions * i + j ][ spatialDimensions * I + k ] +=  dEnergydChi[ spatialDimensions * i + k ] * eye[ spatialDimensions * j + I ] / detChi;

                        for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                            dCauchyStressdChi[ spatialDimensions * i + j ][ spatialDimensions * k + K ] += - ( 1 / detChi ) * invChi[ spatialDimensions * K + k ] * dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ]
                                                                                                         + d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * I + spatialDimensions * k + K ] * chi[ spatialDimensions * j + I ] / detChi;

                        }

                    }

                }
            }
        }

        return NULL;

    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi,
                             floatVector &d2EnergydChi2, floatMatrix &d2CauchyStressdChi2 ){
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
         * \param &d2EnergydChi2: The second gradient of the energy w.r.t. the micro deformation
         * \param &d2CauchyStressdChi2: The second gradient of the Cauchy stress w.r.t. the micro deformation
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
        floatVector invChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            dDetChidChi = vectorTools::computeDDetAdJ( chi, spatialDimensions, spatialDimensions );

            invChi = vectorTools::inverse( chi, spatialDimensions, spatialDimensions );

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
        floatType trEsq = 0;
        floatVector dtrEdChi( spatialDimensions * spatialDimensions, 0. );
        floatVector dtrEsqdChi( spatialDimensions * spatialDimensions, 0. );
        floatVector d2trEdChi2( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0. );
        vectorTools::eye( d2trEdChi2 );
        floatVector d2trEsqdChi2( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0. );
        floatVector d3trEsqdChi3( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0. );
        energy = 0;
        dEnergydChi = floatVector( spatialDimensions * spatialDimensions, 0. );

        floatVector eye( spatialDimensions * spatialDimensions, 0 );
        vectorTools::eye( eye );


        return NULL;

    }

}

