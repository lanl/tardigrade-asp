#include<constraint_equations.h>
#include<sstream>

namespace constraintEquations{

    errorOut tractionConstraint( const floatVector &cauchyStress, const floatVector &n, const floatVector &traction, const floatType &P, floatType &C ){
        /*!
         * Compute the traction constraint where
         * 
         * \f$ C = P \left( \sigma_{ij} n_i - t_j \right) \left( \sigma_{kj} n_k - t_j \right) \f$
         * 
         * \param &cauchyStress: The cauchy stress in row-major vector format
         * \param &n: The normal vector in the current configuration
         * \param &traction: The traction vector in the current configuration
         * \param &P: The penalty parameter or a lagrange multiplier depending on how it will be used
         * \param &C: The resulting error
         */

        floatVector error;
        error = -traction;

        for ( unsigned int i = 0; i < n.size( ); i++ ){

            for ( unsigned int j = 0; j < n.size( ); j++ ){

                error[ i ] += cauchyStress[ n.size( ) * j + i ] * n[ j ];

            }

        }

        C = P * vectorTools::dot( error, error );

        return NULL;

    }

    errorOut tractionConstraint( const floatVector &cauchyStress, const floatVector &n, const floatVector &traction, const floatType &P, floatType &C,
                                 floatVector &dCdCauchyStress, floatVector &dCdNormal, floatVector &dCdTraction, floatType &dCdP ){
        /*!
         * Compute the traction constraint where
         * 
         * \f$ C = P \left( \sigma_{ij} n_i - t_j \right) \left( \sigma_{kj} n_k - t_j \right) \f$
         * 
         * \param &cauchyStress: The cauchy stress in row-major vector format
         * \param &n: The normal vector in the current configuration
         * \param &traction: The traction vector in the current configuration
         * \param &P: The penalty parameter or a lagrange multiplier depending on how it will be used
         * \param &C: The resulting error
         * \param &dCdCauchyStress: The derivative of the constraint w.r.t. the Cauchy stress
         * \param &dCdNormal: The derivative of the constraint w.r.t. the normal
         * \param &dCdTraction: The derivative of the constraint w.r.t. the traction
         * \param &dCdP: The derivative of the constraint w.r.t. the penalty parameter
         */

        floatVector error;

        error = -traction;

        dCdCauchyStress = floatVector( cauchyStress.size( ), 0 );

        dCdNormal = floatVector( n.size( ), 0 );

        dCdTraction = floatVector( traction.size( ), 0 );

        for ( unsigned int i = 0; i < n.size( ); i++ ){

            dCdTraction[ i ] += 2 * P * traction[ i ];

            for ( unsigned int j = 0; j < n.size( ); j++ ){

                error[ i ] += cauchyStress[ n.size( ) * j + i ] * n[ j ];

                dCdCauchyStress[ n.size( ) * i + j ] += -2 * P * n[ i ] * traction[ j ];

                dCdNormal[ i ] += -2 * P * cauchyStress[ n.size( ) * i + j ] * traction[ j ];

                dCdTraction[ i ] += -2 * P * cauchyStress[ n.size( ) * j + i ] * n[ j ];

                for ( unsigned int k = 0; k < n.size( ); k++ ){

                    dCdCauchyStress[ n.size( ) * i + j ] += 2 * P * n[ i ] * cauchyStress[ n.size( ) * k + j ] * n[ k ];

                    dCdNormal[ i ] += 2 * P * cauchyStress[ n.size( ) * i + j ] * cauchyStress[ n.size( ) * k + j ] * n[ k ];

                }

            }

        }

        C = P * vectorTools::dot( error, error );

        dCdP = vectorTools::dot( error, error );

        return NULL;

    }

    errorOut tractionConstraint( const floatVector &cauchyStress, const floatVector &n, const floatVector &traction, const floatType &P, floatType &C,
                                 floatVector &dCdCauchyStress, floatVector &dCdNormal, floatVector &dCdTraction, floatType &dCdP,
                                 floatVector &d2CdCauchyStressdNormal, floatVector &d2CdCauchyStressdP,
                                 floatVector &d2CdNormaldP,            floatVector &d2CdTractiondP ){
        /*!
         * Compute the traction constraint where
         * 
         * \f$ C = P \left( \sigma_{ij} n_i - t_j \right) \left( \sigma_{kj} n_k - t_j \right) \f$
         * 
         * \param &cauchyStress: The cauchy stress in row-major vector format
         * \param &n: The normal vector in the current configuration
         * \param &traction: The traction vector in the current configuration
         * \param &P: The penalty parameter or a lagrange multiplier depending on how it will be used
         * \param &C: The resulting error
         * \param &dCdCauchyStress: The derivative of the constraint w.r.t. the Cauchy stress
         * \param &dCdNormal: The derivative of the constraint w.r.t. the normal
         * \param &dCdTraction: The derivative of the constraint w.r.t. the traction
         * \param &dCdP: The derivative of the constraint w.r.t. the penalty parameter
         * \param &dCdCauchyStressdNormal: The second derivative of the constraint w.r.t. the Cauchy stress and normal vector
         * \param &dCdCauchyStressdP: The second derivative of the constraint w.r.t. the Cauchy stress and the penalty parameter
         * \param &dCdNormaldP: The second derivative of the constraint w.r.t. the normal and the penalty parameter
         * \param &dCdTractiondP: The second derivative of the constraint w.r.t. the traction and the penalty parameter
         */

        floatVector error;

        error = -traction;

        d2CdCauchyStressdNormal = floatVector( cauchyStress.size( ) * n.size( ), 0 );

        d2CdCauchyStressdP = floatVector( cauchyStress.size( ), 0 );

        d2CdNormaldP = floatVector( n.size( ), 0 );

        d2CdTractiondP = floatVector( traction.size( ), 0 );

        floatVector eye( n.size( ) * n.size( ), 0 );

        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < n.size( ); i++ ){

            d2CdTractiondP[ i ] += 2 * traction[ i ];

            for ( unsigned int j = 0; j < n.size( ); j++ ){

                error[ i ] += cauchyStress[ n.size( ) * j + i ] * n[ j ];

                d2CdCauchyStressdP[ n.size( ) * i + j ] += -2 * n[ i ] * traction[ j ];

                d2CdNormaldP[ i ] += -2 * cauchyStress[ n.size( ) * i + j ] * traction[ j ];

                d2CdTractiondP[ i ] += -2 * cauchyStress[ n.size( ) * j + i ] * n[ j ];

                for ( unsigned int k = 0; k < n.size( ); k++ ){

                    d2CdCauchyStressdP[ n.size( ) * i + j ] += 2 * n[ i ] * cauchyStress[ n.size( ) * k + j ] * n[ k ];

                    d2CdNormaldP[ i ] += 2 * cauchyStress[ n.size( ) * i + j ] * cauchyStress[ n.size( ) * k + j ] * n[ k ];

                    d2CdCauchyStressdNormal[ n.size( ) * n.size( ) * i + n.size( ) * j + k ] += -2 * P * ( eye[ n.size( ) * i + k ] * traction[ j ] - n[ i ] * cauchyStress[ n.size( ) * k + j ] );

                    for ( unsigned int l = 0; l < n.size( ); l++ ){

                        d2CdCauchyStressdNormal[ n.size( ) * n.size( ) * i + n.size( ) * j + k ] += 2 * P * eye[ n.size( ) * i + k ] * cauchyStress[ n.size( ) * l + j ] * n[ l ];

                    }

                }

            }

        }

//        C = P ( \sigma_{ij} n_i - t_i )( \sigma_{kj} n_k - t_j )
//
//        dCdsigma_ab = 2 P n_i ( \sigma_{kj} n_k - t_j )
//        dCdsigma_abc = 2 P \delta_{ic} ( \sigma_{lj} n_l - t_j ) + 2 P n_i \sigma_{cj}
//                     = 2 P \delta_{ik} \sigma_{lj} n_l - 2 P \delta_{ik} t_j + 2 P n_i \sigma_{kj}

        C = P * vectorTools::dot( error, error );

        dCdP = vectorTools::dot( error, error );

        dCdCauchyStress = P * d2CdCauchyStressdP;

        dCdNormal = P * d2CdNormaldP;

        dCdTraction = P * d2CdTractiondP;

        return NULL;

    }

}
