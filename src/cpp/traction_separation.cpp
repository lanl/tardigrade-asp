/**
  ******************************************************************************
  * \file traction_separation.cpp
  ******************************************************************************
  * The source file for an implementation of various traction separation laws
  * for use in ASP. We will start with linear laws and, if required, progress
  * from there.
  ******************************************************************************
  */

#include<traction_separation.h>

namespace tractionSeparation{

    errorOut computeCurrentDistance( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &grad_chi,
                                     floatVector &d ){
        /*!
         * Compute the distance in the current configuration where
         * 
         * \f$d_i = dx_i - \xi_i^1 + \xi_i^2 \f$
         * 
         * \f$dx_i = F_{iI} dX_I \f$
         * 
         * \f$\xi_i^1 = \chi_{iI} \Xi_I^1 \f$
         * 
         * \f$\xi_i^2 = \left(\chi_{iI} + \chi_{iI,J} dX_J \right) \Xi_I^2 \f$
         * 
         * \f$dX_I = \Xi_I^1 + D_I - \Xi_I^2 \f$
         * 
         * \param &Xi_1: The micro-position vector for the local particle
         * \param &Xi_2: The micro-position vector for the non-local particle
         * \param &D: The initial separation between the particles
         * \param &F: The deformation gradient \f$\frac{dx_i}{dX_I}\f$
         * \param &chi: The micro-deformation
         * \param &grad_chi: the gradient of the micro deformation w.r.t. the reference spatial variable
         * \param &d: The current separation between the particles
         */

        floatVector dX = Xi_1 + D - Xi_2;

        floatVector dx( dX.size( ), 0 );
        floatVector chi_2 = chi;

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                dx[ i ] += F[ dX.size( ) * i + I ] * dX[ I ];

                for ( unsigned int J = 0; J < dX.size( ); J++ ){

                    chi_2[ dX.size( ) * i + I ] += grad_chi[ dX.size( ) * dX.size( ) * i + dX.size( ) * I + J ] * dX[ J ];

                }

            }

        }

        floatVector xi_1( dX.size( ), 0 );
        floatVector xi_2( dX.size( ), 0 );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                xi_1[ i ] += chi[ dX.size( ) * i + I ] * Xi_1[ I ];
                xi_2[ i ] += chi_2[ dX.size( ) * i + I ] * Xi_2[ I ];

            }

        }

        d = dx - xi_1 + xi_2;

        return NULL;

    }

    errorOut computeCurrentDistance( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &grad_chi,
                                     floatVector &d, floatMatrix &dddF, floatMatrix &dddChi, floatMatrix &dddGradChi ){
        /*!
         * Compute the distance in the current configuration where
         * 
         * \f$d_i = dx_i - \xi_i^1 + \xi_i^2 \f$
         * 
         * \f$dx_i = F_{iI} dX_I \f$
         * 
         * \f$\xi_i^1 = \chi_{iI} \Xi_I^1 \f$
         * 
         * \f$\xi_i^2 = \left(\chi_{iI} + \chi_{iI,J} dX_J \right) \Xi_I^2 \f$
         * 
         * \f$dX_I = \Xi_I^1 + D_I - \Xi_I^2 \f$
         * 
         * \param &Xi_1: The micro-position vector for the local particle
         * \param &Xi_2: The micro-position vector for the non-local particle
         * \param &D: The initial separation between the particles
         * \param &F: The deformation gradient \f$\frac{dx_i}{dX_I}\f$
         * \param &chi: The micro-deformation
         * \param &grad_chi: the gradient of the micro deformation w.r.t. the reference spatial variable
         * \param &d: The current separation between the particles
         * \param &dddF: The gradient of the separation w.r.t. the deformation gradient
         * \param &dddChi: The gradient of the separation w.r.t. the micro deformation
         * \param &dddGradChi: The gradient of the separation w.r.t. the gradient of the micro deformation
         */

        floatVector dX = Xi_1 + D - Xi_2;

        floatVector dx( dX.size( ), 0 );
        floatVector chi_2 = chi;

        floatMatrix dxdF( dX.size( ), floatVector( F.size( ), 0 ) );

        floatMatrix dchi_2dGradChi( chi.size( ), floatVector( grad_chi.size( ), 0 ) );
        floatMatrix dchi_2dX( chi.size( ), floatVector( dX.size( ), 0 ) );

        floatVector eye( dX.size( ) * dX.size( ), 0 );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                dx[ i ] += F[ dX.size( ) * i + I ] * dX[ I ];

                for ( unsigned int J = 0; J < dX.size( ); J++ ){

                    chi_2[ dX.size( ) * i + I ] += grad_chi[ dX.size( ) * dX.size( ) * i + dX.size( ) * I + J ] * dX[ J ];

                    dxdF[ i ][ dX.size( ) * I + J ] += eye[ dX.size( ) * i + I ] * dX[ J ];

                    for ( unsigned int k = 0; k < dX.size( ); k++ ){

                        for ( unsigned int K = 0; K < dX.size( ); K++ ){

                            dchi_2dGradChi[ dX.size( ) * i + I ][ dX.size( ) * dX.size( ) * k + dX.size( ) * K + J ] += eye[ dX.size( ) * i + k ] * eye[ dX.size( ) * I + K ] * dX[ J ];

                        } 

                    }

                }

            }

        }

        floatVector xi_1( dX.size( ), 0 );
        floatVector xi_2( dX.size( ), 0 );

        floatMatrix dxi_1dChi( dX.size( ), floatVector( chi.size( ), 0 ) );
        floatMatrix dxi_2dChi_2( dX.size( ), floatVector( chi.size( ), 0 ) );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                xi_1[ i ] += chi[ dX.size( ) * i + I ] * Xi_1[ I ];
                xi_2[ i ] += chi_2[ dX.size( ) * i + I ] * Xi_2[ I ];

                for ( unsigned int K = 0; K < dX.size( ); K++ ){

                    dxi_1dChi[ i ][ dX.size( ) * I + K ] += eye[ dX.size( ) * i + I ] * Xi_1[ K ];
                    dxi_2dChi_2[ i ][ dX.size( ) * I + K ] += eye[ dX.size( ) * i + I ] * Xi_2[ K ];

                }

            }

        }

        d = dx - xi_1 + xi_2;

        dddF = dxdF;

        dddChi = -dxi_1dChi + dxi_2dChi_2;

        dddGradChi = vectorTools::dot( dxi_2dChi_2, dchi_2dGradChi );

        return NULL;

    }

    errorOut computeCurrentDistance( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &grad_chi,
                                     floatVector &d, floatMatrix &dddF, floatMatrix &dddChi, floatMatrix &dddGradChi,
                                     floatMatrix &d2ddFdF, floatMatrix &d2ddFdChi, floatMatrix &d2ddFdGradChi,
                                     floatMatrix &d2ddChidChi, floatMatrix &d2ddChidGradChi,
                                     floatMatrix &d2ddGradChidGradChi );

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt );

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt,
                              floatMatrix &ddndd, floatMatrix &ddndn );

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt,
                              floatMatrix &ddndd, floatMatrix &ddndn,
                              floatMatrix &d2dndddd, floatMatrix &d2dndddn,
                              floatMatrix &d2dndndn,
                              floatMatrix &d2dtdddd, floatMatrix &d2dtdddn,
                              floatMatrix &d2dtdndn );

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy );


    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt );

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt,
                                          floatVector &d2energyddnddn, floatVector &d2energyddnddt,
                                          floatVector &d2energyddtddt );

}
