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
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
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
         * \param &gradChi: the gradient of the micro deformation w.r.t. the reference spatial variable
         * \param &d: The current separation between the particles
         */

        floatVector dX = Xi_1 + D - Xi_2;

        floatVector dx( dX.size( ), 0 );
        floatVector chi_2 = chi;

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                dx[ i ] += F[ dX.size( ) * i + I ] * dX[ I ];

                for ( unsigned int J = 0; J < dX.size( ); J++ ){

                    chi_2[ dX.size( ) * i + I ] += gradChi[ dX.size( ) * dX.size( ) * i + dX.size( ) * I + J ] * dX[ J ];

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
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &d,
                                     floatMatrix &dddXi_1, floatMatrix &dddXi_2, floatMatrix &dddD,
                                     floatMatrix &dddF, floatMatrix &dddChi, floatMatrix &dddGradChi ){
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
         * \param &gradChi: the gradient of the micro deformation w.r.t. the reference spatial variable
         * \param &d: The current separation between the particles
         * \param &dddXi_1: The gradient of the separation w.r.t. the local reference micro-position vector
         * \param &dddXi_2: The gradient of the separation w.r.t. the non-local reference micro-position vector
         * \param &dddD: The gradient of the separation w.r.t. the local reference separation distance
         * \param &dddF: The gradient of the separation w.r.t. the deformation gradient
         * \param &dddChi: The gradient of the separation w.r.t. the micro deformation
         * \param &dddGradChi: The gradient of the separation w.r.t. the gradient of the micro deformation
         */

        floatVector dX = Xi_1 + D - Xi_2;

        floatVector dx( dX.size( ), 0 );
        floatVector chi_2 = chi;

        floatMatrix dxdF( dX.size( ), floatVector( F.size( ), 0 ) );
        floatMatrix dxdX( dX.size( ), floatVector( dX.size( ), 0 ) );

        floatMatrix dchi_2dGradChi( chi.size( ), floatVector( gradChi.size( ), 0 ) );
        floatMatrix dchi_2dX( chi.size( ), floatVector( dX.size( ), 0 ) );

        floatVector eye( dX.size( ) * dX.size( ), 0 );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                dx[ i ] += F[ dX.size( ) * i + I ] * dX[ I ];

                dxdX[ i ][ I ] += F[ dX.size( ) * i + I ];

                for ( unsigned int J = 0; J < dX.size( ); J++ ){

                    chi_2[ dX.size( ) * i + I ] += gradChi[ dX.size( ) * dX.size( ) * i + dX.size( ) * I + J ] * dX[ J ];

                    dxdF[ i ][ dX.size( ) * I + J ] += eye[ dX.size( ) * i + I ] * dX[ J ];

                    dchi_2dX[ dX.size( ) * i + I ][ J ] += gradChi[ dX.size( ) * dX.size( ) * i + dX.size( ) * I + J ];

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

        floatMatrix dxi_1dXi_1( dX.size( ), floatVector( Xi_1.size( ), 0 ) );
        floatMatrix dxi_2dXi_2( dX.size( ), floatVector( Xi_2.size( ), 0 ) );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                xi_1[ i ] += chi[ dX.size( ) * i + I ] * Xi_1[ I ];
                xi_2[ i ] += chi_2[ dX.size( ) * i + I ] * Xi_2[ I ];

                dxi_1dXi_1[ i ][ I ] += chi[ dX.size( ) * i + I ];
                dxi_2dXi_2[ i ][ I ] += chi_2[ dX.size( ) * i + I ];

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

        dddXi_1 = dxdX - dxi_1dXi_1 + vectorTools::dot( dxi_2dChi_2, dchi_2dX );

        dddXi_2 =  -dxdX + dxi_2dXi_2 - vectorTools::dot( dxi_2dChi_2, dchi_2dX );

        dddD    =  dxdX + vectorTools::dot( dxi_2dChi_2, dchi_2dX );

        return NULL;

    }

    errorOut computeCurrentDistance( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &d,
                                     floatMatrix &dddXi_1, floatMatrix &dddXi_2, floatMatrix &dddD,
                                     floatMatrix &dddF, floatMatrix &dddChi, floatMatrix &dddGradChi,
                                     floatMatrix &d2ddFdXi_1,       floatMatrix &d2ddFdXi_2,       floatMatrix &d2ddFdD,
                                     floatMatrix &d2ddChidXi_1,     floatMatrix &d2ddChidXi_2,     floatMatrix &d2ddChidD,
                                     floatMatrix &d2ddGradChidXi_1, floatMatrix &d2ddGradChidXi_2, floatMatrix &d2ddGradChidD ){
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
         * \param &gradChi: the gradient of the micro deformation w.r.t. the reference spatial variable
         * \param &d: The current separation between the particles
         * \param &dddF: The gradient of the separation w.r.t. the deformation gradient
         * \param &dddChi: The gradient of the separation w.r.t. the micro deformation
         * \param &dddGradChi: The gradient of the separation w.r.t. the gradient of the micro deformation
         * \param &d2ddFdXi_1: The second derivative of d w.r.t. F and \f$\Xi_1\f$
         * \param &d2ddFdXi_2: The second derivative of d w.r.t. F and \f$\Xi_2\f$
         * \param &d2ddFdD: The second derivative of d w.r.t. F and D
         * \param &d2ddChidXi_1: The second derivative of d w.r.t. \f$\chi\f$ and \f$\Xi_1\f$
         * \param &d2ddChidXi_2: The second derivative of d w.r.t. \f$\chi\f$ and \f$\Xi_2\f$
         * \param &d2ddChidD: The second derivative of d w.r.t. \f$\chi\f$ and D
         * \param &d2ddGradChidXi_1: The second derivative of d w.r.t. \f$\nabla_X \chi\f$ and \f$\Xi_1\f$
         * \param &d2ddGradChidXi_2: The second derivative of d w.r.t. \f$\nabla_X \chi\f$ and \f$\Xi_2\f$
         * \param &d2ddGradChidD: The second derivative of d w.r.t. \f$\nabla_X \chi\f$ and D
         */

        floatVector dX = Xi_1 + D - Xi_2;

        floatVector dx( dX.size( ), 0 );
        floatVector chi_2 = chi;

        floatMatrix dxdF( dX.size( ), floatVector( F.size( ), 0 ) );
        floatMatrix dxdX( dX.size( ), floatVector( dX.size( ), 0 ) );

        floatMatrix d2xdFdX( dX.size( ), floatVector( F.size( ) * dX.size( ), 0 ) );

        floatMatrix dchi_2dGradChi( chi.size( ), floatVector( gradChi.size( ), 0 ) );
        floatMatrix dchi_2dX( chi.size( ), floatVector( dX.size( ), 0 ) );

        floatMatrix d2chi_2dGradChidX( chi.size( ), floatVector( gradChi.size( ) * dX.size( ), 0 ) );

        floatVector eye( dX.size( ) * dX.size( ), 0 );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                dx[ i ] += F[ dX.size( ) * i + I ] * dX[ I ];

                dxdX[ i ][ I ] += F[ dX.size( ) * i + I ];

                for ( unsigned int J = 0; J < dX.size( ); J++ ){

                    chi_2[ dX.size( ) * i + I ] += gradChi[ dX.size( ) * dX.size( ) * i + dX.size( ) * I + J ] * dX[ J ];

                    dxdF[ i ][ dX.size( ) * I + J ] += eye[ dX.size( ) * i + I ] * dX[ J ];

                    dchi_2dX[ dX.size( ) * i + I ][ J ] += gradChi[ dX.size( ) * dX.size( ) * i + dX.size( ) * I + J ];

                    for ( unsigned int k = 0; k < dX.size( ); k++ ){

                        d2xdFdX[ i ][ dX.size( ) * dX.size( ) * I + dX.size( ) * J + k ] += eye[ dX.size( ) * i + I ] * eye[ dX.size( ) * J + k ];

                        for ( unsigned int K = 0; K < dX.size( ); K++ ){

                            dchi_2dGradChi[ dX.size( ) * i + I ][ dX.size( ) * dX.size( ) * k + dX.size( ) * K + J ] += eye[ dX.size( ) * i + k ] * eye[ dX.size( ) * I + K ] * dX[ J ];

                            for ( unsigned int L = 0; L < dX.size( ); L++ ){

                                d2chi_2dGradChidX[ dX.size( ) * i + I ][ dX.size( ) * dX.size( ) * dX.size( ) * k + dX.size( ) * dX.size( ) * K + dX.size( ) * J + L ]
                                    += eye[ dX.size( ) * i + k ] * eye[ dX.size( ) * I + K ] * eye[ dX.size( ) * J + L ];

                            }

                        } 

                    }

                }

            }

        }

        floatVector xi_1( dX.size( ), 0 );
        floatVector xi_2( dX.size( ), 0 );

        floatMatrix dxi_1dChi( dX.size( ), floatVector( chi.size( ), 0 ) );
        floatMatrix dxi_2dChi_2( dX.size( ), floatVector( chi.size( ), 0 ) );

        floatMatrix d2xi_1dChidXi_1( dX.size( ), floatVector( chi.size( ) * Xi_1.size( ), 0 ) );
        floatMatrix d2xi_2dChi_2dXi_2( dX.size( ), floatVector( chi.size( ) * Xi_1.size( ), 0 ) );

        floatMatrix dxi_1dXi_1( dX.size( ), floatVector( Xi_1.size( ), 0 ) );
        floatMatrix dxi_2dXi_2( dX.size( ), floatVector( Xi_2.size( ), 0 ) );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){

            for ( unsigned int I = 0; I < dX.size( ); I++ ){

                xi_1[ i ] += chi[ dX.size( ) * i + I ] * Xi_1[ I ];
                xi_2[ i ] += chi_2[ dX.size( ) * i + I ] * Xi_2[ I ];

                dxi_1dXi_1[ i ][ I ] += chi[ dX.size( ) * i + I ];
                dxi_2dXi_2[ i ][ I ] += chi_2[ dX.size( ) * i + I ];

                for ( unsigned int K = 0; K < dX.size( ); K++ ){

                    dxi_1dChi[ i ][ dX.size( ) * I + K ] += eye[ dX.size( ) * i + I ] * Xi_1[ K ];
                    dxi_2dChi_2[ i ][ dX.size( ) * I + K ] += eye[ dX.size( ) * i + I ] * Xi_2[ K ];

                    for ( unsigned int L = 0; L < dX.size( ); L++ ){

                        d2xi_1dChidXi_1[ i ][ dX.size( ) * dX.size( ) * I + dX.size( ) * K + L ] += eye[ dX.size( ) * i + I ] * eye[ dX.size( ) * K + L ];
                        d2xi_2dChi_2dXi_2[ i ][ dX.size( ) * dX.size( ) * I + dX.size( ) * K + L ] += eye[ dX.size( ) * i + I ] * eye[ dX.size( ) * K + L ];

                    }

                }

            }

        }

        d = dx - xi_1 + xi_2;

        dddF = dxdF;

        dddChi = -dxi_1dChi + dxi_2dChi_2;

        dddGradChi = vectorTools::dot( dxi_2dChi_2, dchi_2dGradChi );

        dddXi_1 = dxdX - dxi_1dXi_1 + vectorTools::dot( dxi_2dChi_2, dchi_2dX );

        dddXi_2 =  -dxdX + dxi_2dXi_2 - vectorTools::dot( dxi_2dChi_2, dchi_2dX );

        dddD    =  dxdX + vectorTools::dot( dxi_2dChi_2, dchi_2dX );

        d2ddFdXi_1       = d2xdFdX;
        d2ddChidXi_1     = -d2xi_1dChidXi_1;
        d2ddGradChidXi_1 = vectorTools::dot( dxi_2dChi_2, d2chi_2dGradChidX );

        d2ddFdXi_2       = -d2xdFdX;
        d2ddChidXi_2     = d2xi_2dChi_2dXi_2;
        d2ddGradChidXi_2 = -vectorTools::dot( dxi_2dChi_2, d2chi_2dGradChidX );

        d2ddFdD          = d2xdFdX;
        d2ddChidD        = floatMatrix( d.size( ), floatVector( chi.size( ) * D.size( ), 0 ) );
        d2ddGradChidD    = vectorTools::dot( dxi_2dChi_2, d2chi_2dGradChidX );

        for ( unsigned int i = 0; i < dX.size( ); i++ ){
            for ( unsigned int j = 0; j < dX.size( ); j++ ){
                for ( unsigned int J = 0; J < dX.size( ); J++ ){
                    for ( unsigned int K = 0; K < dX.size( ); K++ ){
                        for ( unsigned int L = 0; L < dX.size( ); L++ ){
                            for ( unsigned int a = 0; a < dX.size( ); a++ ){
                                for ( unsigned int A = 0; A < dX.size( ); A++ ){
                                     d2ddGradChidXi_2[ i ][ dX.size( ) * dX.size( ) * dX.size( ) * j + dX.size( ) * dX.size( ) * J + dX.size( ) * K + L ]
                                         += d2xi_2dChi_2dXi_2[ i ][ dX.size( ) * dX.size( ) * a + dX.size( ) * A + L ] * dchi_2dGradChi[ dX.size( ) * a + A ][ dX.size( ) * dX.size( ) * j + dX.size( ) * J + K ];
                                }

                            }

                        }

                    }

                }

            }

        }

        return NULL;

    }

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt ){
        /*!
         * Decompose a vector into a normal part and a tangential part via
         * 
         * \f$ d^n_i = d_a n_a n_i \f$
         * 
         * \f$ d^t_i = d_i - d^n_i \f$
         * 
         * \param &d: The vector to decompose
         * \param &n: The normal vector ( assumed to be unit length ).
         * \param &dn: The part of the vector in the normal direction
         * \param &dt: The tangential part of the vector
         */

        if ( !vectorTools::fuzzyEquals( vectorTools::l2norm( n ), 1. ) ){
            return new errorNode( __func__, "The normal vector isn't a unit vector!" );
        }

        dn = vectorTools::dot( d, n ) * n;

        dt = d - dn;

        return NULL;

    }

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt,
                              floatMatrix &ddndd, floatMatrix &ddndn,
                              floatMatrix &ddtdd, floatMatrix &ddtdn ){
        /*!
         * Decompose a vector into a normal part and a tangential part via
         * 
         * \f$ d^n_i = d_a n_a n_i \f$
         * 
         * \f$ d^t_i = d_i - d^n_i \f$
         * 
         * \param &d: The vector to decompose
         * \param &n: The normal vector ( assumed to be unit length ).
         * \param &dn: The part of the vector in the normal direction
         * \param &dt: The tangential part of the vector
         * \param &ddndd: The derivative of the normal part of the vector w.r.t. the original vector
         * \param &ddndn: The derivative of the normal part of the vector w.r.t. the normal vector
         * \param &ddtdd: The derivative of the tangential part of the vector w.r.t. the original vector
         * \param &ddtdn: The derivative of the tangential part of the vector w.r.t. the normal vector
         */

        if ( !vectorTools::fuzzyEquals( vectorTools::l2norm( n ), 1. ) ){
            return new errorNode( __func__, "The normal vector isn't a unit vector!" );
        }

        dn = vectorTools::dot( d, n ) * n;

        ddndd = vectorTools::dyadic( n, n );

        ddndn = vectorTools::dyadic( n, d ) + vectorTools::dot( d, n ) * vectorTools::eye<floatType>( d.size( ) );

        dt = d - dn;

        ddtdd = vectorTools::eye<floatType>( d.size( ) ) - ddndd;

        ddtdn = -ddndn;

        return NULL;

    }

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt,
                              floatMatrix &ddndd, floatMatrix &ddndn,
                              floatMatrix &ddtdd, floatMatrix &ddtdn,
                              floatMatrix &d2dndddd, floatMatrix &d2dndddn,
                              floatMatrix &d2dndndn,
                              floatMatrix &d2dtdddd, floatMatrix &d2dtdddn,
                              floatMatrix &d2dtdndn ){
        /*!
         * Decompose a vector into a normal part and a tangential part via
         * 
         * \f$ d^n_i = d_a n_a n_i \f$
         * 
         * \f$ d^t_i = d_i - d^n_i \f$
         * 
         * \param &d: The vector to decompose
         * \param &n: The normal vector ( assumed to be unit length ).
         * \param &dn: The part of the vector in the normal direction
         * \param &dt: The tangential part of the vector
         * \param &ddndd: The derivative of the normal part of the vector w.r.t. the original vector
         * \param &ddndn: The derivative of the normal part of the vector w.r.t. the normal vector
         * \param &ddtdd: The derivative of the tangential part of the vector w.r.t. the original vector
         * \param &ddtdn: The derivative of the tangential part of the vector w.r.t. the normal vector
         * \param &d2dndddd: The second derivative of the normal part of the vector w.r.t. the original vector
         * \param &d2dndddn: The second derivative of the normal part of the vector w.r.t. the original vector and the normal vector
         * \param &d2dndndn: The second derivative of the normal part of the vector w.r.t. the normal vector
         * \param &d2dtdddd: The second derivative of the tangential part of the vector w.r.t. the original vector
         * \param &d2dtdddn: The second derivative of the tangential part of the vector w.r.t. the original vector and the normal vector
         * \param &d2dtdndn: The second derivative of the tangential part of the vector w.r.t. the normal vector
         */

        if ( !vectorTools::fuzzyEquals( vectorTools::l2norm( n ), 1. ) ){
            return new errorNode( __func__, "The normal vector isn't a unit vector!" );
        }

        dn = vectorTools::dot( d, n ) * n;

        ddndd = vectorTools::dyadic( n, n );

        ddndn = vectorTools::dyadic( n, d ) + vectorTools::dot( d, n ) * vectorTools::eye<floatType>( d.size( ) );

        dt = d - dn;

        ddtdd = vectorTools::eye<floatType>( d.size( ) ) - ddndd;

        ddtdn = -ddndn;

        d2dndddd = floatMatrix( d.size( ), floatVector( d.size( ) * d.size( ), 0 ) );
        d2dndddn = floatMatrix( d.size( ), floatVector( d.size( ) * d.size( ), 0 ) );
        d2dndndn = floatMatrix( d.size( ), floatVector( d.size( ) * d.size( ), 0 ) );

        d2dtdddd = floatMatrix( d.size( ), floatVector( d.size( ) * d.size( ), 0 ) );
        d2dtdddn = floatMatrix( d.size( ), floatVector( d.size( ) * d.size( ), 0 ) );
        d2dtdndn = floatMatrix( d.size( ), floatVector( d.size( ) * d.size( ), 0 ) );

        floatVector eye( d.size( ) * d.size( ) );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < d.size( ); i++ ){

            for ( unsigned int j = 0; j < d.size( ); j++ ){

                for ( unsigned int k = 0; k < d.size( ); k++ ){

                    d2dndddn[ i ][ d.size( ) * j + k ] += eye[ d.size( ) * i + k ] * n[ j ] + n[ i ] * eye[ d.size( ) * j + k ];
                    d2dndndn[ i ][ d.size( ) * j + k ] += eye[ d.size( ) * i + k ] * d[ j ] + d[ k ] * eye[ d.size( ) * i + j ];

                }

            }

        }

        d2dtdddn = -d2dndddn;
        d2dtdndn = -d2dndndn;

        return NULL;

    }

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy ){
        /*!
         * Compute the linear traction-separation energy
         * 
         * \f$ e^t = \frac{1}{2} \left[ E^n d^n_i d^n_i + E^t d^t_i d^t_i\right]\f$
         * 
         * \param &normalDeformationMeasure: The normal deformation measure \f$ d^n_i \f$
         * \param &tangentialDeformationMeasure: The tangential deformation measure \f$ d^t_i \f$
         * \param &parameters: The material parameters \f$ E^n \f$ and \f$ E^t \f$.
         * \param &energy: The returned energy value \f$ e^t \f$
         */

        if ( parameters.size( ) != 2 ){

            return new errorNode( __func__, "Two parameters are required for the traction separation law. " + std::to_string( parameters.size( ) ) + " are provided." );

        }

        floatType En = parameters[ 0 ];

        floatType Et = parameters[ 1 ];

        energy = 0.5 * ( En * vectorTools::dot( normalDeformationMeasure, normalDeformationMeasure ) + Et * vectorTools::dot( tangentialDeformationMeasure, tangentialDeformationMeasure ) );

        return NULL;

    }

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt ){
        /*!
         * Compute the linear traction-separation energy
         * 
         * \f$ e^t = \frac{1}{2} \left[ E^n d^n_i d^n_i + E^t d^t_i d^t_i\right]\f$
         * 
         * \param &normalDeformationMeasure: The normal deformation measure \f$ d^n_i \f$
         * \param &tangentialDeformationMeasure: The tangential deformation measure \f$ d^t_i \f$
         * \param &parameters: The material parameters \f$ E^n \f$ and \f$ E^t \f$.
         * \param &energy: The returned energy value \f$ e^t \f$
         * \param &denergyddn: The derivative of the energy w.r.t. the normal deformation measure
         * \param &denergyddt: The derivative of the energy w.r.t. the tangential deformation measure
         */

        if ( parameters.size( ) != 2 ){

            return new errorNode( __func__, "Two parameters are required for the traction separation law. " + std::to_string( parameters.size( ) ) + " are provided." );

        }

        floatType En = parameters[ 0 ];

        floatType Et = parameters[ 1 ];

        energy = 0.5 * ( En * vectorTools::dot( normalDeformationMeasure, normalDeformationMeasure ) + Et * vectorTools::dot( tangentialDeformationMeasure, tangentialDeformationMeasure ) );

        denergyddn = En * normalDeformationMeasure;

        denergyddt = Et * tangentialDeformationMeasure;

        return NULL;

    }

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt,
                                          floatVector &d2energyddnddn, floatVector &d2energyddnddt,
                                          floatVector &d2energyddtddt ){
        /*!
         * Compute the linear traction-separation energy
         * 
         * \f$ e^t = \frac{1}{2} \left[ E^n d^n_i d^n_i + E^t d^t_i d^t_i\right]\f$
         * 
         * \param &normalDeformationMeasure: The normal deformation measure \f$ d^n_i \f$
         * \param &tangentialDeformationMeasure: The tangential deformation measure \f$ d^t_i \f$
         * \param &parameters: The material parameters \f$ E^n \f$ and \f$ E^t \f$.
         * \param &energy: The returned energy value \f$ e^t \f$
         * \param &denergyddn: The derivative of the energy w.r.t. the normal deformation measure
         * \param &denergyddt: The derivative of the energy w.r.t. the tangential deformation measure
         * \param &d2energyddnddn: The second derivative of the energy w.r.t. the normal deformation measure
         * \param &d2energyddnddt: The second derivative of the energy w.r.t. the normal and tangential deformation measures
         * \param &d2energyddtddt: The second derivative of the energy w.r.t. the tangential deformation measure
         */

        if ( parameters.size( ) != 2 ){

            return new errorNode( __func__, "Two parameters are required for the traction separation law. " + std::to_string( parameters.size( ) ) + " are provided." );

        }

        floatType En = parameters[ 0 ];

        floatType Et = parameters[ 1 ];

        energy = 0.5 * ( En * vectorTools::dot( normalDeformationMeasure, normalDeformationMeasure ) + Et * vectorTools::dot( tangentialDeformationMeasure, tangentialDeformationMeasure ) );

        denergyddn = En * normalDeformationMeasure;

        denergyddt = Et * tangentialDeformationMeasure;

        floatVector eye( normalDeformationMeasure.size( ) * normalDeformationMeasure.size( ) );

        vectorTools::eye( eye );

        d2energyddnddn = En * eye;

        d2energyddnddt = floatVector( eye.size( ), 0 );

        d2energyddtddt = Et * eye;

        return NULL;

    }

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt, floatVector &denergydParameters ){
        /*!
         * Compute the linear traction-separation energy
         * 
         * \f$ e^t = \frac{1}{2} \left[ E^n d^n_i d^n_i + E^t d^t_i d^t_i\right]\f$
         * 
         * \param &normalDeformationMeasure: The normal deformation measure \f$ d^n_i \f$
         * \param &tangentialDeformationMeasure: The tangential deformation measure \f$ d^t_i \f$
         * \param &parameters: The material parameters \f$ E^n \f$ and \f$ E^t \f$.
         * \param &energy: The returned energy value \f$ e^t \f$
         * \param &denergyddn: The derivative of the energy w.r.t. the normal deformation measure
         * \param &denergyddt: The derivative of the energy w.r.t. the tangential deformation measure
         * \param &denergydParameters: The derivative of the energy w.r.t. the parameters
         */

        if ( parameters.size( ) != 2 ){

            return new errorNode( __func__, "Two parameters are required for the traction separation law. " + std::to_string( parameters.size( ) ) + " are provided." );

        }

        floatType En = parameters[ 0 ];

        floatType Et = parameters[ 1 ];

        energy = 0.5 * ( En * vectorTools::dot( normalDeformationMeasure, normalDeformationMeasure ) + Et * vectorTools::dot( tangentialDeformationMeasure, tangentialDeformationMeasure ) );

        denergyddn = En * normalDeformationMeasure;

        denergyddt = Et * tangentialDeformationMeasure;

        denergydParameters = { 0.5 * vectorTools::dot( normalDeformationMeasure, normalDeformationMeasure ),
                               0.5 * vectorTools::dot( tangentialDeformationMeasure, tangentialDeformationMeasure ) };

        return NULL;
    }

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt, floatVector &denergydParameters,
                                          floatVector &d2energyddnddn, floatVector &d2energyddnddt, floatVector &d2energyddndParameters,
                                          floatVector &d2energyddtddt, floatVector &d2energyddtdParameters,
                                          floatVector &d2energydParametersdParameters ){
        /*!
         * Compute the linear traction-separation energy
         * 
         * \f$ e^t = \frac{1}{2} \left[ E^n d^n_i d^n_i + E^t d^t_i d^t_i\right]\f$
         * 
         * \param &normalDeformationMeasure: The normal deformation measure \f$ d^n_i \f$
         * \param &tangentialDeformationMeasure: The tangential deformation measure \f$ d^t_i \f$
         * \param &parameters: The material parameters \f$ E^n \f$ and \f$ E^t \f$.
         * \param &energy: The returned energy value \f$ e^t \f$
         * \param &denergyddn: The derivative of the energy w.r.t. the normal deformation measure
         * \param &denergyddt: The derivative of the energy w.r.t. the tangential deformation measure
         * \param &denergydParameters: The derivative of the energy w.r.t. the parameters
         * \param &d2energyddnddn: The second derivative of the energy w.r.t. the normal deformation measure
         * \param &d2energyddnddt: The second derivative of the energy w.r.t. the normal and tangential deformation measures
         * \param &d2energyddnddt: The second derivative of the energy w.r.t. the normal deformation measure and the parameters
         * \param &d2energyddtdParameters: The second derivative of the energy w.r.t. the tangential deformation measure and the parameters
         * \param &d2energyddtdParameters: The second derivative of the energy w.r.t. the parameters
         */

        if ( parameters.size( ) != 2 ){

            return new errorNode( __func__, "Two parameters are required for the traction separation law. " + std::to_string( parameters.size( ) ) + " are provided." );

        }

        floatType En = parameters[ 0 ];

        floatType Et = parameters[ 1 ];

        energy = 0.5 * ( En * vectorTools::dot( normalDeformationMeasure, normalDeformationMeasure ) + Et * vectorTools::dot( tangentialDeformationMeasure, tangentialDeformationMeasure ) );

        denergyddn = En * normalDeformationMeasure;

        denergyddt = Et * tangentialDeformationMeasure;

        denergydParameters = { 0.5 * vectorTools::dot( normalDeformationMeasure, normalDeformationMeasure ),
                               0.5 * vectorTools::dot( tangentialDeformationMeasure, tangentialDeformationMeasure ) };

        floatVector eye( normalDeformationMeasure.size( ) * normalDeformationMeasure.size( ) );

        vectorTools::eye( eye );

        d2energyddnddn = En * eye;

        d2energyddnddt = floatVector( eye.size( ), 0 );

        d2energyddndParameters = floatVector( normalDeformationMeasure.size( ) * parameters.size( ), 0 );

        for ( unsigned int i = 0; i < normalDeformationMeasure.size( ); i++ ){

            d2energyddndParameters[ parameters.size( ) * i ] += normalDeformationMeasure[ i ];

        }

        d2energyddtddt = Et * eye;

        d2energyddtdParameters = floatVector( tangentialDeformationMeasure.size( ) * parameters.size( ), 0 );

        for ( unsigned int i = 0; i < tangentialDeformationMeasure.size( ); i++ ){

            d2energyddtdParameters[ parameters.size( ) * i + 1 ] += tangentialDeformationMeasure[ i ];

        }

        d2energydParametersdParameters = floatVector( parameters.size( ) * parameters.size( ), 0 );

        return NULL;

    }

    errorOut computeNansonsRelation( const floatVector &deformationGradient, const floatVector &dAN, floatVector &dan ){
        /*!
         * Compute Nanson's relation
         * 
         * \f$ da n_i = J dA N_I F_{Ii}^{-1} \f$
         * 
         * where \f$ J \f$ is the determinant of the deformation gradient \f$ F \f$
         * 
         * \param &deformationGradient: The deformation gradient \f$ F \f$
         * \param &dAN: The product of the reference area \f$ dA \f$ and the reference normal \f$ N \f$
         * \param &dan: The mapping of \f$ dA N_I \f$ to the current configuration.
         */

        if ( deformationGradient.size( ) != ( dAN.size( ) * dAN.size( ) ) ){

            return new errorNode( __func__, "The deformation gradient must have " + std::to_string( dAN.size( ) * dAN.size( ) ) + " terms and has " + std::to_string( deformationGradient.size( ) ) );

        }

        floatVector Finv = vectorTools::inverse( deformationGradient, dAN.size( ), dAN.size( ) );

        floatType J = vectorTools::determinant( deformationGradient, dAN.size( ), dAN.size( ) );

        dan = floatVector( dAN.size( ), 0 );

        for ( unsigned int i = 0; i < dAN.size( ); i++ ){

            for ( unsigned int I = 0; I < dAN.size( ); I++ ){

                dan[ i ] += J * dAN[ I ] * Finv[ dAN.size( ) * I + i ];

            }

        }

        return NULL;

    }

    errorOut computeNansonsRelation( const floatVector &deformationGradient, const floatVector &dAN, floatVector &dan,
                                     floatMatrix &ddandF, floatMatrix &ddanddAN ){
        /*!
         * Compute Nanson's relation
         * 
         * \f$ da n_i = J dA N_I F_{Ii}^{-1} \f$
         * 
         * where \f$ J \f$ is the determinant of the deformation gradient \f$ F \f$
         * 
         * \param &deformationGradient: The deformation gradient \f$ F \f$
         * \param &dAN: The product of the reference area \f$ dA \f$ and the reference normal \f$ N \f$
         * \param &dan: The mapping of \f$ dA N_I \f$ to the current configuration.
         * \param &ddandF: The gradient of the current surface area w.r.t. the deformation gradient
         * \param &ddanddAN: The gradient of the current surface area w.r.t. the reference surface area
         */

        if ( deformationGradient.size( ) != ( dAN.size( ) * dAN.size( ) ) ){

            return new errorNode( __func__, "The deformation gradient must have " + std::to_string( dAN.size( ) * dAN.size( ) ) + " terms and has " + std::to_string( deformationGradient.size( ) ) );

        }

        floatVector Finv = vectorTools::inverse( deformationGradient, dAN.size( ), dAN.size( ) );

        floatType J = vectorTools::determinant( deformationGradient, dAN.size( ), dAN.size( ) );

        dan = floatVector( dAN.size( ), 0 );

        ddandF = floatMatrix( dan.size( ), floatVector( deformationGradient.size( ), 0 ) );

        ddanddAN = floatMatrix( dan.size( ), floatVector( dAN.size( ), 0 ) );

        for ( unsigned int i = 0; i < dAN.size( ); i++ ){

            for ( unsigned int I = 0; I < dAN.size( ); I++ ){

                dan[ i ] += J * dAN[ I ] * Finv[ dAN.size( ) * I + i ];

                ddanddAN[ i ][ I ] += J * Finv[ dAN.size( ) * I + i ];

                for ( unsigned int a = 0; a < dAN.size( ); a++ ){

                    for ( unsigned int A = 0; A < dAN.size( ); A++ ){

                        ddandF[ i ][ dAN.size( ) * a + A ] += J * dAN[ I ] * ( Finv[ dAN.size( ) * I + i ] * Finv[ dAN.size( ) * A + a ] - Finv[ dAN.size( ) * A + i ] * Finv[ dAN.size( ) * I + a ] );

                    }

                }

            }

        }

        return NULL;

    }

    errorOut computeNansonsRelation( const floatVector &deformationGradient, const floatVector &dAN, floatVector &dan,
                                     floatMatrix &ddandF, floatMatrix &ddanddAN,
                                     floatMatrix &d2dandFdF, floatMatrix &d2dandFddAN ){
        /*!
         * Compute Nanson's relation
         * 
         * \f$ da n_i = J dA N_I F_{Ii}^{-1} \f$
         * 
         * where \f$ J \f$ is the determinant of the deformation gradient \f$ F \f$
         * 
         * \param &deformationGradient: The deformation gradient \f$ F \f$
         * \param &dAN: The product of the reference area \f$ dA \f$ and the reference normal \f$ N \f$
         * \param &dan: The mapping of \f$ dA N_I \f$ to the current configuration.
         * \param &ddandF: The gradient of the current surface area w.r.t. the deformation gradient
         * \param &ddanddAN: The gradient of the current surface area w.r.t. the reference surface area
         * \param &d2dandFdF: The second derivative of the current surface area w.r.t. the deformation gradient
         * \param &d2dandFddAN: The second derivative of the current surface area w.r.t. the deformation gradient and reference surface area
         */

        if ( deformationGradient.size( ) != ( dAN.size( ) * dAN.size( ) ) ){

            return new errorNode( __func__, "The deformation gradient must have " + std::to_string( dAN.size( ) * dAN.size( ) ) + " terms and has " + std::to_string( deformationGradient.size( ) ) );

        }

        floatVector Finv = vectorTools::inverse( deformationGradient, dAN.size( ), dAN.size( ) );

        floatType J = vectorTools::determinant( deformationGradient, dAN.size( ), dAN.size( ) );

        dan = floatVector( dAN.size( ), 0 );

        ddandF = floatMatrix( dan.size( ), floatVector( deformationGradient.size( ), 0 ) );

        ddanddAN = floatMatrix( dan.size( ), floatVector( dAN.size( ), 0 ) );

        d2dandFdF = floatMatrix( dan.size( ), floatVector( deformationGradient.size( ) * deformationGradient.size( ), 0 ) );

        d2dandFddAN = floatMatrix( dan.size( ), floatVector( deformationGradient.size( ) * dAN.size( ), 0 ) );

        for ( unsigned int i = 0; i < dAN.size( ); i++ ){

            for ( unsigned int I = 0; I < dAN.size( ); I++ ){

                dan[ i ] += J * dAN[ I ] * Finv[ dAN.size( ) * I + i ];

                ddanddAN[ i ][ I ] += J * Finv[ dAN.size( ) * I + i ];

                for ( unsigned int a = 0; a < dAN.size( ); a++ ){

                    for ( unsigned int A = 0; A < dAN.size( ); A++ ){

                        ddandF[ i ][ dAN.size( ) * a + A ] += J * dAN[ I ] * ( Finv[ dAN.size( ) * I + i ] * Finv[ dAN.size( ) * A + a ] - Finv[ dAN.size( ) * A + i ] * Finv[ dAN.size( ) * I + a ] );

                        d2dandFddAN[ i ][ dAN.size( ) * dAN.size( ) * a + dAN.size( ) * A + I ] += J * ( Finv[ dAN.size( ) * I + i ] * Finv[ dAN.size( ) * A + a ] - Finv[ dAN.size( ) * A + i ] * Finv[ dAN.size( ) * I + a ] );

                        for ( unsigned int b = 0; b < dAN.size( ); b++ ){

                            for ( unsigned int B = 0; B < dAN.size( ); B++ ){

                                d2dandFdF[ i ][ dAN.size( ) * dAN.size( ) * dAN.size( ) * a + dAN.size( ) * dAN.size( ) * A + dAN.size( ) * b + B ]
                                    += J * dAN[ I ] * Finv[ dAN.size( ) * B + b ] *  ( Finv[ dAN.size( ) * I + i ] * Finv[ dAN.size( ) * A + a ] - Finv[ dAN.size( ) * A + i ] * Finv[ dAN.size( ) * I + a ] )
                                     - J * dAN[ I ] * ( Finv[ dAN.size( ) * I + b ] * Finv[ dAN.size( ) * B + i ] * Finv[ dAN.size( ) * A + a ] + Finv[ dAN.size( ) * I + i ] * Finv[ dAN.size( ) * A + b ] * Finv[ dAN.size( ) * B + a ]
                                                      - Finv[ dAN.size( ) * A + b ] * Finv[ dAN.size( ) * B + i ] * Finv[ dAN.size( ) * I + a ] - Finv[ dAN.size( ) * A + i ] * Finv[ dAN.size( ) * I + b ] * Finv[ dAN.size( ) * B + a ] );

                            }

                        }

                    }

                }

            }

        }

        return NULL;

    }

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatType &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &overlap ){
        /*!
         * Compute the amount that a point on the local particle overlaps with the non-local particle. For now, we assume
         * a micromorphic theory of degree 1 meaning that for the local particle
         * 
         * \f$ \xi_i = \chi_{iI} \Xi_I\f$
         * 
         * and for the non-local particle
         * 
         * \f$ \xi_i^{NL} = \chi_{iI}^{NL} \Xi_I = \left(\chi_{iI} + \chi_{iI,J} dX_J\right) \Xi_I\f$
         * 
         * where
         * 
         * \f$ dX_I = Xi_I^1 + D_I - Xi_I^2 \f$
         * 
         * So the overlap vector can be computed as
         * 
         * 
         * 
         * \param &Xi_1: The local micro relative position vector.
         * \param &Xi_2: The non-local micro relative position vector.
         * \param &D: The initial spacing between the particles
         * \param &F: The deformation gradient
         * \param &chi: The micro deformation tensor
         * \param &gradChi: The gradient of the micro deformation tensor w.r.t. the reference spatial position
         * \param &overlap: The overlap vector
         */

        return NULL;

    }

    errorOut computeOverlapDistanceLagrangian( const floatVector &X, const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatType &L ){
        /*!
         * Compute the Lagrangian for the computation of the amount that the point \f$\xi_t\f$ in the local particle
         * is overlapping it's non-local neighbor.
         * 
         * \f$ \mathcal{L} = \frac{1}{2} \left( \chi^{nl}_{iI} \Xi_I - \xi_i \right)\left( \chi^{nl}_{iJ} \Xi_J - \xi_i\right) - \lambda \left( \Xi_I \Xi_I - 1\right)\f$
         * 
         * \param &X: The vector of unknowns. The first values are the current estimate of \f$\Xi\f$ and the final value is the Lagrange multiplier \f$\lambda\f$.
         * \param &chi_nl: The non-local micro-deformation tensor
         * \param &xi_t: The target point in the local particle
         * \param &R_nl: The non-local particle radius in the reference configuration
         * \param &L: The value of the Lagrangian
         */

        unsigned int Xsize = X.size( );

        if ( Xsize < ( xi_t.size( ) + 1 ) ){

            return new errorNode( __func__, "X has a size of " + std::to_string( Xsize ) + " and should have a size of " + std::to_string( xi_t.size( ) + 1 ) );

        }

        if ( chi_nl.size( ) != ( xi_t.size( ) * xi_t.size( ) ) ){

            return new errorNode( __func__, "chi_nl has a size of " + std::to_string( chi_nl.size( ) ) + " and should have a size of " + std::to_string( xi_t.size( ) * xi_t.size( ) ) );

        }

        floatVector Xi( X.begin( ), X.begin( ) + xi_t.size( ) );
        floatType lambda = X.back( );

        floatVector d = -xi_t;
        for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

            for ( unsigned int I = 0; I < xi_t.size( ); I++ ){

                d[ i ] += chi_nl[ xi_t.size( ) * i + I ] * Xi[ I ];

            }

        }
        
        L = 0.5 * vectorTools::dot( d, d ) - lambda * ( vectorTools::dot( Xi, Xi ) - R_nl );

        return NULL;

    }

    errorOut computeOverlapDistanceLagrangian( const floatVector &X, const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatType &L,
                                               floatVector &dLdX, floatVector &dLdchi_nl, floatVector &dLdxi_t, floatType &dLdR_nl ){
        /*!
         * Compute the Lagrangian for the computation of the amount that the point \f$\xi_t\f$ in the local particle
         * is overlapping it's non-local neighbor.
         * 
         * \f$ \mathcal{L} = \frac{1}{2} \left( \chi^{nl}_{iI} \Xi_I - \xi_i \right)\left( \chi^{nl}_{iJ} \Xi_J - \xi_i\right) - \lambda \left( \Xi_I \Xi_I - 1\right)\f$
         * 
         * \param &X: The vector of unknowns. The first values are the current estimate of \f$\Xi\f$ and the final value is the Lagrange multiplier \f$\lambda\f$.
         * \param &chi_nl: The non-local micro-deformation tensor
         * \param &xi_t: The target point in the local particle
         * \param &R_nl: The radius of the non-local particle in the reference configuration
         * \param &L: The value of the Lagrangian
         * \param &dLdX: The gradient of the Lagrangian w.r.t. the unknown vector
         * \param &dLdchi_nl: The gradient of the Lagrangian w.r.t. the non-local micro deformation tensor
         * \param &dLdxi_t: The gradient of the Lagrangian w.r.t. the target point
         * \param &dLdR_nl: The gradient of the Lagrangian w.r.t. the reference configuration non-local particle radius
         */

        unsigned int Xsize = X.size( );

        if ( Xsize < ( xi_t.size( ) + 1 ) ){

            return new errorNode( __func__, "X has a size of " + std::to_string( Xsize ) + " and should have a size of " + std::to_string( xi_t.size( ) + 1 ) );

        }

        if ( chi_nl.size( ) != ( xi_t.size( ) * xi_t.size( ) ) ){

            return new errorNode( __func__, "chi_nl has a size of " + std::to_string( chi_nl.size( ) ) + " and should have a size of " + std::to_string( xi_t.size( ) * xi_t.size( ) ) );

        }

        floatVector Xi( X.begin( ), X.begin( ) + xi_t.size( ) );
        floatType lambda = X.back( );

        floatVector d = -xi_t;
        for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

            for ( unsigned int I = 0; I < xi_t.size( ); I++ ){

                d[ i ] += chi_nl[ xi_t.size( ) * i + I ] * Xi[ I ];

            }

        }
        
        L = 0.5 * vectorTools::dot( d, d ) - lambda * ( vectorTools::dot( Xi, Xi ) - R_nl );

        dLdX = floatVector( X.size( ), 0 );

        dLdchi_nl = floatVector( chi_nl.size( ), 0 );

        dLdxi_t = -d;

        dLdR_nl = lambda;

        for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

            for ( unsigned int I = 0; I < xi_t.size( ); I++ ){

                dLdX[ I ] += chi_nl[ xi_t.size( ) * i + I ] * d[ i ];

                dLdchi_nl[ xi_t.size( ) * i + I ] += Xi[ I ] * d[ i ];

            }

            dLdX[ i ] -= 2 * lambda * Xi[ i ];

        }

        dLdX[ Xi.size( ) ] -= ( vectorTools::dot( Xi, Xi ) - R_nl );

        return NULL;

    }

    errorOut computeOverlapDistanceLagrangian( const floatVector &X, const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatType &L,
                                               floatVector &dLdX, floatVector &dLdchi_nl, floatVector &dLdxi_t, floatType &dLdR_nl,
                                               floatVector &d2LdXdX, floatVector &d2LdXdchi_nl, floatVector &d2LdXdxi_t, floatVector &d2LdXdR_nl,
                                               floatVector &d2Ldchi_nldchi_nl, floatVector &d2Ldchi_nldxi_t, floatVector &d2Ldchi_nldR_nl,
                                               floatVector &d2Ldxi_tdxi_t, floatVector &d2Ldxi_tdR_nl,
                                               floatType &d2LdR_nldR_nl ){
        /*!
         * Compute the Lagrangian for the computation of the amount that the point \f$\xi_t\f$ in the local particle
         * is overlapping it's non-local neighbor.
         * 
         * \f$ \mathcal{L} = \frac{1}{2} \left( \chi^{nl}_{iI} \Xi_I - \xi_i \right)\left( \chi^{nl}_{iJ} \Xi_J - \xi_i\right) - \lambda \left( \Xi_I \Xi_I - 1\right)\f$
         * 
         * \param &X: The vector of unknowns. The first values are the current estimate of \f$\Xi\f$ and the final value is the Lagrange multiplier \f$\lambda\f$.
         * \param &chi_nl: The non-local micro-deformation tensor
         * \param &xi_t: The target point in the local particle
         * \param &R_nl: The radius of the non-local particle in the reference configuration
         * \param &L: The value of the Lagrangian
         * \param &dLdX: The gradient of the Lagrangian w.r.t. the unknown vector
         * \param &dLdchi_nl: The gradient of the Lagrangian w.r.t. the non-local micro deformation tensor
         * \param &dLdxi_t: The gradient of the Lagrangian w.r.t. the target point
         * \param &dLdR_nl: The gradient of the Lagrangian w.r.t. the reference configuration non-local particle radius
         * \param &d2LdXdX: The second gradient of the Lagrangian w.r.t. the solution vector
         * \param &d2LdXdchi_nl: The second gradient of the Lagrangian w.r.t. the solution vector and the non-local micro deformation tensor
         * \param &d2LdXdxi_t: The second gradient of the Lagrangian w.r.t. the solution vector and the target point
         * \param &d2LdXdR_nl: The second gradient of the Lagrangina w.r.t. the solution vector and the non-local radius in the reference configuration
         * \param &d2Ldchi_nldchi_nl: The second gradient of the Lagrangian w.r.t. the non-local micro deformation tensor
         * \param &d2Ldchi_nldxi_t: The second gradient of the Lagrangian w.r.t. the non-local micro deformation tensor and the target point
         * \param &d2Ldchi_nldR_nl: The second gradient of the Lagrangina w.r.t. the non-local micro deformation tensor and the non-local radius in the reference configuration
         * \param &d2Ldxi_tdxi_t: The second gradient of the Lagrangian w.r.t. the target point
         * \param &d2Ldxi_tdR_nl: The second gradient of the Lagrangina w.r.t. the target point and the non-local radius in the reference configuration
         * \param &d2LdR_nldR_nl: The second gradient of the Lagrangina w.r.t. the non-local radius in the reference configuration
         */

        unsigned int Xsize = X.size( );

        if ( Xsize < ( xi_t.size( ) + 1 ) ){

            return new errorNode( __func__, "X has a size of " + std::to_string( Xsize ) + " and should have a size of " + std::to_string( xi_t.size( ) + 1 ) );

        }

        if ( chi_nl.size( ) != ( xi_t.size( ) * xi_t.size( ) ) ){

            return new errorNode( __func__, "chi_nl has a size of " + std::to_string( chi_nl.size( ) ) + " and should have a size of " + std::to_string( xi_t.size( ) * xi_t.size( ) ) );

        }

        floatVector Xi( X.begin( ), X.begin( ) + xi_t.size( ) );
        floatType lambda = X.back( );

        floatVector d = -xi_t;
        for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

            for ( unsigned int I = 0; I < xi_t.size( ); I++ ){

                d[ i ] += chi_nl[ xi_t.size( ) * i + I ] * Xi[ I ];

            }

        }
        
        L = 0.5 * vectorTools::dot( d, d ) - lambda * ( vectorTools::dot( Xi, Xi ) - R_nl );

        dLdX = floatVector( X.size( ), 0 );

        dLdchi_nl = floatVector( chi_nl.size( ), 0 );

        dLdxi_t = -d;

        dLdR_nl = lambda;

        d2LdXdX = floatVector( X.size( ) * X.size( ), 0 );

        d2LdXdchi_nl = floatVector( X.size( ) * chi_nl.size( ), 0 );

        d2LdXdxi_t = floatVector( X.size( ) * xi_t.size( ), 0 );

        d2Ldchi_nldchi_nl = floatVector( chi_nl.size( ) * chi_nl.size( ), 0 );

        d2Ldchi_nldxi_t = floatVector( chi_nl.size( ) * xi_t.size( ), 0 );

        d2Ldxi_tdxi_t = floatVector( xi_t.size( ) * xi_t.size( ), 0 );

        d2LdXdR_nl = floatVector( X.size( ), 0 );

        d2Ldchi_nldR_nl = floatVector( chi_nl.size( ), 0 );

        d2Ldxi_tdR_nl = floatVector( xi_t.size( ), 0 );

        d2LdR_nldR_nl = 0;

        floatVector eye( chi_nl.size( ), 0 );

        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

            d2Ldxi_tdxi_t[ xi_t.size( ) * i + i ] = 1;

            for ( unsigned int I = 0; I < xi_t.size( ); I++ ){

                dLdX[ I ] += chi_nl[ xi_t.size( ) * i + I ] * d[ i ];

                dLdchi_nl[ xi_t.size( ) * i + I ] += Xi[ I ] * d[ i ];

                d2LdXdxi_t[ xi_t.size( ) * i + I ] -= chi_nl[ xi_t.size( ) * I + i ];

                for ( unsigned int J = 0; J < xi_t.size( ); J++ ){

                    d2LdXdX[ X.size( ) * I + J ] += chi_nl[ xi_t.size( ) * i + I ] * chi_nl[ xi_t.size( ) * i + J ];

                    d2LdXdchi_nl[ xi_t.size( ) * xi_t.size( ) * I + xi_t.size( ) * i + J ] += d[ i ] * eye[ xi_t.size( ) * I + J ] + chi_nl[ xi_t.size( ) * i + I ] * Xi[ J ];

                    d2Ldchi_nldxi_t[ xi_t.size( ) * xi_t.size( ) * i + xi_t.size( ) * I + J ] -= Xi[ I ] * eye[ xi_t.size( ) * i + J ];

                    for ( unsigned int j = 0; j < xi_t.size( ); j++ ){

                        d2Ldchi_nldchi_nl[ xi_t.size( ) * xi_t.size( ) * xi_t.size( ) * i + xi_t.size( ) * xi_t.size( ) * I + j * xi_t.size( ) + J ] += Xi[ I ] * eye[ xi_t.size( ) * i + j ] * Xi[ J ];

                    }

                }

            }

            dLdX[ i ] -= 2 * lambda * Xi[ i ];

            d2LdXdX[ X.size( ) * i + i ] -= 2 * lambda;

            d2LdXdX[ X.size( ) * i + X.size( ) - 1 ] -= 2 * Xi[ i ];
            d2LdXdX[ X.size( ) * ( X.size( ) - 1 ) + i ] -= 2 * Xi[ i ];

        }

        dLdX[ Xi.size( ) ] -= ( vectorTools::dot( Xi, Xi ) - R_nl );

        d2LdXdR_nl[ X.size( ) - 1 ] = 1;

        return NULL;

    }

    errorOut solveOverlapDistance( const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatVector &d,
                                   const floatType tolr, const floatType tola, const unsigned int max_iteration,
                                   const unsigned int max_ls, const floatType alpha_ls ){
        /*!
         * Solve for the overlap distance where \f$\xi_t\f$ is known to be inside of the non-local particle
         * 
         * \param &chi_nl: The non-local micro-deformation tensor
         * \param &xi_t: The position inside of the non-local particle
         * \param &R_nl: The non-local particle radius in the reference configuration
         * \param &d: The distance vector going from the solved point on the surface of the non-local particle
         *     to \f$\xi_t\f$.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param max_iteration: The maximum number of iterations
         * \param max_ls: The maximum number of line-search iterations
         * \param alpha_ls: The alpha parameter for the line-search
         */

        floatVector inv_chi_nl = vectorTools::inverse( chi_nl, xi_t.size( ), xi_t.size( ) );

        floatVector Xi_t( xi_t.size( ), 0 );

        for ( unsigned int I = 0; I < xi_t.size( ); I++ ){

            for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

                Xi_t[ I ] += inv_chi_nl[ xi_t.size( ) * I + i ] * xi_t[ i ];

            }

        }

        floatVector lagrange( 1, 1 );

        floatVector X = vectorTools::appendVectors( { Xi_t, lagrange } );

        floatType L;

        floatVector dLdX, dLdchi_nl, dLdxi_t;

        floatType dLdR_nl;

        floatVector d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t, d2LdXdR_nl, d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t, d2Ldchi_nldR_nl, d2Ldxi_tdxi_t, d2Ldxi_tdR_nl;

        floatType d2LdR_nldR_nl;

        errorOut error = computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl, L,
                                                           dLdX, dLdchi_nl, dLdxi_t, dLdR_nl,
                                                           d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t, d2LdXdR_nl,
                                                           d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t, d2Ldchi_nldR_nl,
                                                           d2Ldxi_tdxi_t, d2Ldxi_tdR_nl,
                                                           d2LdR_nldR_nl );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the initial Lagrangian" );

            result->addNext( error );

            return result;

        }

        floatType R = vectorTools::l2norm( dLdX );

        floatType Rp = R;

        floatType tol = tolr * R + tola;

        unsigned int num_iteration = 0;

        floatVector dX;

        while ( ( num_iteration < max_iteration ) && ( R > tol ) ){

            std::cout << "  dLdX: "; vectorTools::print( dLdX );
            std::cout << "  d2LdXdX:\n"; vectorTools::print( vectorTools::inflate( d2LdXdX, dLdX.size( ), dLdX.size( ) ) );

            unsigned int rank;

            dX = -vectorTools::solveLinearSystem( d2LdXdX, dLdX, dLdX.size( ), dLdX.size( ), rank );

            std::cout << "  dX: "; vectorTools::print( dX );
            std::cout << "  rank: " << rank << "\n";

            floatType lambda = 1;

            error = computeOverlapDistanceLagrangian( X + lambda * dX, chi_nl, xi_t, R_nl, L,
                                                      dLdX, dLdchi_nl, dLdxi_t, dLdR_nl,
                                                      d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t, d2LdXdR_nl,
                                                      d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t, d2Ldchi_nldR_nl,
                                                      d2Ldxi_tdxi_t, d2Ldxi_tdR_nl,
                                                      d2LdR_nldR_nl );

            if ( error ){
    
                errorOut result = new errorNode( __func__, "Error in the computation of iteration " + std::to_string( num_iteration ) + " Lagrangian" );
    
                result->addNext( error );
    
                return result;
    
            }

            unsigned int num_ls = 0;

            R = vectorTools::l2norm( dLdX );

            while ( ( num_ls < max_ls ) && ( R > ( 1 - alpha_ls ) * Rp ) ){ 

                lambda *= 0.5;

                error = computeOverlapDistanceLagrangian( X + lambda * dX, chi_nl, xi_t, R_nl, L,
                                                          dLdX, dLdchi_nl, dLdxi_t, dLdR_nl,
                                                          d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t, d2LdXdR_nl,
                                                          d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t, d2Ldchi_nldR_nl,
                                                          d2Ldxi_tdxi_t, d2Ldxi_tdR_nl,
                                                          d2LdR_nldR_nl );
    
                if ( error ){
        
                    errorOut result = new errorNode( __func__, "Error in the " + std::to_string( num_ls ) + " line search of iteration " + std::to_string( num_iteration ) + " initial Lagrangian" );
        
                    result->addNext( error );
        
                    return result;
        
                }

                R = vectorTools::l2norm( dLdX );

                num_ls++;

            }

            X += lambda * dX;

            Rp = R;

            num_iteration++;

        }

        if ( R > tol ){

            return new errorNode( __func__, "The optimizer did not converge" );

        }

        return NULL;

    }

    errorOut solveOverlapDistance( const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatVector &d,
                                   floatMatrix &dddchi_nl, floatMatrix &dddxi_t, floatVector &dddR_nl,
                                   const floatType tolr, const floatType tola, const unsigned int max_iteration,
                                   const unsigned int max_ls, const floatType alpha_ls );

    errorOut solveOverlapDistance( const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatVector &d,
                                   floatMatrix &dddchi_nl, floatMatrix &dddxi_t, floatVector &dddR_nl,
                                   floatMatrix &d2ddchi_nldchi_nl, floatMatrix &d2ddchi_nldxi_t, floatMatrix &d2ddchi_nldR_nl,
                                   floatMatrix &d2ddxi_tdxi_t, floatMatrix &d2ddxi_tdR_nl,
                                   floatVector &d2ddR_nldR_nl,
                                   const floatType tolr, const floatType tola, const unsigned int max_iteration,
                                   const unsigned int max_ls, const floatType alpha_ls );


}
