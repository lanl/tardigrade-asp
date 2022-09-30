/**
  * \file test_traction_separation.cpp
  *
  * Tests for the traction separation support module
  */

#include<traction_separation.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_traction_separation
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tractionSeparation::floatType floatType; //!< Redefinition for the float type
typedef tractionSeparation::floatVector floatVector; //1< Redefinition for the float vector
typedef tractionSeparation::floatMatrix floatMatrix; //1< Redefinition for the float matrix

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( test_computeCurrentDistance ){
    /*!
     * Tests for the computation of the current distance
     */

    floatVector Xi_1 = { 0.69646919, 0.28613933, 0.22685145 };
    floatVector Xi_2 = { 0.55131477, 0.71946897, 0.42310646 };
    floatVector D    = { 0.9807642 , 0.68482974, 0.4809319  };
    
    floatVector F    = { 0.39211752, 0.34317802, 0.72904971,
                         0.43857224, 0.0596779 , 0.39804426,
                         0.73799541, 0.18249173, 0.17545176 };

    floatVector chi  = { 0.53155137, 0.53182759, 0.63440096,
                         0.84943179, 0.72445532, 0.61102351,
                         0.72244338, 0.32295891, 0.36178866 };

    floatVector gradChi = { 0.22826323, 0.29371405, 0.63097612,
                            0.09210494, 0.43370117, 0.43086276,
                            0.4936851 , 0.42583029, 0.31226122,
                            0.42635131, 0.89338916, 0.94416002,
                            0.50183668, 0.62395295, 0.1156184 ,
                            0.31728548, 0.41482621, 0.86630916,
                            0.25045537, 0.48303426, 0.98555979,
                            0.51948512, 0.61289453, 0.12062867,
                            0.8263408 , 0.60306013, 0.54506801 };

    floatVector d_answer = { 1.85403834, 2.3121126 , 2.48987739 };

    floatVector d;

    BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi, gradChi, d ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d, d_answer ) );

    d.clear( );

    floatMatrix dddXi_1, dddXi_2, dddD, dddF, dddChi, dddGradChi;

    BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi, gradChi, d, dddXi_1, dddXi_2, dddD, dddF, dddChi, dddGradChi ) );

    floatVector d_2;
    floatMatrix dddXi_1_2, dddXi_2_2, dddD_2, dddF_2, dddChi_2, dddGradChi_2;
    floatMatrix d2ddFdXi_1,       d2ddFdXi_2,       d2ddFdD;
    floatMatrix d2ddChidXi_1,     d2ddChidXi_2,     d2ddChidD;
    floatMatrix d2ddGradChidXi_1, d2ddGradChidXi_2, d2ddGradChidD;

    BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi, gradChi, d_2, dddXi_1_2, dddXi_2_2, dddD_2, dddF_2, dddChi_2, dddGradChi_2,
                                                              d2ddFdXi_1,       d2ddFdXi_2,       d2ddFdD,
                                                              d2ddChidXi_1,     d2ddChidXi_2,     d2ddChidD,
                                                              d2ddGradChidXi_1, d2ddGradChidXi_2, d2ddGradChidD ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d_2, d_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_1,    dddXi_1_2    ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_2,    dddXi_2_2    ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( dddD,       dddD_2       ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( dddF,       dddF_2       ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( dddChi,     dddChi_2     ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( dddGradChi, dddGradChi_2 ) );

    floatMatrix dddXi_1_answer( d_answer.size( ), floatVector( Xi_1.size( ), 0 ) );
    floatMatrix dddXi_2_answer( d_answer.size( ), floatVector( Xi_2.size( ), 0 ) );
    floatMatrix dddD_answer( d_answer.size( ), floatVector( D.size( ), 0 ) );
    floatMatrix dddF_answer( d_answer.size( ), floatVector( F.size( ), 0 ) );
    floatMatrix dddChi_answer( d_answer.size( ), floatVector( chi.size( ), 0 ) );
    floatMatrix dddGradChi_answer( d_answer.size( ), floatVector( gradChi.size( ), 0 ) );
    floatMatrix d2ddFdXi_1_answer( d_answer.size( ), floatVector( F.size( ) * Xi_1.size( ), 0 ) );
    floatMatrix d2ddFdXi_2_answer( d_answer.size( ), floatVector( F.size( ) * Xi_2.size( ), 0 ) );
    floatMatrix d2ddFdD_answer( d_answer.size( ), floatVector( F.size( ) * D.size( ), 0 ) );
    floatMatrix d2ddChidXi_1_answer( d_answer.size( ), floatVector( chi.size( ) * Xi_1.size( ), 0 ) );
    floatMatrix d2ddChidXi_2_answer( d_answer.size( ), floatVector( chi.size( ) * Xi_2.size( ), 0 ) );
    floatMatrix d2ddChidD_answer( d_answer.size( ), floatVector( chi.size( ) * D.size( ), 0 ) );
    floatMatrix d2ddGradChidXi_1_answer( d_answer.size( ), floatVector( gradChi.size( ) * Xi_1.size( ), 0 ) );
    floatMatrix d2ddGradChidXi_2_answer( d_answer.size( ), floatVector( gradChi.size( ) * Xi_2.size( ), 0 ) );
    floatMatrix d2ddGradChidD_answer( d_answer.size( ), floatVector( gradChi.size( ) * D.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){
        
        floatVector delta( Xi_1.size( ), 0 );
        delta[ i ] = eps * std::abs( Xi_1[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1 + delta, Xi_2, D, F, chi, gradChi, dp ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1 - delta, Xi_2, D, F, chi, gradChi, dm ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            dddXi_1_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dddXi_1p, dddXi_1m;
        floatMatrix dddXi_2p, dddXi_2m;
        floatMatrix dddDp, dddDm;
        floatMatrix dddFp, dddFm;
        floatMatrix dddChip, dddChim;
        floatMatrix dddGradChip, dddGradChim; 

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1 + delta, Xi_2, D, F, chi, gradChi, dp, dddXi_1p, dddXi_2p, dddDp, dddFp, dddChip, dddGradChip ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1 - delta, Xi_2, D, F, chi, gradChi, dm, dddXi_1m, dddXi_2m, dddDm, dddFm, dddChim, dddGradChim ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_1_answer[ j ][ i ], ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < d_answer.size( ); k++ ){

                for ( unsigned int K = 0; K < d_answer.size( ); K++ ){

                    d2ddFdXi_1_answer[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + i ]
                        += ( dddFp[ j ][ d_answer.size( ) * k + K ] - dddFm[ j ][ d_answer.size( ) * k + K ] ) / ( 2 * delta [ i ] );
                    d2ddChidXi_1_answer[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + i ]
                        += ( dddChip[ j ][ d_answer.size( ) * k + K ] - dddChim[ j ][ d_answer.size( ) * k + K ] ) / ( 2 * delta [ i ] );

                    for ( unsigned int L = 0; L < d_answer.size( ); L++ ){

                        d2ddGradChidXi_1_answer[ j ][ d_answer.size( ) * d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * d_answer.size( ) * K + d_answer.size( ) * L + i ]
                            += ( dddGradChip[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + L ] - dddGradChim[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + L ] ) / ( 2 * delta[ i ] );

                    }

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_1,          dddXi_1_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXi_1,       d2ddFdXi_1_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidXi_1,     d2ddChidXi_1_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidXi_1, d2ddGradChidXi_1_answer ) );

    for ( unsigned int i = 0; i < Xi_2.size( ); i++ ){
        
        floatVector delta( Xi_2.size( ), 0 );
        delta[ i ] = eps * std::abs( Xi_2[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2 + delta, D, F, chi, gradChi, dp ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2 - delta, D, F, chi, gradChi, dm ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            dddXi_2_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dddXi_1p, dddXi_1m;
        floatMatrix dddXi_2p, dddXi_2m;
        floatMatrix dddDp, dddDm;
        floatMatrix dddFp, dddFm;
        floatMatrix dddChip, dddChim;
        floatMatrix dddGradChip, dddGradChim; 

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2 + delta, D, F, chi, gradChi, dp, dddXi_1p, dddXi_2p, dddDp, dddFp, dddChip, dddGradChip ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2 - delta, D, F, chi, gradChi, dm, dddXi_1m, dddXi_2m, dddDm, dddFm, dddChim, dddGradChim ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_2_answer[ j ][ i ], ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < d_answer.size( ); k++ ){

                for ( unsigned int K = 0; K < d_answer.size( ); K++ ){

                    d2ddFdXi_2_answer[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + i ]
                        += ( dddFp[ j ][ d_answer.size( ) * k + K ] - dddFm[ j ][ d_answer.size( ) * k + K ] ) / ( 2 * delta [ i ] );
                    d2ddChidXi_2_answer[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + i ]
                        += ( dddChip[ j ][ d_answer.size( ) * k + K ] - dddChim[ j ][ d_answer.size( ) * k + K ] ) / ( 2 * delta [ i ] );

                    for ( unsigned int L = 0; L < d_answer.size( ); L++ ){

                        d2ddGradChidXi_2_answer[ j ][ d_answer.size( ) * d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * d_answer.size( ) * K + d_answer.size( ) * L + i ]
                            += ( dddGradChip[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + L ] - dddGradChim[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + L ] ) / ( 2 * delta[ i ] );

                    }

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_2,          dddXi_2_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXi_2,       d2ddFdXi_2_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidXi_2,     d2ddChidXi_2_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidXi_2, d2ddGradChidXi_2_answer ) );

    for ( unsigned int i = 0; i < D.size( ); i++ ){
        
        floatVector delta( D.size( ), 0 );
        delta[ i ] = eps * std::abs( D[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D + delta, F, chi, gradChi, dp ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D - delta, F, chi, gradChi, dm ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            dddD_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dddXi_1p, dddXi_1m;
        floatMatrix dddXi_2p, dddXi_2m;
        floatMatrix dddDp, dddDm;
        floatMatrix dddFp, dddFm;
        floatMatrix dddChip, dddChim;
        floatMatrix dddGradChip, dddGradChim; 

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D + delta, F, chi, gradChi, dp, dddXi_1p, dddXi_2p, dddDp, dddFp, dddChip, dddGradChip ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D - delta, F, chi, gradChi, dm, dddXi_1m, dddXi_2m, dddDm, dddFm, dddChim, dddGradChim ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( dddD_answer[ j ][ i ], ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < d_answer.size( ); k++ ){

                for ( unsigned int K = 0; K < d_answer.size( ); K++ ){

                    d2ddFdD_answer[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + i ]
                        += ( dddFp[ j ][ d_answer.size( ) * k + K ] - dddFm[ j ][ d_answer.size( ) * k + K ] ) / ( 2 * delta [ i ] );
                    d2ddChidD_answer[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + i ]
                        += ( dddChip[ j ][ d_answer.size( ) * k + K ] - dddChim[ j ][ d_answer.size( ) * k + K ] ) / ( 2 * delta [ i ] );

                    for ( unsigned int L = 0; L < d_answer.size( ); L++ ){

                        d2ddGradChidD_answer[ j ][ d_answer.size( ) * d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * d_answer.size( ) * K + d_answer.size( ) * L + i ]
                            += ( dddGradChip[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + L ] - dddGradChim[ j ][ d_answer.size( ) * d_answer.size( ) * k + d_answer.size( ) * K + L ] ) / ( 2 * delta[ i ] );

                    }

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddD,          dddD_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdD,       d2ddFdD_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddChidD,     d2ddChidD_answer ) );
    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddGradChidD, d2ddGradChidD_answer ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){
        
        floatVector delta( F.size( ), 0 );
        delta[ i ] = eps * std::abs( F[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F + delta, chi, gradChi, dp ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F - delta, chi, gradChi, dm ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            dddF_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dddXi_1p, dddXi_1m;
        floatMatrix dddXi_2p, dddXi_2m;
        floatMatrix dddDp, dddDm;
        floatMatrix dddFp, dddFm;
        floatMatrix dddChip, dddChim;
        floatMatrix dddGradChip, dddGradChim; 

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F + delta, chi, gradChi, dp, dddXi_1p, dddXi_2p, dddDp, dddFp, dddChip, dddGradChip ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F - delta, chi, gradChi, dm, dddXi_1m, dddXi_2m, dddDm, dddFm, dddChim, dddGradChim ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( dddF_answer[ j ][ i ], ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        BOOST_CHECK( vectorTools::fuzzyEquals( dddFp, dddFm ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( dddChip, dddChim ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( dddGradChip, dddGradChim ) );

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddF, dddF_answer ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){
        
        floatVector delta( chi.size( ), 0 );
        delta[ i ] = eps * std::abs( chi[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi + delta, gradChi, dp ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi - delta, gradChi, dm ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            dddChi_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dddXi_1p, dddXi_1m;
        floatMatrix dddXi_2p, dddXi_2m;
        floatMatrix dddDp, dddDm;
        floatMatrix dddFp, dddFm;
        floatMatrix dddChip, dddChim;
        floatMatrix dddGradChip, dddGradChim;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi + delta, gradChi, dp, dddXi_1p, dddXi_2p, dddDp, dddFp, dddChip, dddGradChip ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi - delta, gradChi, dm, dddXi_1m, dddXi_2m, dddDm, dddFm, dddChim, dddGradChim ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( dddChi_answer[ j ][ i ], ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        BOOST_CHECK( vectorTools::fuzzyEquals( dddFp, dddFm ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( dddChip, dddChim ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( dddGradChip, dddGradChim ) );

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddChi, dddChi_answer ) );

    for ( unsigned int i = 0; i < gradChi.size( ); i++ ){
        
        floatVector delta( gradChi.size( ), 0 );
        delta[ i ] = eps * std::abs( gradChi[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi, gradChi + delta, dp ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi, gradChi - delta, dm ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            dddGradChi_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dddXi_1p, dddXi_1m;
        floatMatrix dddXi_2p, dddXi_2m;
        floatMatrix dddDp, dddDm;
        floatMatrix dddFp, dddFm;
        floatMatrix dddChip, dddChim;
        floatMatrix dddGradChip, dddGradChim;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi, gradChi + delta, dp, dddXi_1p, dddXi_2p, dddDp, dddFp, dddChip, dddGradChip ) );
        BOOST_CHECK( !tractionSeparation::computeCurrentDistance( Xi_1, Xi_2, D, F, chi, gradChi - delta, dm, dddXi_1m, dddXi_2m, dddDm, dddFm, dddChim, dddGradChim ) );

        for ( unsigned int j = 0; j < d_answer.size( ); j++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( dddGradChi_answer[ j ][ i ], ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        BOOST_CHECK( vectorTools::fuzzyEquals( dddFp, dddFm ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( dddChip, dddChim ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( dddGradChip, dddGradChim ) );

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddGradChi, dddGradChi_answer ) );

}
