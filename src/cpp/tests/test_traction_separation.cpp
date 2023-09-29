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

BOOST_AUTO_TEST_CASE( test_computeCurrentDistanceGeneral ){
    /*!
     * Test the computation of a distance between two points in a more generalized fashion
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

    floatVector chiNL = { 0.88594794, 0.07791236, 0.97964616,
                          0.24767146, 0.75288472, 0.52667564,
                          0.90755375, 0.8840703 , 0.08926896 };

    floatVector d_answer = { 1.02803094, 0.58567184, 1.42330265 };

    floatVector d;

    BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi, chiNL, d ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d, d_answer ) );

    floatVector d_2;

    floatMatrix dddXi_1, dddXi_2, dddD, dddF, dddchi, dddchiNL;

    BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi, chiNL, d_2,
                                                                     dddXi_1, dddXi_2, dddD, dddF, dddchi, dddchiNL ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d_2, d_answer ) );

    floatVector d_3;

    floatMatrix dddXi_1_3, dddXi_2_3, dddD_3, dddF_3, dddchi_3, dddchiNL_3,
                d2ddFdXi_1, d2ddchidXi_1, d2ddFdXi_2, d2ddchiNLdXi_2, d2ddFdD;

    BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi, chiNL, d_3,
                                                                     dddXi_1_3, dddXi_2_3, dddD_3, dddF_3, dddchi_3, dddchiNL_3,
                                                                     d2ddFdXi_1, d2ddchidXi_1, d2ddFdXi_2, d2ddchiNLdXi_2, d2ddFdD ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d_3, d_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_1_3, dddXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_2_3, dddXi_2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddD_3, dddD ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddF_3, dddF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddchi_3, dddchi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dddchiNL_3, dddchiNL ) );

    floatType eps = 1e-6;

    // Test Jacobians w.r.t. the local reference relative position vector
    floatMatrix dddXi_1_answer( d_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix d2ddFdXi_1_answer( d_answer.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d2ddchidXi_1_answer( d_answer.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){

        floatVector deltas( Xi_1.size( ), 0 );

        deltas[ i ] += eps * std::fabs( Xi_1[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1 + deltas, Xi_2, D, F, chi, chiNL, dp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1 - deltas, Xi_2, D, F, chi, chiNL, dm ) );

        floatVector gradCol = ( dp - dm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradCol.size( ); j++ ){

            dddXi_1_answer[ j ][ i ] = gradCol[ j ];

        }

        floatMatrix dddXi_1p, dddXi_1m;

        floatMatrix dddXi_2p, dddXi_2m;

        floatMatrix dddDp, dddDm;

        floatMatrix dddFp, dddFm;

        floatMatrix dddchip, dddchim;

        floatMatrix dddchiNLp, dddchiNLm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1 + deltas, Xi_2, D, F, chi, chiNL, dp,
                                                                         dddXi_1p, dddXi_2p, dddDp, dddFp, dddchip, dddchiNLp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1 - deltas, Xi_2, D, F, chi, chiNL, dm,
                                                                         dddXi_1m, dddXi_2m, dddDm, dddFm, dddchim, dddchiNLm ) );

        floatMatrix gradMat = ( dddXi_1p - dddXi_1m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddXi_2p - dddXi_2m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddDp - dddDm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddFp - dddFm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddFdXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddchip - dddchim ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddchidXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddchiNLp - dddchiNLm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_1, dddXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXi_1, d2ddFdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddchidXi_1, d2ddchidXi_1_answer ) );

    // Test Jacobians w.r.t. the non-local reference relative position vector
    floatMatrix dddXi_2_answer( d_answer.size( ), floatVector( Xi_2.size( ), 0 ) );

    floatMatrix d2ddFdXi_2_answer( d_answer.size( ), floatVector( Xi_2.size( ) * F.size( ), 0 ) );

    floatMatrix d2ddchiNLdXi_2_answer( d_answer.size( ), floatVector( Xi_2.size( ) * chiNL.size( ), 0 ) );

    for ( unsigned int i = 0; i < Xi_2.size( ); i++ ){

        floatVector deltas( Xi_2.size( ), 0 );

        deltas[ i ] += eps * std::fabs( Xi_2[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2 + deltas, D, F, chi, chiNL, dp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2 - deltas, D, F, chi, chiNL, dm ) );

        floatVector gradCol = ( dp - dm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradCol.size( ); j++ ){

            dddXi_2_answer[ j ][ i ] = gradCol[ j ];

        }

        floatMatrix dddXi_1p, dddXi_1m;

        floatMatrix dddXi_2p, dddXi_2m;

        floatMatrix dddDp, dddDm;

        floatMatrix dddFp, dddFm;

        floatMatrix dddchip, dddchim;

        floatMatrix dddchiNLp, dddchiNLm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2 + deltas, D, F, chi, chiNL, dp,
                                                                         dddXi_1p, dddXi_2p, dddDp, dddFp, dddchip, dddchiNLp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2 - deltas, D, F, chi, chiNL, dm,
                                                                         dddXi_1m, dddXi_2m, dddDm, dddFm, dddchim, dddchiNLm ) );

        floatMatrix gradMat = ( dddXi_1p - dddXi_1m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddXi_2p - dddXi_2m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddDp - dddDm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddFp - dddFm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddFdXi_2_answer[ j ][ Xi_2.size( ) * k + i ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddchip - dddchim ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchiNLp - dddchiNLm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddchiNLdXi_2_answer[ j ][ Xi_2.size( ) * k + i ] = gradMat[ j ][ k ];

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddXi_2, dddXi_2_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXi_2, d2ddFdXi_2_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddchiNLdXi_2, d2ddchiNLdXi_2_answer ) );

    // Test Jacobians w.r.t. the non-local reference separation
    floatMatrix dddD_answer( d_answer.size( ), floatVector( D.size( ), 0 ) );

    floatMatrix d2ddFdD_answer( d_answer.size( ), floatVector( F.size( ) * D.size( ), 0 ) );

    for ( unsigned int i = 0; i < D.size( ); i++ ){

        floatVector deltas( D.size( ), 0 );

        deltas[ i ] += eps * std::fabs( D[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D + deltas, F, chi, chiNL, dp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D - deltas, F, chi, chiNL, dm ) );

        floatVector gradCol = ( dp - dm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradCol.size( ); j++ ){

            dddD_answer[ j ][ i ] = gradCol[ j ];

        }

        floatMatrix dddXi_1p, dddXi_1m;

        floatMatrix dddXi_2p, dddXi_2m;

        floatMatrix dddDp, dddDm;

        floatMatrix dddFp, dddFm;

        floatMatrix dddchip, dddchim;

        floatMatrix dddchiNLp, dddchiNLm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D + deltas, F, chi, chiNL, dp,
                                                                         dddXi_1p, dddXi_2p, dddDp, dddFp, dddchip, dddchiNLp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D - deltas, F, chi, chiNL, dm,
                                                                         dddXi_1m, dddXi_2m, dddDm, dddFm, dddchim, dddchiNLm ) );

        floatMatrix gradMat = ( dddXi_1p - dddXi_1m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddXi_2p - dddXi_2m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddDp - dddDm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddFp - dddFm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddFdD_answer[ j ][ D.size( ) * k + i ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddchip - dddchim ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchiNLp - dddchiNLm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddD, dddD_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdD, d2ddFdD_answer ) );

    // Test Jacobians w.r.t. the deformation gradient
    floatMatrix dddF_answer( d_answer.size( ), floatVector( F.size( ), 0 ) );

    d2ddFdXi_1_answer = floatMatrix( d_answer.size( ), floatVector( F.size( ) * Xi_1.size( ), 0 ) );

    d2ddFdXi_2_answer = floatMatrix( d_answer.size( ), floatVector( F.size( ) * Xi_2.size( ), 0 ) );

    d2ddFdD_answer    = floatMatrix( d_answer.size( ), floatVector( F.size( ) * D.size( ), 0 ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){

        floatVector deltas( F.size( ), 0 );

        deltas[ i ] += eps * std::fabs( F[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F + deltas, chi, chiNL, dp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F - deltas, chi, chiNL, dm ) );

        floatVector gradCol = ( dp - dm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradCol.size( ); j++ ){

            dddF_answer[ j ][ i ] = gradCol[ j ];

        }

        floatMatrix dddXi_1p, dddXi_1m;

        floatMatrix dddXi_2p, dddXi_2m;

        floatMatrix dddDp, dddDm;

        floatMatrix dddFp, dddFm;

        floatMatrix dddchip, dddchim;

        floatMatrix dddchiNLp, dddchiNLm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F + deltas, chi, chiNL, dp,
                                                                         dddXi_1p, dddXi_2p, dddDp, dddFp, dddchip, dddchiNLp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F - deltas, chi, chiNL, dm,
                                                                         dddXi_1m, dddXi_2m, dddDm, dddFm, dddchim, dddchiNLm ) );

        floatMatrix gradMat = ( dddXi_1p - dddXi_1m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddFdXi_1_answer[ j ][ Xi_1.size( ) * i + k ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddXi_2p - dddXi_2m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddFdXi_2_answer[ j ][ Xi_2.size( ) * i + k ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddDp - dddDm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddFdD_answer[ j ][ D.size( ) * i + k ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddFp - dddFm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchip - dddchim ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchiNLp - dddchiNLm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddF, dddF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXi_1, d2ddFdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdXi_2, d2ddFdXi_2_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddFdD, d2ddFdD_answer ) );

    // Test Jacobians w.r.t. the local micro-deformation gradient
    floatMatrix dddchi_answer( d_answer.size( ), floatVector( chi.size( ), 0 ) );

    d2ddchidXi_1_answer = floatMatrix( d_answer.size( ), floatVector( chi.size( ) * Xi_1.size( ), 0 ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector deltas( chi.size( ), 0 );

        deltas[ i ] += eps * std::fabs( chi[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi + deltas, chiNL, dp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi - deltas, chiNL, dm ) );

        floatVector gradCol = ( dp - dm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradCol.size( ); j++ ){

            dddchi_answer[ j ][ i ] = gradCol[ j ];

        }

        floatMatrix dddXi_1p, dddXi_1m;

        floatMatrix dddXi_2p, dddXi_2m;

        floatMatrix dddDp, dddDm;

        floatMatrix dddFp, dddFm;

        floatMatrix dddchip, dddchim;

        floatMatrix dddchiNLp, dddchiNLm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi + deltas, chiNL, dp,
                                                                         dddXi_1p, dddXi_2p, dddDp, dddFp, dddchip, dddchiNLp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi - deltas, chiNL, dm,
                                                                         dddXi_1m, dddXi_2m, dddDm, dddFm, dddchim, dddchiNLm ) );

        floatMatrix gradMat = ( dddXi_1p - dddXi_1m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddchidXi_1_answer[ j ][ Xi_1.size( ) * i + k ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddXi_2p - dddXi_2m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddDp - dddDm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddFp - dddFm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchip - dddchim ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchiNLp - dddchiNLm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddchi, dddchi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddchidXi_1, d2ddchidXi_1_answer ) );

    // Test Jacobians w.r.t. the non-local micro-deformation gradient
    floatMatrix dddchiNL_answer( d_answer.size( ), floatVector( chiNL.size( ), 0 ) );

    d2ddchiNLdXi_2 = floatMatrix( d_answer.size( ), floatVector( chiNL.size( ) * Xi_2.size( ), 0 ) );

    for ( unsigned int i = 0; i < chiNL.size( ); i++ ){

        floatVector deltas( chiNL.size( ), 0 );

        deltas[ i ] += eps * std::fabs( chiNL[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi, chiNL + deltas, dp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi, chiNL - deltas, dm ) );

        floatVector gradCol = ( dp - dm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradCol.size( ); j++ ){

            dddchiNL_answer[ j ][ i ] = gradCol[ j ];

        }

        floatMatrix dddXi_1p, dddXi_1m;

        floatMatrix dddXi_2p, dddXi_2m;

        floatMatrix dddDp, dddDm;

        floatMatrix dddFp, dddFm;

        floatMatrix dddchip, dddchim;

        floatMatrix dddchiNLp, dddchiNLm;

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi, chiNL + deltas, dp,
                                                                         dddXi_1p, dddXi_2p, dddDp, dddFp, dddchip, dddchiNLp ) );

        BOOST_CHECK( !tractionSeparation::computeCurrentDistanceGeneral( Xi_1, Xi_2, D, F, chi, chiNL - deltas, dm,
                                                                         dddXi_1m, dddXi_2m, dddDm, dddFm, dddchim, dddchiNLm ) );

        floatMatrix gradMat = ( dddXi_1p - dddXi_1m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddXi_2p - dddXi_2m ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                d2ddchiNLdXi_2[ j ][ Xi_2.size( ) * i + k ] = gradMat[ j ][ k ];

            }

        }

        gradMat = ( dddDp - dddDm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddFp - dddFm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchip - dddchim ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

        gradMat = ( dddchiNLp - dddchiNLm ) / ( 2 * deltas[ i ] );

        for ( unsigned int j = 0; j < gradMat.size( ); j++ ){

            for ( unsigned int k = 0; k < gradMat[ j ].size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( 0., gradMat[ j ][ k ] ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dddchiNL, dddchiNL_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2ddchiNLdXi_2, d2ddchiNLdXi_2_answer ) );

}

BOOST_AUTO_TEST_CASE( test_decomposeVector ){
    /*!
     * Test the decomposition of a vector into normal and tangential parts
     */

    floatVector d = { 0.69646919, 0.28613933, 0.22685145 };

    floatVector n = { 0.55114872, 0.71925227, 0.42297903 };

    floatVector dn_answer = { 0.37787741, 0.49313221, 0.29000198 };

    floatVector dt_answer = { 0.31859177, -0.20699288, -0.06315053 };

    floatVector dn, dt;

    BOOST_CHECK( !tractionSeparation::decomposeVector( d, n, dn, dt ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dn, dn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dt, dt_answer ) );

    floatMatrix ddndn, ddndd, ddtdn, ddtdd;

    BOOST_CHECK( !tractionSeparation::decomposeVector( d, n, dn, dt, ddndd, ddndn, ddtdd, ddtdn ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dn, dn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dt, dt_answer ) );

    floatVector dn_2, dt_2;

    floatMatrix ddndn_2, ddndd_2, ddtdn_2, ddtdd_2;
    floatMatrix d2dndddd, d2dndddn, d2dndndn;
    floatMatrix d2dtdddd, d2dtdddn, d2dtdndn;

    BOOST_CHECK( !tractionSeparation::decomposeVector( d, n, dn_2, dt_2, ddndd_2, ddndn_2, ddtdd_2, ddtdn_2,
                                                       d2dndddd, d2dndddn, d2dndndn, d2dtdddd, d2dtdddn, d2dtdndn ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dn_2, dn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dt_2, dt_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ddndd, ddndd_2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ddndn, ddndn_2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ddtdd, ddtdd_2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ddtdn, ddtdn_2 ) );

    floatMatrix ddndd_answer( dn.size( ), floatVector( d.size( ), 0 ) );
    floatMatrix ddndn_answer( dn.size( ), floatVector( n.size( ), 0 ) );
    floatMatrix ddtdd_answer( dt.size( ), floatVector( d.size( ), 0 ) );
    floatMatrix ddtdn_answer( dt.size( ), floatVector( n.size( ), 0 ) );

    floatMatrix d2dndddd_answer( dn.size( ), floatVector( d.size( ) * d.size( ), 0 ) );
    floatMatrix d2dndddn_answer( dn.size( ), floatVector( d.size( ) * n.size( ), 0 ) );
    floatMatrix d2dndndn_answer( dn.size( ), floatVector( n.size( ) * n.size( ), 0 ) );
    floatMatrix d2dtdddd_answer( dt.size( ), floatVector( d.size( ) * d.size( ), 0 ) );
    floatMatrix d2dtdddn_answer( dt.size( ), floatVector( d.size( ) * n.size( ), 0 ) );
    floatMatrix d2dtdndn_answer( dt.size( ), floatVector( n.size( ) * n.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < d.size( ); i++ ){

        floatVector delta( d.size( ), 0 );

        delta[ i ] += eps * std::abs( d[ i ] ) + eps;

        floatVector dnp, dnm;

        floatVector dtp, dtm;

        BOOST_CHECK( !tractionSeparation::decomposeVector( d + delta, n, dnp, dtp ) );

        BOOST_CHECK( !tractionSeparation::decomposeVector( d - delta, n, dnm, dtm ) );

        for ( unsigned int j = 0; j < d.size( ); j++ ){

            ddndd_answer[ j ][ i ] = ( dnp[ j ] - dnm[ j ] ) / ( 2 * delta[ i ] );

            ddtdd_answer[ j ][ i ] = ( dtp[ j ] - dtm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix ddnddp, ddnddm;
        floatMatrix ddndnp, ddndnm;
        floatMatrix ddtddp, ddtddm;
        floatMatrix ddtdnp, ddtdnm;

        BOOST_CHECK( !tractionSeparation::decomposeVector( d + delta, n, dnp, dtp, ddnddp, ddndnp, ddtddp, ddtdnp ) );

        BOOST_CHECK( !tractionSeparation::decomposeVector( d - delta, n, dnm, dtm, ddnddm, ddndnm, ddtddm, ddtdnm ) );

        for ( unsigned int j = 0; j < d.size( ); j++ ){

            for ( unsigned int k = 0; k < d.size( ); k++ ){

                d2dndddd_answer[ j ][ d.size( ) * k + i ] += ( ddnddp[ j ][ k ] - ddnddm[ j ][ k ] ) / ( 2 * delta[ i ] );

                d2dtdddd_answer[ j ][ d.size( ) * k + i ] += ( ddtddp[ j ][ k ] - ddtddm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( ddndd, ddndd_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ddtdd, ddtdd_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dtdddd, d2dtdddd_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dndddd, d2dndddd_answer ) );

    for ( unsigned int i = 0; i < n.size( ); i++ ){

        floatVector delta( n.size( ), 0 );

        delta[ i ] += eps * std::abs( n[ i ] ) + eps;

        floatVector dnp, dnm;

        floatVector dtp, dtm;

        BOOST_CHECK( !tractionSeparation::decomposeVector( d, n + delta, dnp, dtp ) );

        BOOST_CHECK( !tractionSeparation::decomposeVector( d, n - delta, dnm, dtm ) );

        for ( unsigned int j = 0; j < d.size( ); j++ ){

            ddndn_answer[ j ][ i ] = ( dnp[ j ] - dnm[ j ] ) / ( 2 * delta[ i ] );

            ddtdn_answer[ j ][ i ] = ( dtp[ j ] - dtm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix ddnddp, ddnddm;
        floatMatrix ddndnp, ddndnm;
        floatMatrix ddtddp, ddtddm;
        floatMatrix ddtdnp, ddtdnm;

        BOOST_CHECK( !tractionSeparation::decomposeVector( d, n + delta, dnp, dtp, ddnddp, ddndnp, ddtddp, ddtdnp ) );

        BOOST_CHECK( !tractionSeparation::decomposeVector( d, n - delta, dnm, dtm, ddnddm, ddndnm, ddtddm, ddtdnm ) );

        for ( unsigned int j = 0; j < d.size( ); j++ ){

            for ( unsigned int k = 0; k < d.size( ); k++ ){

                d2dndddn_answer[ j ][ d.size( ) * k + i ] += ( ddnddp[ j ][ k ] - ddnddm[ j ][ k ] ) / ( 2 * delta[ i ] );

                d2dtdddn_answer[ j ][ d.size( ) * k + i ] += ( ddtddp[ j ][ k ] - ddtddm[ j ][ k ] ) / ( 2 * delta[ i ] );

                d2dndndn_answer[ j ][ d.size( ) * k + i ] += ( ddndnp[ j ][ k ] - ddndnm[ j ][ k ] ) / ( 2 * delta[ i ] );

                d2dtdndn_answer[ j ][ d.size( ) * k + i ] += ( ddtdnp[ j ][ k ] - ddtdnm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( ddndn, ddndn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ddtdn, ddtdn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dtdddn, d2dtdddn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dndddn, d2dndddn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dtdndn, d2dtdndn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dndndn, d2dndndn_answer ) );

}

BOOST_AUTO_TEST_CASE( test_computeLinearTractionEnergy ){
    /*!
     * Compute the linear traction energy
     */

    floatVector dt = { 0.69646919, 0.28613933, 0.22685145 };

    floatVector dn = { 0.55131477, 0.71946897, 0.42310646 };

    floatVector parameters = { 0.9807641983846155, 0.6848297385848633 };

    floatType energy_answer = 0.702429252330725;

    floatType energy;

    BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters, energy ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy_answer, energy ) );

    floatVector denergydn, denergydt;

    energy = 0.;

    BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters, energy, denergydn, denergydt ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy_answer, energy ) );

    floatType energy_2;

    floatVector denergydn_2, denergydt_2;

    floatVector d2energydndn, d2energydndt, d2energydtdt;

    BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters, energy_2, denergydn_2, denergydt_2,
                                                                   d2energydndn, d2energydndt, d2energydtdt ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy_answer, energy_2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydn_2, denergydn ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydt_2, denergydt ) );

    floatType energy_3;

    floatVector denergydn_3, denergydt_3;

    floatVector denergydParameters;

    BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters, energy_3, denergydn_3, denergydt_3, denergydParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy_answer, energy_3 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydn_3, denergydn ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydt_3, denergydt ) );

    floatType energy_4;

    floatVector denergydn_4, denergydt_4, denergydParameters_4;

    floatVector d2energydndn_4, d2energydndt_4, d2energyddndParameters;
    floatVector d2energydtdt_4, d2energyddtdParameters;
    floatVector d2energydParametersdParameters;

    BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters, energy_4, denergydn_4, denergydt_4, denergydParameters_4,
                                                                   d2energydndn_4, d2energydndt_4, d2energyddndParameters,
                                                                   d2energydtdt_4, d2energyddtdParameters,
                                                                   d2energydParametersdParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy_answer, energy_4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydn_4, denergydn ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydt_4, denergydt ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydParameters_4, denergydParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydndn_4, d2energydndn ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydndt_4, d2energydndt ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydtdt_4, d2energydtdt ) );

    floatType eps = 1e-6;

    floatVector denergydn_answer( dn.size( ), 0 );

    floatVector denergydt_answer( dt.size( ), 0 );

    floatVector denergydParameters_answer( parameters.size( ), 0 );

    floatVector d2energydndn_answer( dn.size( ) * dn.size( ), 0 );

    floatVector d2energydndt_answer( dn.size( ) * dt.size( ), 0 );

    floatVector d2energydtdt_answer( dt.size( ) * dt.size( ), 0 );

    floatVector d2energyddndParameters_answer( dn.size( ) * parameters.size( ), 0 );

    floatVector d2energyddtdParameters_answer( dn.size( ) * parameters.size( ), 0 );

    floatVector d2energydParametersdParameters_answer( parameters.size( ) * parameters.size( ), 0 );

    for ( unsigned int i = 0; i < dn.size( ); i++ ){

        floatVector delta( dn.size( ), 0 );

        delta[ i ] = eps * std::abs( dn[ i ] ) + eps;

        floatType energyp, energym;

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn + delta, dt, parameters, energyp ) );

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn - delta, dt, parameters, energym ) );

        denergydn_answer[ i ] += ( energyp - energym ) / ( 2 * delta[ i ] );

        floatVector dednp, dednm;

        floatVector dedtp, dedtm;

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn + delta, dt, parameters, energyp, dednp, dedtp ) );

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn - delta, dt, parameters, energym, dednm, dedtm ) );

        for ( unsigned int j = 0; j < dn.size( ); j++ ){

            d2energydndn_answer[ dn.size( ) * j + i ] += ( dednp[ j ] - dednm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydn, denergydn_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydndn, d2energydndn_answer ) );

    for ( unsigned int i = 0; i < dt.size( ); i++ ){

        floatVector delta( dt.size( ), 0 );

        delta[ i ] = eps * std::abs( dt[ i ] ) + eps;

        floatType energyp, energym;

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt + delta, parameters, energyp ) );

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt - delta, parameters, energym ) );

        denergydt_answer[ i ] += ( energyp - energym ) / ( 2 * delta[ i ] );

        floatVector dednp, dednm;

        floatVector dedtp, dedtm;

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt + delta, parameters, energyp, dednp, dedtp ) );

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt - delta, parameters, energym, dednm, dedtm ) );

        for ( unsigned int j = 0; j < dt.size( ); j++ ){

            d2energydndt_answer[ dt.size( ) * j + i ] += ( dednp[ j ] - dednm[ j ] ) / ( 2 * delta[ i ] );

            d2energydtdt_answer[ dt.size( ) * j + i ] += ( dedtp[ j ] - dedtm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydt, denergydt_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydndt, d2energydndt_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydtdt, d2energydtdt_answer ) );

    for ( unsigned int i = 0; i < parameters.size( ); i++ ){

        floatVector delta( parameters.size( ), 0 );

        delta[ i ] = eps * std::abs( parameters[ i ] ) + eps;

        floatType energyp, energym;

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters + delta, energyp ) );

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters - delta, energym ) );

        denergydParameters_answer[ i ] += ( energyp - energym ) / ( 2 * delta[ i ] );

        floatVector dednp, dednm;

        floatVector dedtp, dedtm;

        floatVector dedPp, dedPm;

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters + delta, energyp, dednp, dedtp, dedPp ) );

        BOOST_CHECK( !tractionSeparation::computeLinearTractionEnergy( dn, dt, parameters - delta, energym, dednm, dedtm, dedPm ) );

        for ( unsigned int j = 0; j < dt.size( ); j++ ){

            d2energyddndParameters_answer[ parameters.size( ) * j + i ] += ( dednp[ j ] - dednm[ j ] ) / ( 2 * delta[ i ] );

            d2energyddtdParameters_answer[ parameters.size( ) * j + i ] += ( dedtp[ j ] - dedtm[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < parameters.size( ); j++ ){

            d2energydParametersdParameters_answer[ parameters.size( ) * j + i ] += ( dedPp[ j ] - dedPm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydParameters, denergydParameters_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energyddndParameters, d2energyddndParameters_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energyddtdParameters, d2energyddtdParameters_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydParametersdParameters, d2energydParametersdParameters_answer ) );

}

BOOST_AUTO_TEST_CASE( test_computeNansonsRelation ){
    /*!
     * Compute Nanson's relation
     */

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector dAN = { 0.02650791,  0.03853289, -0.05628004 };

    floatVector dan_answer = { 0.01713406,  0.04519858, -0.00390936 };

    floatVector dan;

    BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient, dAN, dan ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dan, dan_answer ) );

    floatMatrix ddandF, ddanddAN;

    BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient, dAN, dan, ddandF, ddanddAN ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dan, dan_answer ) );

    floatVector dan_2;

    floatMatrix ddandF_2, ddanddAN_2;

    floatMatrix d2dandFdF, d2dandFddAN;

    BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient, dAN, dan_2, ddandF_2, ddanddAN_2, d2dandFdF, d2dandFddAN ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dan_2, dan_answer ) );

    floatMatrix ddandF_answer( dan_answer.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix ddanddAN_answer( dan_answer.size( ), floatVector( dAN.size( ), 0 ) );

    floatMatrix d2dandFdF_answer( dan_answer.size( ), floatVector( deformationGradient.size( ) * deformationGradient.size( ), 0 ) );

    floatMatrix d2dandFddAN_answer( dan_answer.size( ), floatVector( deformationGradient.size( ) * dAN.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::abs( deformationGradient[ i ] ) + eps;

        floatVector danp, danm;

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient + delta, dAN, danp ) );

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient - delta, dAN, danm ) );

        for ( unsigned int j = 0; j < dan_answer.size( ); j++ ){

            ddandF_answer[ j ][ i ] += ( danp[ j ] - danm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix ddandFp, ddandFm;

        floatMatrix ddanddANp, ddanddANm;

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient + delta, dAN, danp, ddandFp, ddanddANp ) );

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient - delta, dAN, danm, ddandFm, ddanddANm ) );

        for ( unsigned int j = 0; j < dan_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < dan_answer.size( ); k++ ){

                for ( unsigned int K = 0; K < dan_answer.size( ); K++ ){

                    d2dandFdF_answer[ j ][ dan_answer.size( ) * dan_answer.size( ) * dan_answer.size( ) * k + dan_answer.size( ) * dan_answer.size( ) * K + i ] += ( ddandFp[ j ][ dan_answer.size( ) * k + K ] - ddandFm[ j ][ dan_answer.size( ) * k + K ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( ddandF, ddandF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dandFdF, d2dandFdF_answer ) );

    for ( unsigned int i = 0; i < dAN.size( ); i++ ){

        floatVector delta( dAN.size( ), 0 );

        delta[ i ] += eps * std::abs( dAN[ i ] ) + eps;

        floatVector danp, danm;

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient, dAN + delta, danp ) );

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient, dAN - delta, danm ) );

        for ( unsigned int j = 0; j < dan_answer.size( ); j++ ){

            ddanddAN_answer[ j ][ i ] += ( danp[ j ] - danm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix ddandFp, ddandFm;

        floatMatrix ddanddANp, ddanddANm;

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient, dAN + delta, danp, ddandFp, ddanddANp ) );

        BOOST_CHECK( !tractionSeparation::computeNansonsRelation( deformationGradient, dAN - delta, danm, ddandFm, ddanddANm ) );

        for ( unsigned int j = 0; j < dan_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < dan_answer.size( ); k++ ){

                for ( unsigned int K = 0; K < dan_answer.size( ); K++ ){

                    d2dandFddAN_answer[ j ][ dan_answer.size( ) * dan_answer.size( ) * k + dan_answer.size( ) * K + i ] += ( ddandFp[ j ][ dan_answer.size( ) * k + K ] - ddandFm[ j ][ dan_answer.size( ) * k + K ] ) / ( 2 * delta[ i ] );

                }

            }

        }

        BOOST_CHECK( vectorTools::fuzzyEquals( ( ddanddANp - ddanddANm ) / ( 2 * delta[ i ] ), floatMatrix( dan_answer.size( ), floatVector( dAN.size( ), 0 ) ) ) );

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( ddanddAN, ddanddAN_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2dandFddAN, d2dandFddAN_answer ) );

}

BOOST_AUTO_TEST_CASE( test_computeOverlapDistanceLagrangian ){

    floatVector X = { 0.39293837, -0.42772133, -0.54629709,  0.10262954 };

    floatVector chi_nl = { 0.43893794, -0.15378708,  0.9615284 ,
                           0.36965948, -0.0381362 , -0.21576496,
                          -0.31364397,  0.45809941, -0.12285551 };

    floatVector xi_t = { -0.88064421, -0.20391149,  0.47599081 };

    floatType R_nl = 0.7453818629935041;

    floatType L_answer = 0.5498070003214279;

    floatType L;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl, L ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L, L_answer ) );

    floatVector dLdX, dLdchi_nl, dLdxi_t;
    
    floatType dLdR_nl;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl, L, dLdX, dLdchi_nl, dLdxi_t, dLdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L, L_answer ) );

    floatType L_2;

    floatVector dLdX_2, dLdchi_nl_2, dLdxi_t_2;

    floatType dLdR_nl_2;

    floatVector d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t, d2LdXdR_nl, d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t, d2Ldchi_nldR_nl, d2Ldxi_tdxi_t, d2Ldxi_tdR_nl;

    floatType d2LdR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl, L_2, dLdX_2, dLdchi_nl_2, dLdxi_t_2, dLdR_nl_2,
                                                                        d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t, d2LdXdR_nl,
                                                                        d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t, d2Ldchi_nldR_nl,
                                                                        d2Ldxi_tdxi_t, d2Ldxi_tdR_nl,
                                                                        d2LdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L_2, L_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdX_2, dLdX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdchi_nl_2, dLdchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdxi_t_2, dLdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdR_nl_2, dLdR_nl ) );

    floatType L_3;

    floatVector dLdX_3, dLdchi_nl_3, dLdxi_t_3;

    floatType dLdR_nl_3;

    floatVector d2LdXdX_3, d2LdXdchi_nl_3, d2LdXdxi_t_3, d2LdXdR_nl_3, d2Ldchi_nldchi_nl_3, d2Ldchi_nldxi_t_3, d2Ldchi_nldR_nl_3, d2Ldxi_tdxi_t_3, d2Ldxi_tdR_nl_3;

    floatType d2LdR_nldR_nl_3;

    floatVector d3LdXdXdX, d3LdXdXdchi_nl, d3LdXdchi_nldchi_nl, d3LdXdXdxi_t, d3LdXdchi_nldxi_t, d3LdXdXdR_nl, d3LdXdchi_nldR_nl, d3LdXdxi_tdxi_t, d3LdXdxi_tdR_nl, d3LdXdR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl, L_3, dLdX_3, dLdchi_nl_3, dLdxi_t_3, dLdR_nl_3,
                                                                        d2LdXdX_3, d2LdXdchi_nl_3, d2LdXdxi_t_3, d2LdXdR_nl_3,
                                                                        d2Ldchi_nldchi_nl_3, d2Ldchi_nldxi_t_3, d2Ldchi_nldR_nl_3,
                                                                        d2Ldxi_tdxi_t_3, d2Ldxi_tdR_nl_3,
                                                                        d2LdR_nldR_nl_3,
                                                                        d3LdXdXdX, d3LdXdXdchi_nl, d3LdXdchi_nldchi_nl,
                                                                        d3LdXdXdxi_t, d3LdXdchi_nldxi_t,
                                                                        d3LdXdXdR_nl, d3LdXdchi_nldR_nl,
                                                                        d3LdXdxi_tdxi_t, d3LdXdxi_tdR_nl,
                                                                        d3LdXdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L_3, L_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdX_3, dLdX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdchi_nl_3, dLdchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdxi_t_3, dLdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdR_nl_3, dLdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdX_3, d2LdXdX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdchi_nl_3, d2LdXdchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdxi_t_3, d2LdXdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdR_nl_3, d2LdXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldchi_nl_3, d2Ldchi_nldchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldxi_t_3, d2Ldchi_nldxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldR_nl_3, d2Ldchi_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldxi_tdxi_t_3, d2Ldxi_tdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldxi_tdR_nl_3, d2Ldxi_tdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdR_nldR_nl_3, d2LdR_nldR_nl ) );

    floatType L_4;

    floatVector dLdX_4, dLdchi_nl_4, dLdxi_t_4;

    floatType dLdR_nl_4;

    floatVector d2LdXdX_4, d2LdXdchi_nl_4, d2LdXdxi_t_4, d2LdXdR_nl_4, d2Ldchi_nldchi_nl_4, d2Ldchi_nldxi_t_4, d2Ldchi_nldR_nl_4, d2Ldxi_tdxi_t_4, d2Ldxi_tdR_nl_4;

    floatType d2LdR_nldR_nl_4;

    floatVector d3LdXdXdX_4, d3LdXdXdchi_nl_4, d3LdXdchi_nldchi_nl_4, d3LdXdXdxi_t_4, d3LdXdchi_nldxi_t_4, d3LdXdXdR_nl_4, d3LdXdchi_nldR_nl_4, d3LdXdxi_tdxi_t_4,
                d3LdXdxi_tdR_nl_4, d3LdXdR_nldR_nl_4;

    floatVector d4LdXdXdchi_nldchi_nl;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl, L_4, dLdX_4, dLdchi_nl_4, dLdxi_t_4, dLdR_nl_4,
                                                                        d2LdXdX_4, d2LdXdchi_nl_4, d2LdXdxi_t_4, d2LdXdR_nl_4,
                                                                        d2Ldchi_nldchi_nl_4, d2Ldchi_nldxi_t_4, d2Ldchi_nldR_nl_4,
                                                                        d2Ldxi_tdxi_t_4, d2Ldxi_tdR_nl_4,
                                                                        d2LdR_nldR_nl_4,
                                                                        d3LdXdXdX_4, d3LdXdXdchi_nl_4, d3LdXdchi_nldchi_nl_4,
                                                                        d3LdXdXdxi_t_4, d3LdXdchi_nldxi_t_4,
                                                                        d3LdXdXdR_nl_4, d3LdXdchi_nldR_nl_4,
                                                                        d3LdXdxi_tdxi_t_4, d3LdXdxi_tdR_nl_4,
                                                                        d3LdXdR_nldR_nl_4,
                                                                        d4LdXdXdchi_nldchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L_4, L_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdX_4, dLdX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdchi_nl_4, dLdchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdxi_t_4, dLdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdR_nl_4, dLdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdX_4, d2LdXdX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdchi_nl_4, d2LdXdchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdxi_t_4, d2LdXdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdR_nl_4, d2LdXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldchi_nl_4, d2Ldchi_nldchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldxi_t_4, d2Ldchi_nldxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldR_nl_4, d2Ldchi_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldxi_tdxi_t_4, d2Ldxi_tdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldxi_tdR_nl_4, d2Ldxi_tdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdR_nldR_nl_4, d2LdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdX_4, d3LdXdXdX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdchi_nl_4, d3LdXdXdchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdchi_nldchi_nl_4, d3LdXdchi_nldchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdxi_t_4, d3LdXdXdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdchi_nldxi_t_4, d3LdXdchi_nldxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdR_nl_4, d3LdXdXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdchi_nldR_nl_4, d3LdXdchi_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdxi_tdxi_t_4, d3LdXdxi_tdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdxi_tdR_nl_4, d3LdXdxi_tdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdR_nldR_nl_4, d3LdXdR_nldR_nl ) );

    floatType eps = 1e-6;

    floatVector dLdX_answer( X.size( ), 0 );

    floatVector dLdchi_nl_answer( chi_nl.size( ), 0 );

    floatVector dLdxi_t_answer( xi_t.size( ), 0 );

    floatType dLdR_nl_answer = 0;

    floatVector d2LdXdX_answer( X.size( ) * X.size( ), 0 );

    floatVector d2LdXdchi_nl_answer( X.size( ) * chi_nl.size( ), 0 );

    floatVector d2LdXdxi_t_answer( X.size( ) * xi_t.size( ), 0 );

    floatVector d2LdXdR_nl_answer( X.size( ), 0 );

    floatVector d2Ldchi_nldchi_nl_answer( chi_nl.size( ) * chi_nl.size( ), 0 );

    floatVector d2Ldchi_nldxi_t_answer( chi_nl.size( ) * xi_t.size( ), 0 );

    floatVector d2Ldchi_nldR_nl_answer( chi_nl.size( ), 0 );

    floatVector d2Ldxi_tdxi_t_answer( xi_t.size( ) * xi_t.size( ), 0 );

    floatVector d2Ldxi_tdR_nl_answer( xi_t.size( ), 0 );

    floatType d2LdR_nldR_nl_answer = 0;

    floatVector d3LdXdXdX_answer( X.size( ) * X.size( ) * X.size( ), 0 );

    floatVector d3LdXdXdchi_nl_answer( X.size( ) * X.size( ) * chi_nl.size( ), 0 );

    floatVector d3LdXdXdxi_t_answer( X.size( ) * X.size( ) * xi_t.size( ), 0 );

    floatVector d3LdXdXdR_nl_answer( X.size( ) * X.size( ), 0 );

    floatVector d3LdXdchi_nldchi_nl_answer( X.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 );

    floatVector d3LdXdchi_nldxi_t_answer( X.size( ) * chi_nl.size( ) * xi_t.size( ), 0 );

    floatVector d3LdXdchi_nldR_nl_answer( X.size( ) * chi_nl.size( ), 0 );

    floatVector d3LdXdxi_tdxi_t_answer( X.size( ) * xi_t.size( ) * xi_t.size( ), 0 );

    floatVector d3LdXdxi_tdR_nl_answer( X.size( ) * xi_t.size( ), 0 );

    floatVector d3LdXdR_nldR_nl_answer( X.size( ), 0 );

    floatVector d4LdXdXdXdX_answer( X.size( ) * X.size( ) * X.size( ) * X.size( ), 0 );

    floatVector d4LdXdXdXdchi_nl_answer( X.size( ) * X.size( ) * X.size( ) * chi_nl.size( ), 0 );

    floatVector d4LdXdXdchi_nldchi_nl_answer( X.size( ) * X.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 );

    floatVector d4LdXdchi_nldchi_nldchi_nl_answer( X.size( ) * chi_nl.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 );

    floatVector d4LdXdXdXdxi_t_answer( X.size( ) * X.size( ) * X.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdXdchi_nldxi_t_answer( X.size( ) * X.size( ) * chi_nl.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdchi_nldchi_nldxi_t_answer( X.size( ) * chi_nl.size( ) * chi_nl.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdXdxi_tdxi_t_answer( X.size( ) * X.size( ) * xi_t.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdchi_nldxi_tdxi_t_answer( X.size( ) * chi_nl.size( ) * xi_t.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdxi_tdxi_tdxi_t_answer( X.size( ) * xi_t.size( ) * xi_t.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdXdXdR_nl_answer( X.size( ) * X.size( ) * X.size( ), 0 );

    floatVector d4LdXdXdchi_nldR_nl_answer( X.size( ) * X.size( ) * chi_nl.size( ), 0 );

    floatVector d4LdXdXdR_nldR_nl_answer( X.size( ) * X.size( ), 0 );

    floatVector d4LdXdchi_nldchi_nldR_nl_answer( X.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 );

    floatVector d4LdXdchi_nldR_nldR_nl_answer( X.size( ) * chi_nl.size( ), 0 );

    floatVector d4LdXdR_nldR_nldR_nl_answer( X.size( ), 0 );

    floatVector d4LdXdXdxi_tdR_nl_answer( X.size( ) * X.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdchi_nldxi_tdR_nl_answer( X.size( ) * chi_nl.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdxi_tdxi_tdR_nl_answer( X.size( ) * xi_t.size( ) * xi_t.size( ), 0 );

    floatVector d4LdXdxi_tdR_nldR_nl_answer( X.size( ) * xi_t.size( ), 0 );

    for ( unsigned int i = 0; i < X.size( ); i++ ){

        floatVector delta( X.size( ), 0 );

        delta[ i ] += eps * std::abs( X[ i ] ) + eps;

        floatType Lp, Lm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X + delta, chi_nl, xi_t, R_nl, Lp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X - delta, chi_nl, xi_t, R_nl, Lm ) );

        dLdX_answer[ i ] += ( Lp - Lm ) / ( 2 * delta[ i ] );

        floatVector dLdXp, dLdXm;

        floatVector dLdchi_nlp, dLdchi_nlm;

        floatVector dLdxi_tp, dLdxi_tm;

        floatType dLdR_nlp, dLdR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X + delta, chi_nl, xi_t, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X - delta, chi_nl, xi_t, R_nl, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            d2LdXdX_answer[ X.size( ) * j + i ] = ( dLdXp[ j ] - dLdXm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatVector d2LdXdXp, d2LdXdXm;

        floatVector d2LdXdchi_nlp, d2LdXdchi_nlm;

        floatVector d2LdXdxi_tp, d2LdXdxi_tm;

        floatVector d2LdXdR_nlp, d2LdXdR_nlm;

        floatVector d2Ldchi_nldchi_nlp, d2Ldchi_nldchi_nlm;

        floatVector d2Ldchi_nldxi_tp, d2Ldchi_nldxi_tm;

        floatVector d2Ldchi_nldR_nlp, d2Ldchi_nldR_nlm;

        floatVector d2Ldxi_tdxi_tp, d2Ldxi_tdxi_tm;

        floatVector d2Ldxi_tdR_nlp, d2Ldxi_tdR_nlm;

        floatType d2LdR_nldR_nlp, d2LdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X + delta, chi_nl, xi_t, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X - delta, chi_nl, xi_t, R_nl, Lp, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                d3LdXdXdX_answer[ X.size( ) * X.size( ) * j + X.size( ) * k + i ] += ( d2LdXdXp[ X.size( ) * j + k ] - d2LdXdXm[ X.size( ) * j + k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatVector d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp, d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp, d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp, d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp, d3LdXdR_nldR_nlp;
        floatVector d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm, d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm, d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm, d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm, d3LdXdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X + delta, chi_nl, xi_t, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp,
                                                                            d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp,
                                                                            d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp,
                                                                            d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp,
                                                                            d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp,
                                                                            d3LdXdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X - delta, chi_nl, xi_t, R_nl, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm,
                                                                            d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm,
                                                                            d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm,
                                                                            d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm,
                                                                            d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm,
                                                                            d3LdXdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                for ( unsigned int l = 0; l < X.size( ); l++ ){

                    d4LdXdXdXdX_answer[ X.size( ) * X.size( ) * X.size( ) * j + X.size( ) * X.size( ) * k + X.size( ) * l + i ]
                        = ( d3LdXdXdXp[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] - d3LdXdXdXm[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdX, dLdX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdX, d2LdXdX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdX, d3LdXdXdX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdXdX_answer, floatVector( X.size( ) * X.size( ) * X.size( ) * X.size( ), 0 ) ) );

    for ( unsigned int i = 0; i < chi_nl.size( ); i++ ){

        floatVector delta( chi_nl.size( ), 0 );

        delta[ i ] += eps * std::abs( chi_nl[ i ] ) + eps;

        floatType Lp, Lm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl + delta, xi_t, R_nl, Lp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl - delta, xi_t, R_nl, Lm ) );

        dLdchi_nl_answer[ i ] += ( Lp - Lm ) / ( 2 * delta[ i ] );

        floatVector dLdXp, dLdXm;

        floatVector dLdchi_nlp, dLdchi_nlm;

        floatVector dLdxi_tp, dLdxi_tm;

        floatType dLdR_nlp, dLdR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl + delta, xi_t, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl - delta, xi_t, R_nl, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            d2LdXdchi_nl_answer[ chi_nl.size( ) * j + i ] = ( dLdXp[ j ] - dLdXm[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < chi_nl.size( ); j++ ){

            d2Ldchi_nldchi_nl_answer[ chi_nl.size( ) * j + i ] = ( dLdchi_nlp[ j ] - dLdchi_nlm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatVector d2LdXdXp, d2LdXdXm;

        floatVector d2LdXdchi_nlp, d2LdXdchi_nlm;

        floatVector d2LdXdxi_tp, d2LdXdxi_tm;

        floatVector d2LdXdR_nlp, d2LdXdR_nlm;

        floatVector d2Ldchi_nldchi_nlp, d2Ldchi_nldchi_nlm;

        floatVector d2Ldchi_nldxi_tp, d2Ldchi_nldxi_tm;

        floatVector d2Ldchi_nldR_nlp, d2Ldchi_nldR_nlm;

        floatVector d2Ldxi_tdxi_tp, d2Ldxi_tdxi_tm;

        floatVector d2Ldxi_tdR_nlp, d2Ldxi_tdR_nlm;

        floatType d2LdR_nldR_nlp, d2LdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl + delta, xi_t, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl - delta, xi_t, R_nl, Lp, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                d3LdXdXdchi_nl_answer[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + i ] += ( d2LdXdXp[ X.size( ) * j + k ] - d2LdXdXm[ X.size( ) * j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d3LdXdchi_nldchi_nl_answer[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + i ] += ( d2LdXdchi_nlp[ chi_nl.size( ) * j + k ] - d2LdXdchi_nlm[ chi_nl.size( ) * j + k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatVector d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp, d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp, d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp, d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp, d3LdXdR_nldR_nlp;
        floatVector d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm, d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm, d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm, d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm, d3LdXdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl + delta, xi_t, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp,
                                                                            d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp,
                                                                            d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp,
                                                                            d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp,
                                                                            d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp,
                                                                            d3LdXdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl - delta, xi_t, R_nl, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm,
                                                                            d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm,
                                                                            d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm,
                                                                            d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm,
                                                                            d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm,
                                                                            d3LdXdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                for ( unsigned int l = 0; l < X.size( ); l++ ){

                    d4LdXdXdXdchi_nl_answer[ X.size( ) * X.size( ) * chi_nl.size( ) * j + X.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d3LdXdXdXp[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] - d3LdXdXdXm[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d4LdXdXdchi_nldchi_nl_answer[ X.size( ) * chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d3LdXdXdchi_nlp[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] - d3LdXdXdchi_nlm[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d4LdXdchi_nldchi_nldchi_nl_answer[ chi_nl.size( ) * chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d3LdXdchi_nldchi_nlp[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] - d3LdXdchi_nldchi_nlm[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdchi_nl, dLdchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdchi_nl, d2LdXdchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldchi_nl, d2Ldchi_nldchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdchi_nl, d3LdXdXdchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdchi_nldchi_nl, d3LdXdchi_nldchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdXdchi_nl_answer, floatVector( X.size( ) * X.size( ) * X.size( ) * chi_nl.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdchi_nldchi_nl, d4LdXdXdchi_nldchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdchi_nldchi_nldchi_nl_answer, floatVector( X.size( ) * chi_nl.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) ) );

    for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

        floatVector delta( xi_t.size( ), 0 );

        delta[ i ] += eps * std::abs( xi_t[ i ] ) + eps;

        floatType Lp, Lm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t + delta, R_nl, Lp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t - delta, R_nl, Lm ) );

        dLdxi_t_answer[ i ] += ( Lp - Lm ) / ( 2 * delta[ i ] );

        floatVector dLdXp, dLdXm;

        floatVector dLdchi_nlp, dLdchi_nlm;

        floatVector dLdxi_tp, dLdxi_tm;

        floatType dLdR_nlp, dLdR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t + delta, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t - delta, R_nl, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            d2LdXdxi_t_answer[ xi_t.size( ) * j + i ] = ( dLdXp[ j ] - dLdXm[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < chi_nl.size( ); j++ ){

            d2Ldchi_nldxi_t_answer[ xi_t.size( ) * j + i ] = ( dLdchi_nlp[ j ] - dLdchi_nlm[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < xi_t.size( ); j++ ){

            d2Ldxi_tdxi_t_answer[ xi_t.size( ) * j + i ] = ( dLdxi_tp[ j ] - dLdxi_tm[ j ] ) / ( 2 * delta[ i ] );
        }

        floatVector d2LdXdXp, d2LdXdXm;

        floatVector d2LdXdchi_nlp, d2LdXdchi_nlm;

        floatVector d2LdXdxi_tp, d2LdXdxi_tm;

        floatVector d2LdXdR_nlp, d2LdXdR_nlm;

        floatVector d2Ldchi_nldchi_nlp, d2Ldchi_nldchi_nlm;

        floatVector d2Ldchi_nldxi_tp, d2Ldchi_nldxi_tm;

        floatVector d2Ldchi_nldR_nlp, d2Ldchi_nldR_nlm;

        floatVector d2Ldxi_tdxi_tp, d2Ldxi_tdxi_tm;

        floatVector d2Ldxi_tdR_nlp, d2Ldxi_tdR_nlm;

        floatType d2LdR_nldR_nlp, d2LdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t + delta, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t - delta, R_nl, Lp, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                d3LdXdXdxi_t_answer[ X.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + i ] += ( d2LdXdXp[ X.size( ) * j + k ] - d2LdXdXm[ X.size( ) * j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d3LdXdchi_nldxi_t_answer[ chi_nl.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + i ] += ( d2LdXdchi_nlp[ chi_nl.size( ) * j + k ] - d2LdXdchi_nlm[ chi_nl.size( ) * j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                d3LdXdxi_tdxi_t_answer[ xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + i ] += ( d2LdXdxi_tp[ xi_t.size( ) * j + k ] - d2LdXdxi_tm[ xi_t.size( ) * j + k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatVector d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp, d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp, d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp, d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp, d3LdXdR_nldR_nlp;
        floatVector d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm, d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm, d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm, d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm, d3LdXdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t + delta, R_nl, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp,
                                                                            d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp,
                                                                            d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp,
                                                                            d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp,
                                                                            d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp,
                                                                            d3LdXdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t - delta, R_nl, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm,
                                                                            d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm,
                                                                            d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm,
                                                                            d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm,
                                                                            d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm,
                                                                            d3LdXdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                for ( unsigned int l = 0; l < X.size( ); l++ ){

                    d4LdXdXdXdxi_t_answer[ X.size( ) * X.size( ) * xi_t.size( ) * j + X.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        = ( d3LdXdXdXp[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] - d3LdXdXdXm[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d4LdXdXdchi_nldxi_t_answer[ X.size( ) * chi_nl.size( ) * xi_t.size( ) * j + chi_nl.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        = ( d3LdXdXdchi_nlp[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] - d3LdXdXdchi_nlm[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d4LdXdXdxi_tdxi_t_answer[ X.size( ) * xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        = ( d3LdXdXdxi_tp[ X.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] - d3LdXdXdxi_tm[ X.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d4LdXdchi_nldchi_nldxi_t_answer[ chi_nl.size( ) * chi_nl.size( ) * xi_t.size( ) * j + chi_nl.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        = ( d3LdXdchi_nldchi_nlp[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] - d3LdXdchi_nldchi_nlm[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d4LdXdchi_nldxi_tdxi_t_answer[ chi_nl.size( ) * xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        = ( d3LdXdchi_nldxi_tp[ chi_nl.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] - d3LdXdchi_nldxi_tm[ chi_nl.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d4LdXdxi_tdxi_tdxi_t_answer[ xi_t.size( ) * xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        = ( d3LdXdxi_tdxi_tp[ xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] - d3LdXdxi_tdxi_tm[ xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdxi_t, dLdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdxi_t, d2LdXdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldxi_t, d2Ldchi_nldxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldxi_tdxi_t, d2Ldxi_tdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdxi_t, d3LdXdXdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdchi_nldxi_t, d3LdXdchi_nldxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdxi_tdxi_t, d3LdXdxi_tdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdXdxi_t_answer, floatVector( X.size( ) * X.size( ) * X.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdchi_nldxi_t_answer, floatVector( X.size( ) * X.size( ) * chi_nl.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdxi_tdxi_t_answer, floatVector( X.size( ) * X.size( ) * xi_t.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdchi_nldchi_nldxi_t_answer, floatVector( X.size( ) * chi_nl.size( ) * chi_nl.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdchi_nldxi_tdxi_t_answer, floatVector( X.size( ) * chi_nl.size( ) * xi_t.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdxi_tdxi_tdxi_t_answer, floatVector( X.size( ) * xi_t.size( ) * xi_t.size( ) * xi_t.size( ), 0 ) ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( R_nl ) + eps;

        floatType Lp, Lm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl + delta, Lp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl - delta, Lm ) );

        dLdR_nl_answer += ( Lp - Lm ) / ( 2 * delta );

        floatVector dLdXp, dLdXm;

        floatVector dLdchi_nlp, dLdchi_nlm;

        floatVector dLdxi_tp, dLdxi_tm;

        floatType dLdR_nlp, dLdR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl + delta, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl - delta, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            d2LdXdR_nl_answer[ j + i ] = ( dLdXp[ j ] - dLdXm[ j ] ) / ( 2 * delta );

        }

        for ( unsigned int j = 0; j < chi_nl.size( ); j++ ){

            d2Ldchi_nldR_nl_answer[ j + i ] = ( dLdchi_nlp[ j ] - dLdchi_nlm[ j ] ) / ( 2 * delta );

        }

        for ( unsigned int j = 0; j < xi_t.size( ); j++ ){

            d2Ldxi_tdR_nl_answer[ j + i ] = ( dLdxi_tp[ j ] - dLdxi_tm[ j ] ) / ( 2 * delta );
        }

        d2LdR_nldR_nl_answer = ( dLdR_nlp - dLdR_nlm ) / ( 2 * delta );

        floatVector d2LdXdXp, d2LdXdXm;

        floatVector d2LdXdchi_nlp, d2LdXdchi_nlm;

        floatVector d2LdXdxi_tp, d2LdXdxi_tm;

        floatVector d2LdXdR_nlp, d2LdXdR_nlm;

        floatVector d2Ldchi_nldchi_nlp, d2Ldchi_nldchi_nlm;

        floatVector d2Ldchi_nldxi_tp, d2Ldchi_nldxi_tm;

        floatVector d2Ldchi_nldR_nlp, d2Ldchi_nldR_nlm;

        floatVector d2Ldxi_tdxi_tp, d2Ldxi_tdxi_tm;

        floatVector d2Ldxi_tdR_nlp, d2Ldxi_tdR_nlm;

        floatType d2LdR_nldR_nlp, d2LdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl + delta, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl - delta, Lp, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                d3LdXdXdR_nl_answer[ X.size( ) * j + k + i ] += ( d2LdXdXp[ X.size( ) * j + k ] - d2LdXdXm[ X.size( ) * j + k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d3LdXdchi_nldR_nl_answer[ chi_nl.size( ) * j + k + i ] += ( d2LdXdchi_nlp[ chi_nl.size( ) * j + k ] - d2LdXdchi_nlm[ chi_nl.size( ) * j + k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                d3LdXdxi_tdR_nl_answer[ xi_t.size( ) * j + k + i ] += ( d2LdXdxi_tp[ xi_t.size( ) * j + k ] - d2LdXdxi_tm[ xi_t.size( ) * j + k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d3LdXdR_nldR_nl_answer[ j + k + i ] += ( d2LdXdR_nlp[ j + k ] - d2LdXdR_nlm[ j + k ] ) / ( 2 * delta );

            }

        }


        floatVector d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp, d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp, d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp, d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp, d3LdXdR_nldR_nlp;
        floatVector d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm, d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm, d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm, d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm, d3LdXdR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl + delta, Lp, dLdXp, dLdchi_nlp, dLdxi_tp, dLdR_nlp,
                                                                            d2LdXdXp, d2LdXdchi_nlp, d2LdXdxi_tp, d2LdXdR_nlp,
                                                                            d2Ldchi_nldchi_nlp, d2Ldchi_nldxi_tp, d2Ldchi_nldR_nlp,
                                                                            d2Ldxi_tdxi_tp, d2Ldxi_tdR_nlp,
                                                                            d2LdR_nldR_nlp,
                                                                            d3LdXdXdXp, d3LdXdXdchi_nlp, d3LdXdchi_nldchi_nlp,
                                                                            d3LdXdXdxi_tp, d3LdXdchi_nldxi_tp,
                                                                            d3LdXdXdR_nlp, d3LdXdchi_nldR_nlp,
                                                                            d3LdXdxi_tdxi_tp, d3LdXdxi_tdR_nlp,
                                                                            d3LdXdR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, R_nl - delta, Lm, dLdXm, dLdchi_nlm, dLdxi_tm, dLdR_nlm,
                                                                            d2LdXdXm, d2LdXdchi_nlm, d2LdXdxi_tm, d2LdXdR_nlm,
                                                                            d2Ldchi_nldchi_nlm, d2Ldchi_nldxi_tm, d2Ldchi_nldR_nlm,
                                                                            d2Ldxi_tdxi_tm, d2Ldxi_tdR_nlm,
                                                                            d2LdR_nldR_nlm,
                                                                            d3LdXdXdXm, d3LdXdXdchi_nlm, d3LdXdchi_nldchi_nlm,
                                                                            d3LdXdXdxi_tm, d3LdXdchi_nldxi_tm,
                                                                            d3LdXdXdR_nlm, d3LdXdchi_nldR_nlm,
                                                                            d3LdXdxi_tdxi_tm, d3LdXdxi_tdR_nlm,
                                                                            d3LdXdR_nldR_nlm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            for ( unsigned int k = 0; k < X.size( ); k++ ){

                for ( unsigned int l = 0; l < X.size( ); l++ ){

                    d4LdXdXdXdR_nl_answer[ X.size( ) * X.size( ) * j + X.size( ) * k + l + i ]
                        = ( d3LdXdXdXp[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] - d3LdXdXdXm[ X.size( ) * X.size( ) * j + X.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d4LdXdXdchi_nldR_nl_answer[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l + i ]
                        = ( d3LdXdXdchi_nlp[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] - d3LdXdXdchi_nlm[ X.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d4LdXdXdxi_tdR_nl_answer[ X.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l + i ]
                        = ( d3LdXdXdxi_tp[ X.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] - d3LdXdXdxi_tm[ X.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d4LdXdXdR_nldR_nl_answer[ X.size( ) * j + k + l + i ]
                        = ( d3LdXdXdR_nlp[ X.size( ) * j + k + l ] - d3LdXdXdR_nlm[ X.size( ) * j + k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d4LdXdchi_nldchi_nldR_nl_answer[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l + i ]
                        = ( d3LdXdchi_nldchi_nlp[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] - d3LdXdchi_nldchi_nlm[ chi_nl.size( ) * chi_nl.size( ) * j + chi_nl.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d4LdXdchi_nldxi_tdR_nl_answer[ chi_nl.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l + i ]
                        = ( d3LdXdchi_nldxi_tp[ chi_nl.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] - d3LdXdchi_nldxi_tm[ chi_nl.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d4LdXdchi_nldR_nldR_nl_answer[ chi_nl.size( ) * j + k + l + i ]
                        = ( d3LdXdchi_nldR_nlp[ chi_nl.size( ) * j + k + l ] - d3LdXdchi_nldR_nlm[ chi_nl.size( ) * j + k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d4LdXdxi_tdxi_tdR_nl_answer[ xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l + i ]
                        = ( d3LdXdxi_tdxi_tp[ xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] - d3LdXdxi_tdxi_tm[ xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d4LdXdxi_tdR_nldR_nl_answer[ xi_t.size( ) * j + k + l + i ]
                        = ( d3LdXdxi_tdR_nlp[ xi_t.size( ) * j + k + l ] - d3LdXdxi_tdR_nlm[ xi_t.size( ) * j + k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d4LdXdR_nldR_nldR_nl_answer[ j + k + l + i ] = ( d3LdXdR_nldR_nlp[ j ] - d3LdXdR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdR_nl, dLdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdR_nl, d2LdXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldR_nl, d2Ldchi_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldxi_tdR_nl, d2Ldxi_tdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdR_nldR_nl, d2LdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdXdR_nl, d3LdXdXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdchi_nldR_nl, d3LdXdchi_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdxi_tdR_nl, d3LdXdxi_tdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3LdXdR_nldR_nl, d3LdXdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdXdR_nl_answer, floatVector( X.size( ) * X.size( ) * X.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdchi_nldR_nl_answer, floatVector( X.size( ) * X.size( ) * chi_nl.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdxi_tdR_nl_answer, floatVector( X.size( ) * X.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdXdR_nldR_nl_answer, floatVector( X.size( ) * X.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdchi_nldchi_nldR_nl_answer, floatVector( X.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdchi_nldxi_tdR_nl_answer, floatVector( X.size( ) * chi_nl.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdchi_nldR_nldR_nl_answer, floatVector( X.size( ) * chi_nl.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdxi_tdxi_tdR_nl_answer, floatVector( X.size( ) * xi_t.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdxi_tdR_nldR_nl_answer, floatVector( X.size( ) * xi_t.size( ), 0 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d4LdXdR_nldR_nldR_nl_answer, floatVector( X.size( ), 0 ) ) );

}

BOOST_AUTO_TEST_CASE( test_solveOverlapDistance ){

    floatVector chi_nl = { 2.0, 0.0, 0.0,
                           0.0, 0.5, 0.0,
                           0.0, 0.0, 0.5 };

    floatVector xi_t = { 0.0, 0.125, 0.0 };

    floatType R_nl = 0.5;

    floatVector distance_answer = { 0.0, 0.125, 0.0 };

    floatVector distance;

    BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl, distance ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( distance, distance_answer ) );

    chi_nl = { 1.69646919, 0.28613933, 0.22685145,
               0.55131477, 1.71946897, 0.42310646,
               0.9807642 , 0.68482974, 1.4809319 };

    xi_t = { 0.39211752, 0.34317802, 0.72904971 };

    R_nl = 2.3;

    distance_answer = { -0.881778, -0.787106, 1.57023 };

    BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl, distance ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( distance, distance_answer ) );

    floatMatrix dDistancedchi_nl, dDistancedxi_t;

    floatVector dDistancedR_nl;

    BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl, distance, dDistancedchi_nl, dDistancedxi_t, dDistancedR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( distance, distance_answer ) );

    floatVector distance_2;

    floatMatrix dDistancedchi_nl_2, dDistancedxi_t_2;

    floatVector dDistancedR_nl_2;

    floatMatrix d2Distancedchi_nldchi_nl, d2Distancedchi_nldxi_t, d2Distancedchi_nldR_nl, d2Distancedxi_tdxi_t, d2Distancedxi_tdR_nl;

    floatVector d2DistancedR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl, distance_2, dDistancedchi_nl_2, dDistancedxi_t_2, dDistancedR_nl_2,
                                                            d2Distancedchi_nldchi_nl, d2Distancedchi_nldxi_t, d2Distancedchi_nldR_nl,
                                                            d2Distancedxi_tdxi_t, d2Distancedxi_tdR_nl,
                                                            d2DistancedR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( distance_2, distance_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedchi_nl_2, dDistancedchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedxi_t_2, dDistancedxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedR_nl_2, dDistancedR_nl ) );

    floatMatrix dDistancedchi_nl_answer( distance_answer.size( ), floatVector( chi_nl.size( ), 0 ) );

    floatMatrix dDistancedxi_t_answer( distance_answer.size( ), floatVector( xi_t.size( ), 0 ) );

    floatVector dDistancedR_nl_answer( distance_answer.size( ), 0 );

    floatMatrix d2Distancedchi_nldchi_nl_answer( distance_answer.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d2Distancedchi_nldxi_t_answer( distance_answer.size( ), floatVector( chi_nl.size( ) * xi_t.size( ), 0 ) );

    floatMatrix d2Distancedchi_nldR_nl_answer( distance_answer.size( ), floatVector( chi_nl.size( ), 0 ) );

    floatMatrix d2Distancedxi_tdxi_t_answer( distance_answer.size( ), floatVector( xi_t.size( ) * xi_t.size( ), 0 ) );

    floatMatrix d2Distancedxi_tdR_nl_answer( distance_answer.size( ), floatVector( xi_t.size( ), 0 ) );

    floatVector d2DistancedR_nldR_nl_answer( distance_answer.size( ), 0 );

    floatVector distance_3;

    floatMatrix dDistancedchi_nl_3, dDistancedxi_t_3;

    floatVector dDistancedR_nl_3;

    floatMatrix d2Distancedchi_nldchi_nl_3, d2Distancedchi_nldxi_t_3, d2Distancedchi_nldR_nl_3, d2Distancedxi_tdxi_t_3, d2Distancedxi_tdR_nl_3;

    floatVector d2DistancedR_nldR_nl_3;

    floatMatrix d3Distancedchi_nldchi_nldchi_nl, d3Distancedchi_nldchi_nldxi_t, d3Distancedchi_nldchi_nldR_nl;

    floatMatrix d3Distancedchi_nldxi_tdxi_t, d3Distancedchi_nldxi_tdR_nl;

    floatMatrix d3Distancedchi_nldR_nldR_nl;

    floatMatrix d3Distancedxi_tdxi_tdxi_t, d3Distancedxi_tdxi_tdR_nl;

    floatMatrix d3Distancedxi_tdR_nldR_nl;

    floatVector d3DistancedR_nldR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl, distance_3, dDistancedchi_nl_3, dDistancedxi_t_3, dDistancedR_nl_3,
                                                            d2Distancedchi_nldchi_nl_3, d2Distancedchi_nldxi_t_3, d2Distancedchi_nldR_nl_3,
                                                            d2Distancedxi_tdxi_t_3, d2Distancedxi_tdR_nl_3,
                                                            d2DistancedR_nldR_nl_3,
                                                            d3Distancedchi_nldchi_nldchi_nl, d3Distancedchi_nldchi_nldxi_t, d3Distancedchi_nldchi_nldR_nl,
                                                            d3Distancedchi_nldxi_tdxi_t, d3Distancedchi_nldxi_tdR_nl,
                                                            d3Distancedchi_nldR_nldR_nl,
                                                            d3Distancedxi_tdxi_tdxi_t, d3Distancedxi_tdxi_tdR_nl,
                                                            d3Distancedxi_tdR_nldR_nl,
                                                            d3DistancedR_nldR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( distance_3, distance_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedchi_nl_3, dDistancedchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedxi_t_3, dDistancedxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedR_nl_3, dDistancedR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedchi_nldchi_nl_3, d2Distancedchi_nldchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedchi_nldxi_t_3, d2Distancedchi_nldxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedchi_nldR_nl_3, d2Distancedchi_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedxi_tdxi_t_3, d2Distancedxi_tdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedxi_tdR_nl_3, d2Distancedxi_tdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2DistancedR_nldR_nl_3, d2DistancedR_nldR_nl ) );

    floatMatrix d3Distancedchi_nldchi_nldchi_nl_answer( distance_answer.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3Distancedchi_nldchi_nldxi_t_answer( distance_answer.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ) * xi_t.size( ), 0 ) )
;
    floatMatrix d3Distancedchi_nldchi_nldR_nl_answer( distance_answer.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3Distancedchi_nldxi_tdxi_t_answer( distance_answer.size( ), floatVector( chi_nl.size( ) * xi_t.size( ) * xi_t.size( ), 0 ) );

    floatMatrix d3Distancedchi_nldxi_tdR_nl_answer( distance_answer.size( ), floatVector( chi_nl.size( ) * xi_t.size( ), 0 ) );

    floatMatrix d3Distancedchi_nldR_nldR_nl_answer( distance_answer.size( ), floatVector( chi_nl.size( ), 0 ) );

    floatMatrix d3Distancedxi_tdxi_tdxi_t_answer( distance_answer.size( ), floatVector( xi_t.size( ) * xi_t.size( ) * xi_t.size( ), 0 ) );

    floatMatrix d3Distancedxi_tdxi_tdR_nl_answer( distance_answer.size( ), floatVector( xi_t.size( ) * xi_t.size( ), 0 ) );

    floatMatrix d3Distancedxi_tdR_nldR_nl_answer( distance_answer.size( ), floatVector( xi_t.size( ), 0 ) );

    floatVector d3DistancedR_nldR_nldR_nl_answer( distance_answer.size( ), 0 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < chi_nl.size( ); i++ ){

        floatVector delta( chi_nl.size( ), 0 );

        delta[ i ] += eps * std::abs( chi_nl[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl + delta, xi_t, R_nl, dp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl - delta, xi_t, R_nl, dm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            dDistancedchi_nl_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dDistancedchi_nlp, dDistancedchi_nlm;

        floatMatrix dDistancedxi_tp, dDistancedxi_tm;

        floatVector dDistancedR_nlp, dDistancedR_nlm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl + delta, xi_t, R_nl, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl - delta, xi_t, R_nl, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d2Distancedchi_nldchi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dDistancedchi_nlp[ j ][ k ] - dDistancedchi_nlm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp, d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp;
        floatMatrix d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm, d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm;

        floatVector d2DistancedR_nldR_nlp, d2DistancedR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl + delta, xi_t, R_nl, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp,
                                                                d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp,
                                                                d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp,
                                                                d2DistancedR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl - delta, xi_t, R_nl, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm,
                                                                d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm,
                                                                d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm,
                                                                d2DistancedR_nldR_nlm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3Distancedchi_nldchi_nldchi_nl_answer[ j ][ chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        += ( d2Distancedchi_nldchi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2Distancedchi_nldchi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedchi_nl, dDistancedchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedchi_nldchi_nl, d2Distancedchi_nldchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedchi_nldchi_nldchi_nl, d3Distancedchi_nldchi_nldchi_nl_answer ) );

    for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

        floatVector delta( xi_t.size( ), 0 );

        delta[ i ] += eps * std::abs( xi_t[ i ] ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t + delta, R_nl, dp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t - delta, R_nl, dm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            dDistancedxi_t_answer[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dDistancedchi_nlp, dDistancedchi_nlm;

        floatMatrix dDistancedxi_tp, dDistancedxi_tm;

        floatVector dDistancedR_nlp, dDistancedR_nlm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t + delta, R_nl, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t - delta, R_nl, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d2Distancedchi_nldxi_t_answer[ j ][ xi_t.size( ) * k + i ] = ( dDistancedchi_nlp[ j ][ k ] - dDistancedchi_nlm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                d2Distancedxi_tdxi_t_answer[ j ][ xi_t.size( ) * k + i ] = ( dDistancedxi_tp[ j ][ k ] - dDistancedxi_tm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp, d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp;
        floatMatrix d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm, d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm;

        floatVector d2DistancedR_nldR_nlp, d2DistancedR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t + delta, R_nl, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp,
                                                                d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp,
                                                                d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp,
                                                                d2DistancedR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t - delta, R_nl, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm,
                                                                d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm,
                                                                d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm,
                                                                d2DistancedR_nldR_nlm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3Distancedchi_nldchi_nldxi_t_answer[ j ][ chi_nl.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        += ( d2Distancedchi_nldchi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2Distancedchi_nldchi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d3Distancedchi_nldxi_tdxi_t_answer[ j ][ xi_t.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        += ( d2Distancedchi_nldxi_tp[ j ][ xi_t.size( ) * k + l ] - d2Distancedchi_nldxi_tm[ j ][ xi_t.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d3Distancedxi_tdxi_tdxi_t_answer[ j ][ xi_t.size( ) * xi_t.size( ) * k + xi_t.size( ) * l + i ]
                        += ( d2Distancedxi_tdxi_tp[ j ][ xi_t.size( ) * k + l ] - d2Distancedxi_tdxi_tm[ j ][ xi_t.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedxi_t, dDistancedxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedchi_nldxi_t, d2Distancedchi_nldxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedxi_tdxi_t, d2Distancedxi_tdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedchi_nldchi_nldxi_t, d3Distancedchi_nldchi_nldxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedchi_nldxi_tdxi_t, d3Distancedchi_nldxi_tdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedxi_tdxi_tdxi_t, d3Distancedxi_tdxi_tdxi_t_answer ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( R_nl ) + eps;

        floatVector dp, dm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl + delta, dp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl - delta, dm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            dDistancedR_nl_answer[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

        floatMatrix dDistancedchi_nlp, dDistancedchi_nlm;

        floatMatrix dDistancedxi_tp, dDistancedxi_tm;

        floatVector dDistancedR_nlp, dDistancedR_nlm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl + delta, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl - delta, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d2Distancedchi_nldR_nl_answer[ j ][ k + i ] = ( dDistancedchi_nlp[ j ][ k ] - dDistancedchi_nlm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                d2Distancedxi_tdR_nl_answer[ j ][ k + i ] = ( dDistancedxi_tp[ j ][ k ] - dDistancedxi_tm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2DistancedR_nldR_nl_answer[ j + k + i ] = ( dDistancedR_nlp[ j ] - dDistancedR_nlm[ j ] ) / ( 2 * delta );

            }

        }

        floatMatrix d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp, d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp;
        floatMatrix d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm, d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm;

        floatVector d2DistancedR_nldR_nlp, d2DistancedR_nldR_nlm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl + delta, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp,
                                                                d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp,
                                                                d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp,
                                                                d2DistancedR_nldR_nlp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl - delta, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm,
                                                                d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm,
                                                                d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm,
                                                                d2DistancedR_nldR_nlm ) );

        for ( unsigned int j = 0; j < distance_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3Distancedchi_nldchi_nldR_nl_answer[ j ][ chi_nl.size( ) * k + l + i ]
                        += ( d2Distancedchi_nldchi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2Distancedchi_nldchi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d3Distancedchi_nldxi_tdR_nl_answer[ j ][ xi_t.size( ) * k + l + i ]
                        += ( d2Distancedchi_nldxi_tp[ j ][ xi_t.size( ) * k + l ] - d2Distancedchi_nldxi_tm[ j ][ xi_t.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3Distancedchi_nldR_nldR_nl_answer[ j ][ k + l + i ]
                        += ( d2Distancedchi_nldR_nlp[ j ][ k + l ] - d2Distancedchi_nldR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < xi_t.size( ); k++ ){

                for ( unsigned int l = 0; l < xi_t.size( ); l++ ){

                    d3Distancedxi_tdxi_tdR_nl_answer[ j ][ xi_t.size( ) * k + l + i ]
                        += ( d2Distancedxi_tdxi_tp[ j ][ xi_t.size( ) * k + l ] - d2Distancedxi_tdxi_tm[ j ][ xi_t.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3Distancedxi_tdR_nldR_nl_answer[ j ][ k + l + i ]
                        += ( d2Distancedxi_tdR_nlp[ j ][ k + l ] - d2Distancedxi_tdR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3DistancedR_nldR_nldR_nl_answer[ j ]
                        += ( d2DistancedR_nldR_nlp[ j ] - d2DistancedR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dDistancedR_nl, dDistancedR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedchi_nldR_nl, d2Distancedchi_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Distancedxi_tdR_nl, d2Distancedxi_tdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2DistancedR_nldR_nl, d2DistancedR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedchi_nldchi_nldR_nl, d3Distancedchi_nldchi_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedchi_nldxi_tdR_nl, d3Distancedchi_nldxi_tdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedchi_nldR_nldR_nl, d3Distancedchi_nldR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedxi_tdxi_tdR_nl, d3Distancedxi_tdxi_tdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3Distancedxi_tdR_nldR_nl, d3Distancedxi_tdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3DistancedR_nldR_nldR_nl, d3DistancedR_nldR_nldR_nl_answer ) );

}

BOOST_AUTO_TEST_CASE( test_computeParticleOverlap ){

    floatVector Xi_1 = { 1, 0, 0 };

    floatVector dX = { 2, 0, 0 };

    floatType R_nl = 1;

    floatVector F = { 0.75, 0.0, 0.0,
                      0.00, 1.0, 0.0,
                      0.00, 0.0, 1.0 };

    floatVector chi = { 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0 };

    floatVector gradChi( chi.size( ) * dX.size( ) );

    floatVector overlap_answer = { -0.5, 0.0, 0.0 };

    floatVector overlap;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    Xi_1 = { 0.39293837, -0.42772133, -0.54629709 };

    dX = { 0.10262954,  0.43893794, -0.15378708 };

    R_nl = 1.961528396769231;

    F = { 1.36965948, -0.0381362 , -0.21576496,
         -0.31364397,  1.45809941, -0.12285551,
         -0.88064421, -0.20391149,  1.47599081 };

    chi = { 0.36498346, -0.64909649,  0.06310275,
            0.06365517,  1.26880192,  0.69886359,
            0.44891065,  0.22204702,  1.44488677 };

    gradChi = { -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225,
                -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942,
                -0.37547755, -0.14729739,  0.78677833,  0.88832004,  0.00367335,
                 0.2479059 , -0.76876321, -0.36542904, -0.17034758,  0.73261832,
                -0.49908927, -0.03393147,  0.97111957,  0.03897024,  0.22578905,
                -0.75874267,  0.6526816 };

    overlap_answer = floatVector( Xi_1.size( ), 0 );

    overlap.clear( );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    floatMatrix dOverlapdXi_1, dOverlapddX, dOverlapdF, dOverlapdChi, dOverlapdGradChi;

    floatVector dOverlapdR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap,
                                                              dOverlapdXi_1, dOverlapddX, dOverlapdR_nl, dOverlapdF, dOverlapdChi, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    floatVector overlap_2;

    floatMatrix dOverlapdXi_1_2, dOverlapddX_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdGradChi_2;

    floatVector dOverlapdR_nl_2;

    floatMatrix d2OverlapdXi_1dXi_1, d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF, d2OverlapdXi_1dChi, d2OverlapdXi_1dGradChi,
                d2OverlapddXddX, d2OverlapddXdR_nl, d2OverlapddXdF, d2OverlapddXdChi, d2OverlapddXdGradChi,
                d2OverlapdR_nldF, d2OverlapdR_nldChi, d2OverlapdR_nldGradChi,
                d2OverlapdFdF, d2OverlapdFdChi, d2OverlapdFdGradChi,
                d2OverlapdChidChi, d2OverlapdChidGradChi,
                d2OverlapdGradChidGradChi;

    floatVector d2OverlapdR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap_2,
                                                              dOverlapdXi_1_2, dOverlapddX_2, dOverlapdR_nl_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdGradChi_2,
                                                              d2OverlapdXi_1dXi_1,      d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF,       d2OverlapdXi_1dChi,   d2OverlapdXi_1dGradChi,
                                                              d2OverlapddXddX,          d2OverlapddXdR_nl, d2OverlapddXdF,      d2OverlapddXdChi,       d2OverlapddXdGradChi,
                                                              d2OverlapdR_nldR_nl,      d2OverlapdR_nldF,  d2OverlapdR_nldChi,  d2OverlapdR_nldGradChi,
                                                              d2OverlapdFdF,            d2OverlapdFdChi,   d2OverlapdFdGradChi,
                                                              d2OverlapdChidChi,        d2OverlapdChidGradChi,
                                                              d2OverlapdGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_2, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_2, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_2, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_2, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_2, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_2, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_2, dOverlapdGradChi ) );

    floatMatrix dOverlapdXi_1_answer( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix dOverlapddX_answer( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    floatVector dOverlapdR_nl_answer( overlap_answer.size( ), 0 );

    floatMatrix dOverlapdF_answer( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix dOverlapdChi_answer( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix dOverlapdGradChi_answer( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dXi_1_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1ddX_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dR_nl_answer( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dF_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dChi_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dGradChi_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d2OverlapddXddX_answer( overlap_answer.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d2OverlapddXdR_nl_answer( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    floatMatrix d2OverlapddXdF_answer( overlap_answer.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapddXdChi_answer( overlap_answer.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapddXdGradChi_answer( overlap_answer.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    floatVector d2OverlapdR_nldR_nl_answer( overlap_answer.size( ), 0 );

    floatMatrix d2OverlapdR_nldF_answer( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix d2OverlapdR_nldChi_answer( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix d2OverlapdR_nldGradChi_answer( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdFdF_answer( overlap_answer.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapdFdChi_answer( overlap_answer.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdFdGradChi_answer( overlap_answer.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdChidChi_answer( overlap_answer.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdChidGradChi_answer( overlap_answer.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdGradChidGradChi_answer( overlap_answer.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    floatVector overlap_3;

    floatMatrix dOverlapdXi_1_3, dOverlapddX_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdGradChi_3;

    floatVector dOverlapdR_nl_3;

    floatMatrix d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dGradChi_3,
                d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdGradChi_3,
                d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldGradChi_3,
                d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdGradChi_3,
                d2OverlapdChidChi_3, d2OverlapdChidGradChi_3,
                d2OverlapdGradChidGradChi_3;

    floatVector d2OverlapdR_nldR_nl_3;

    floatMatrix d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dGradChi,
                d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdGradChi,
                d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldGradChi,
                d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdGradChi,
                d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidGradChi,
                d3OverlapdXi_1dGradChidGradChi,
                d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdGradChi,
                d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldGradChi,
                d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdGradChi,
                d3OverlapddXdChidChi, d3OverlapddXdChidGradChi,
                d3OverlapddXdGradChidGradChi,
                d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldGradChi,
                d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdGradChi,
                d3OverlapdR_nldChidChi, d3OverlapdR_nldChidGradChi,
                d3OverlapdR_nldGradChidGradChi,
                d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdGradChi,
                d3OverlapdFdChidChi, d3OverlapdFdChidGradChi,
                d3OverlapdFdGradChidGradChi,
                d3OverlapdChidChidChi, d3OverlapdChidChidGradChi,
                d3OverlapdChidGradChidGradChi,
                d3OverlapdGradChidGradChidGradChi;

    floatVector d3OverlapdR_nldR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap_3,
                                                              dOverlapdXi_1_3, dOverlapddX_3, dOverlapdR_nl_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdGradChi_3,
                                                              d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dGradChi_3,
                                                              d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdGradChi_3,
                                                              d2OverlapdR_nldR_nl_3, d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldGradChi_3,
                                                              d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdGradChi_3,
                                                              d2OverlapdChidChi_3, d2OverlapdChidGradChi_3,
                                                              d2OverlapdGradChidGradChi_3,
                                                              d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dGradChi,
                                                              d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdGradChi,
                                                              d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldGradChi,
                                                              d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdGradChi,
                                                              d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidGradChi,
                                                              d3OverlapdXi_1dGradChidGradChi,
                                                              d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdGradChi,
                                                              d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldGradChi,
                                                              d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdGradChi,
                                                              d3OverlapddXdChidChi, d3OverlapddXdChidGradChi,
                                                              d3OverlapddXdGradChidGradChi,
                                                              d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldGradChi,
                                                              d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdGradChi,
                                                              d3OverlapdR_nldChidChi, d3OverlapdR_nldChidGradChi,
                                                              d3OverlapdR_nldGradChidGradChi,
                                                              d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdGradChi,
                                                              d3OverlapdFdChidChi, d3OverlapdFdChidGradChi,
                                                              d3OverlapdFdGradChidGradChi,
                                                              d3OverlapdChidChidChi, d3OverlapdChidChidGradChi,
                                                              d3OverlapdChidGradChidGradChi,
                                                              d3OverlapdGradChidGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_3, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_3, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_3, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_3, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_3, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_3, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_3, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX_3, d2OverlapddXddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl_3, d2OverlapddXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF_3, d2OverlapddXdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_3, d2OverlapddXdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi_3, d2OverlapddXdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl_3, d2OverlapdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF_3, d2OverlapdR_nldF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_3, d2OverlapdR_nldChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi_3, d2OverlapdR_nldGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF_3, d2OverlapdFdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_3, d2OverlapdFdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi_3, d2OverlapdFdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_3, d2OverlapdChidChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi_3, d2OverlapdChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi_3, d2OverlapdGradChidGradChi ) );

    floatMatrix d3OverlapdXi_1dXi_1dXi_1_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1ddX_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXddX_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChidChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChidGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dGradChidGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXddXddX_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdR_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdF_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdChi_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldR_nl_answer( Xi_1.size( ), floatVector( dX.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldF_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdF_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdChi_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdChidChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdChidGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdGradChidGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatVector d3OverlapdR_nldR_nldR_nl_answer( Xi_1.size( ), 0 );

    floatMatrix d3OverlapdR_nldR_nldF_answer( Xi_1.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldR_nldChi_answer( Xi_1.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldR_nldGradChi_answer( Xi_1.size( ), floatVector( gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdF_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChidChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChidGradChi_answer( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldGradChidGradChi_answer( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdF_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdChi_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdFdChidChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdFdChidGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdFdGradChidGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdChidChidChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdChidChidGradChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdChidGradChidGradChi_answer( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdGradChidGradChidGradChi_answer( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    // Tests of the gradients for the non-overlapped case. Everything should be zero!

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){

        floatVector delta( Xi_1.size( ), 0 );

        delta[ i ] += eps * std::fabs( Xi_1[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdXi_1_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dXi_1_answer[ j ][ Xi_1.size( ) * Xi_1.size( ) * k + Xi_1.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1dXi_1_answer ) );

    for ( unsigned int i = 0; i < dX.size( ); i++ ){

        floatVector delta( dX.size( ), 0 );

        delta[ i ] += eps * std::fabs( dX[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapddX_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1ddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1ddX_answer[ j ][ Xi_1.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1ddX, d2OverlapdXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX, d2OverlapddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXddX, d3OverlapddXddXddX_answer ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( R_nl ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdR_nl_answer[ j ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dR_nl_answer[ j ][ k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdR_nl_answer[ j ][ k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldR_nl_answer[ j + k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dR_nl_answer[ j ][ Xi_1.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldR_nl_answer[ j ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdR_nl, d3OverlapddXddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldR_nl_answer ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){

        floatVector delta( F.size( ), 0 );

        delta[ i ] += eps * std::fabs( F[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdF_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dF_answer[ j ][ Xi_1.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdF, d3OverlapddXddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldF, d3OverlapddXdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdF, d3OverlapddXdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdF, d3OverlapdR_nldFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdF, d3OverlapdFdFdF_answer ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi, d3OverlapddXdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi, d3OverlapddXdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi, d3OverlapdFdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi, d3OverlapdFdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi, d3OverlapdChidChidChi_answer ) );

    for ( unsigned int i = 0; i < gradChi.size( ); i++ ){

        floatVector delta( gradChi.size( ), 0 );

        delta[ i ] += eps * std::fabs( gradChi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi + delta, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi - delta, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdGradChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < gradChi.size( ); k++ ){

                d2OverlapdGradChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdGradChip[ j ][ k ] - dOverlapdGradChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi, dOverlapdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dGradChi, d2OverlapdXi_1dGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi, d2OverlapddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi, d2OverlapdR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi, d2OverlapdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi, d2OverlapdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi, d2OverlapdGradChidGradChi_answer ) );

    // Tests of the gradients of the overlapped case. Gradients will be non-zero

    R_nl = 2.5;

    overlap_answer = { -0.0149698, -0.00571608, 0.00344577 };

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap,
                                                              dOverlapdXi_1, dOverlapddX, dOverlapdR_nl, dOverlapdF, dOverlapdChi, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap_2,
                                                              dOverlapdXi_1_2, dOverlapddX_2, dOverlapdR_nl_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdGradChi_2,
                                                              d2OverlapdXi_1dXi_1,      d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF,       d2OverlapdXi_1dChi,   d2OverlapdXi_1dGradChi,
                                                              d2OverlapddXddX,          d2OverlapddXdR_nl, d2OverlapddXdF,      d2OverlapddXdChi,       d2OverlapddXdGradChi,
                                                              d2OverlapdR_nldR_nl,      d2OverlapdR_nldF,  d2OverlapdR_nldChi,  d2OverlapdR_nldGradChi,
                                                              d2OverlapdFdF,            d2OverlapdFdChi,   d2OverlapdFdGradChi,
                                                              d2OverlapdChidChi,        d2OverlapdChidGradChi,
                                                              d2OverlapdGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_2, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_2, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_2, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_2, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_2, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_2, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_2, dOverlapdGradChi ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi, overlap_3,
                                                              dOverlapdXi_1_3, dOverlapddX_3, dOverlapdR_nl_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdGradChi_3,
                                                              d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dGradChi_3,
                                                              d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdGradChi_3,
                                                              d2OverlapdR_nldR_nl_3, d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldGradChi_3,
                                                              d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdGradChi_3,
                                                              d2OverlapdChidChi_3, d2OverlapdChidGradChi_3,
                                                              d2OverlapdGradChidGradChi_3,
                                                              d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dGradChi,
                                                              d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdGradChi,
                                                              d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldGradChi,
                                                              d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdGradChi,
                                                              d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidGradChi,
                                                              d3OverlapdXi_1dGradChidGradChi,
                                                              d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdGradChi,
                                                              d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldGradChi,
                                                              d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdGradChi,
                                                              d3OverlapddXdChidChi, d3OverlapddXdChidGradChi,
                                                              d3OverlapddXdGradChidGradChi,
                                                              d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldGradChi,
                                                              d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdGradChi,
                                                              d3OverlapdR_nldChidChi, d3OverlapdR_nldChidGradChi,
                                                              d3OverlapdR_nldGradChidGradChi,
                                                              d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdGradChi,
                                                              d3OverlapdFdChidChi, d3OverlapdFdChidGradChi,
                                                              d3OverlapdFdGradChidGradChi,
                                                              d3OverlapdChidChidChi, d3OverlapdChidChidGradChi,
                                                              d3OverlapdChidGradChidGradChi,
                                                              d3OverlapdGradChidGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_3, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_3, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_3, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_3, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_3, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_3, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_3, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX_3, d2OverlapddXddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl_3, d2OverlapddXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF_3, d2OverlapddXdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_3, d2OverlapddXdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi_3, d2OverlapddXdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl_3, d2OverlapdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF_3, d2OverlapdR_nldF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_3, d2OverlapdR_nldChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi_3, d2OverlapdR_nldGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF_3, d2OverlapdFdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_3, d2OverlapdFdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi_3, d2OverlapdFdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_3, d2OverlapdChidChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi_3, d2OverlapdChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi_3, d2OverlapdGradChidGradChi ) );

    dOverlapdXi_1_answer    = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    dOverlapddX_answer      = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    dOverlapdR_nl_answer    = floatVector( overlap_answer.size( ), 0 );

    dOverlapdF_answer       = floatMatrix( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    dOverlapdChi_answer     = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    dOverlapdGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    d2OverlapdXi_1dXi_1_answer    = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d2OverlapdXi_1ddX_answer      = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    d2OverlapdXi_1dR_nl_answer    = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    d2OverlapdXi_1dF_answer       = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    d2OverlapdXi_1dChi_answer     = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    d2OverlapdXi_1dGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    d2OverlapddXddX_answer        = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    d2OverlapddXdR_nl_answer      = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    d2OverlapddXdF_answer         = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    d2OverlapddXdChi_answer       = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    d2OverlapddXdGradChi_answer   = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    d2OverlapdR_nldR_nl_answer    = floatVector( overlap_answer.size( ), 0 );

    d2OverlapdR_nldF_answer       = floatMatrix( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    d2OverlapdR_nldChi_answer     = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    d2OverlapdR_nldGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    d2OverlapdFdF_answer          = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    d2OverlapdFdChi_answer        = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    d2OverlapdFdGradChi_answer    = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    d2OverlapdChidChi_answer      = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    d2OverlapdChidGradChi_answer  = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    d2OverlapdGradChidGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dXi_1_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dXi_1ddX_answer            = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dR_nl_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dF_answer             = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dChi_answer           = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dGradChi_answer       = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1ddXddX_answer              = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1ddXdR_nl_answer            = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1ddXdF_answer               = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1ddXdChi_answer             = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1ddXdGradChi_answer         = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dR_nldR_nl_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dR_nldF_answer             = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dR_nldChi_answer           = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dR_nldGradChi_answer       = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dFdF_answer                = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dFdChi_answer              = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dFdGradChi_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dChidChi_answer            = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dChidGradChi_answer        = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dGradChidGradChi_answer    = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXddXddX_answer                = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * dX.size( ), 0 ) );

    d3OverlapddXddXdR_nl_answer              = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    d3OverlapddXddXdF_answer                 = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * F.size( ), 0 ) );

    d3OverlapddXddXdChi_answer               = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi.size( ), 0 ) );

    d3OverlapddXddXdGradChi_answer           = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdR_nldR_nl_answer            = floatMatrix( Xi_1.size( ), floatVector( dX.size( ), 0 ) );

    d3OverlapddXdR_nldF_answer               = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    d3OverlapddXdR_nldChi_answer             = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdR_nldGradChi_answer         = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdFdF_answer                  = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapddXdFdChi_answer                = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdFdGradChi_answer            = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdChidChi_answer              = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdChidGradChi_answer          = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdGradChidGradChi_answer      = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdR_nldR_nldR_nl_answer          = floatVector( Xi_1.size( ), 0 );

    d3OverlapdR_nldR_nldF_answer             = floatMatrix( Xi_1.size( ), floatVector( F.size( ), 0 ) );

    d3OverlapdR_nldR_nldChi_answer           = floatMatrix( Xi_1.size( ), floatVector( chi.size( ), 0 ) );

    d3OverlapdR_nldR_nldGradChi_answer       = floatMatrix( Xi_1.size( ), floatVector( gradChi.size( ), 0 ) );

    d3OverlapdR_nldFdF_answer                = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    d3OverlapdR_nldFdChi_answer              = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    d3OverlapdR_nldFdGradChi_answer          = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdR_nldChidChi_answer            = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdR_nldChidGradChi_answer        = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdR_nldGradChidGradChi_answer    = floatMatrix( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdFdFdF_answer                   = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapdFdFdChi_answer                 = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapdFdFdGradChi_answer             = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdFdChidChi_answer               = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdFdChidGradChi_answer           = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdFdGradChidGradChi_answer       = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdChidChidChi_answer             = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdChidChidGradChi_answer         = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdChidGradChidGradChi_answer     = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdGradChidGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){

        floatVector delta( Xi_1.size( ), 0 );

        delta[ i ] += eps * std::fabs( Xi_1[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdXi_1_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dXi_1_answer[ j ][ Xi_1.size( ) * Xi_1.size( ) * k + Xi_1.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1dXi_1_answer ) );

    for ( unsigned int i = 0; i < dX.size( ); i++ ){

        floatVector delta( dX.size( ), 0 );

        delta[ i ] += eps * std::fabs( dX[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapddX_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1ddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1ddX_answer[ j ][ Xi_1.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1ddX, d2OverlapdXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX, d2OverlapddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXddX, d3OverlapddXddXddX_answer ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( R_nl ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdR_nl_answer[ j ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dR_nl_answer[ j ][ k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdR_nl_answer[ j ][ k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldR_nl_answer[ j + k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dR_nl_answer[ j ][ Xi_1.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldR_nl_answer[ j ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdR_nl, d3OverlapddXddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldR_nl_answer ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){

        floatVector delta( F.size( ), 0 );

        delta[ i ] += eps * std::fabs( F[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdF_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dF_answer[ j ][ Xi_1.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdF, d3OverlapddXddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldF, d3OverlapddXdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdF, d3OverlapddXdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdF, d3OverlapdR_nldFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdF, d3OverlapdFdFdF_answer ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi, d3OverlapddXdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi, d3OverlapddXdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi, d3OverlapdFdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi, d3OverlapdFdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi, d3OverlapdChidChidChi_answer ) );

    for ( unsigned int i = 0; i < gradChi.size( ); i++ ){

        floatVector delta( gradChi.size( ), 0 );

        delta[ i ] += eps * std::fabs( gradChi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi + delta, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi - delta, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdGradChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < gradChi.size( ); k++ ){

                d2OverlapdGradChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdGradChip[ j ][ k ] - dOverlapdGradChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,        d2OverlapdChidGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, gradChi - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,        d2OverlapdChidGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dGradChi_answer[ j ][ Xi_1.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdGradChi_answer[ j ][ dX.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdXi_1dGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdXi_1dGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdGradChi_answer[ j ][ dX.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapddXdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapddXdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ chi.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdR_nldGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdR_nldGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdFdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdFdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdChidGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdChidGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < gradChi.size( ); k++ ){

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdGradChidGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdGradChidGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdGradChidGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi, dOverlapdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dGradChi, d2OverlapdXi_1dGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi, d2OverlapddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi, d2OverlapdR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi, d2OverlapdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi, d2OverlapdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi, d2OverlapdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dGradChi, d3OverlapdXi_1dXi_1dGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdGradChi, d3OverlapdXi_1ddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldGradChi, d3OverlapdXi_1dR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdGradChi, d3OverlapdXi_1dFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidGradChi, d3OverlapdXi_1dChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dGradChidGradChi, d3OverlapdXi_1dGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdGradChi, d3OverlapddXddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldGradChi, d3OverlapddXdR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdGradChi, d3OverlapddXdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidGradChi, d3OverlapddXdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdGradChidGradChi, d3OverlapddXdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldGradChi, d3OverlapdR_nldR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdGradChi, d3OverlapdR_nldFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidGradChi, d3OverlapdR_nldChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldGradChidGradChi, d3OverlapdR_nldGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdGradChi, d3OverlapdFdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidGradChi, d3OverlapdFdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdGradChidGradChi, d3OverlapdFdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidGradChi, d3OverlapdChidChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidGradChidGradChi, d3OverlapdChidGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdGradChidGradChidGradChi, d3OverlapdGradChidGradChidGradChi_answer ) );

}

BOOST_AUTO_TEST_CASE( test_computeParticleOverlapChi_nl ){

    floatVector Xi_1 = { 1, 0, 0 };

    floatVector dX = { 2, 0, 0 };

    floatType R_nl = 1;

    floatVector F = { 0.75, 0.0, 0.0,
                      0.00, 1.0, 0.0,
                      0.00, 0.0, 1.0 };

    floatVector chi = { 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0 };

    floatVector chi_nl = chi;

    floatVector overlap_answer = { -0.5, 0.0, 0.0 };

    floatVector overlap;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    Xi_1 = { 0.39293837, -0.42772133, -0.54629709 };

    dX = { 0.10262954,  0.43893794, -0.15378708 };

    R_nl = 1.961528396769231;

    F = { 1.36965948, -0.0381362 , -0.21576496,
         -0.31364397,  1.45809941, -0.12285551,
         -0.88064421, -0.20391149,  1.47599081 };

    chi = { 0.36498346, -0.64909649,  0.06310275,
            0.06365517,  1.26880192,  0.69886359,
            0.44891065,  0.22204702,  1.44488677 };

    chi_nl = {  0.29089097, -0.45099979, -0.0092573 ,
               -0.09372774,  1.73890107,  0.44306498,
                0.82975541,  0.63883275,  1.03464445 };

    overlap_answer = floatVector( Xi_1.size( ), 0 );

    overlap.clear( );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    floatMatrix dOverlapdXi_1, dOverlapddX, dOverlapdF, dOverlapdChi, dOverlapdChi_nl;

    floatVector dOverlapdR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap,
                                                                    dOverlapdXi_1, dOverlapddX, dOverlapdR_nl, dOverlapdF, dOverlapdChi, dOverlapdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    floatVector overlap_2;

    floatMatrix dOverlapdXi_1_2, dOverlapddX_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdChi_nl_2;

    floatVector dOverlapdR_nl_2;

    floatMatrix d2OverlapdXi_1dXi_1, d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF, d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_nl,
                d2OverlapddXddX, d2OverlapddXdR_nl, d2OverlapddXdF, d2OverlapddXdChi, d2OverlapddXdChi_nl,
                d2OverlapdR_nldF, d2OverlapdR_nldChi, d2OverlapdR_nldChi_nl,
                d2OverlapdFdF, d2OverlapdFdChi, d2OverlapdFdChi_nl,
                d2OverlapdChidChi, d2OverlapdChidChi_nl,
                d2OverlapdChi_nldChi_nl;

    floatVector d2OverlapdR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap_2,
                                                                    dOverlapdXi_1_2, dOverlapddX_2, dOverlapdR_nl_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdChi_nl_2,
                                                                    d2OverlapdXi_1dXi_1,      d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF,       d2OverlapdXi_1dChi,   d2OverlapdXi_1dChi_nl,
                                                                    d2OverlapddXddX,          d2OverlapddXdR_nl, d2OverlapddXdF,      d2OverlapddXdChi,       d2OverlapddXdChi_nl,
                                                                    d2OverlapdR_nldR_nl,      d2OverlapdR_nldF,  d2OverlapdR_nldChi,  d2OverlapdR_nldChi_nl,
                                                                    d2OverlapdFdF,            d2OverlapdFdChi,   d2OverlapdFdChi_nl,
                                                                    d2OverlapdChidChi,        d2OverlapdChidChi_nl,
                                                                    d2OverlapdChi_nldChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_2, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_2, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_2, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_2, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_2, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_2, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_nl_2, dOverlapdChi_nl ) );

    floatMatrix dOverlapdXi_1_answer( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix dOverlapddX_answer( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    floatVector dOverlapdR_nl_answer( overlap_answer.size( ), 0 );

    floatMatrix dOverlapdF_answer( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix dOverlapdChi_answer( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix dOverlapdChi_nl_answer( overlap_answer.size( ), floatVector( chi_nl.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dXi_1_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1ddX_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dR_nl_answer( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dF_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dChi_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dChi_nl_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d2OverlapddXddX_answer( overlap_answer.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d2OverlapddXdR_nl_answer( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    floatMatrix d2OverlapddXdF_answer( overlap_answer.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapddXdChi_answer( overlap_answer.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapddXdChi_nl_answer( overlap_answer.size( ), floatVector( dX.size( ) * chi_nl.size( ), 0 ) );

    floatVector d2OverlapdR_nldR_nl_answer( overlap_answer.size( ), 0 );

    floatMatrix d2OverlapdR_nldF_answer( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix d2OverlapdR_nldChi_answer( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix d2OverlapdR_nldChi_nl_answer( overlap_answer.size( ), floatVector( chi_nl.size( ), 0 ) );

    floatMatrix d2OverlapdFdF_answer( overlap_answer.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapdFdChi_answer( overlap_answer.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdFdChi_nl_answer( overlap_answer.size( ), floatVector( F.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d2OverlapdChidChi_answer( overlap_answer.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdChidChi_nl_answer( overlap_answer.size( ), floatVector( chi.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d2OverlapdChi_nldChi_nl_answer( overlap_answer.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatVector overlap_3;

    floatMatrix dOverlapdXi_1_3, dOverlapddX_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdChi_nl_3;

    floatVector dOverlapdR_nl_3;

    floatMatrix d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dChi_nl_3,
                d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdChi_nl_3,
                d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldChi_nl_3,
                d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdChi_nl_3,
                d2OverlapdChidChi_3, d2OverlapdChidChi_nl_3,
                d2OverlapdChi_nldChi_nl_3;

    floatVector d2OverlapdR_nldR_nl_3;

    floatMatrix d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_nl,
                d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_nl,
                d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_nl,
                d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_nl,
                d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_nl,
                d3OverlapdXi_1dChi_nldChi_nl,
                d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdChi_nl,
                d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_nl,
                d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdChi_nl,
                d3OverlapddXdChidChi, d3OverlapddXdChidChi_nl,
                d3OverlapddXdChi_nldChi_nl,
                d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_nl,
                d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_nl,
                d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_nl,
                d3OverlapdR_nldChi_nldChi_nl,
                d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdChi_nl,
                d3OverlapdFdChidChi, d3OverlapdFdChidChi_nl,
                d3OverlapdFdChi_nldChi_nl,
                d3OverlapdChidChidChi, d3OverlapdChidChidChi_nl,
                d3OverlapdChidChi_nldChi_nl,
                d3OverlapdChi_nldChi_nldChi_nl;

    floatVector d3OverlapdR_nldR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap_3,
                                                                    dOverlapdXi_1_3, dOverlapddX_3, dOverlapdR_nl_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdChi_nl_3,
                                                                    d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dChi_nl_3,
                                                                    d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdChi_nl_3,
                                                                    d2OverlapdR_nldR_nl_3, d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldChi_nl_3,
                                                                    d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdChi_nl_3,
                                                                    d2OverlapdChidChi_3, d2OverlapdChidChi_nl_3,
                                                                    d2OverlapdChi_nldChi_nl_3,
                                                                    d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_nl,
                                                                    d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_nl,
                                                                    d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_nl,
                                                                    d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_nl,
                                                                    d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_nl,
                                                                    d3OverlapdXi_1dChi_nldChi_nl,
                                                                    d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdChi_nl,
                                                                    d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_nl,
                                                                    d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdChi_nl,
                                                                    d3OverlapddXdChidChi, d3OverlapddXdChidChi_nl,
                                                                    d3OverlapddXdChi_nldChi_nl,
                                                                    d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_nl,
                                                                    d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_nl,
                                                                    d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_nl,
                                                                    d3OverlapdR_nldChi_nldChi_nl,
                                                                    d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdChi_nl,
                                                                    d3OverlapdFdChidChi, d3OverlapdFdChidChi_nl,
                                                                    d3OverlapdFdChi_nldChi_nl,
                                                                    d3OverlapdChidChidChi, d3OverlapdChidChidChi_nl,
                                                                    d3OverlapdChidChi_nldChi_nl,
                                                                    d3OverlapdChi_nldChi_nldChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_3, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_3, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_3, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_3, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_3, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_3, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_nl_3, dOverlapdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX_3, d2OverlapddXddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl_3, d2OverlapddXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF_3, d2OverlapddXdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_3, d2OverlapddXdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_nl_3, d2OverlapddXdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl_3, d2OverlapdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF_3, d2OverlapdR_nldF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_3, d2OverlapdR_nldChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_nl_3, d2OverlapdR_nldChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF_3, d2OverlapdFdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_3, d2OverlapdFdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_nl_3, d2OverlapdFdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_3, d2OverlapdChidChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_nl_3, d2OverlapdChidChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_nldChi_nl_3, d2OverlapdChi_nldChi_nl ) );

    floatMatrix d3OverlapdXi_1dXi_1dXi_1_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1ddX_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dChi_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXddX_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdChi_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldChi_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdChi_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChidChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChidChi_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChi_nldChi_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapddXddXddX_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdR_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdF_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdChi_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdChi_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldR_nl_answer( Xi_1.size( ), floatVector( dX.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldF_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldChi_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdF_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdChi_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdChi_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapddXdChidChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdChidChi_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapddXdChi_nldChi_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatVector d3OverlapdR_nldR_nldR_nl_answer( Xi_1.size( ), 0 );

    floatMatrix d3OverlapdR_nldR_nldF_answer( Xi_1.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldR_nldChi_answer( Xi_1.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldR_nldChi_nl_answer( Xi_1.size( ), floatVector( chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdF_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdChi_nl_answer( Xi_1.size( ), floatVector( F.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChidChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChidChi_nl_answer( Xi_1.size( ), floatVector( chi.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChi_nldChi_nl_answer( Xi_1.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdF_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdChi_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdChi_nl_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdFdChidChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdFdChidChi_nl_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdFdChi_nldChi_nl_answer( Xi_1.size( ), floatVector( F.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdChidChidChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdChidChidChi_nl_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdChidChi_nldChi_nl_answer( Xi_1.size( ), floatVector( chi.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    floatMatrix d3OverlapdChi_nldChi_nldChi_nl_answer( Xi_1.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    // Tests of the gradients for the non-overlapped case. Everything should be zero!

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){

        floatVector delta( Xi_1.size( ), 0 );

        delta[ i ] += eps * std::fabs( Xi_1[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 + delta, dX, R_nl, F, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 - delta, dX, R_nl, F, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdXi_1_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 + delta, dX, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 - delta, dX, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 + delta, dX, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 - delta, dX, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dXi_1_answer[ j ][ Xi_1.size( ) * Xi_1.size( ) * k + Xi_1.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1dXi_1_answer ) );

    for ( unsigned int i = 0; i < dX.size( ); i++ ){

        floatVector delta( dX.size( ), 0 );

        delta[ i ] += eps * std::fabs( dX[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX + delta, R_nl, F, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX - delta, R_nl, F, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapddX_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX + delta, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX - delta, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1ddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX + delta, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX - delta, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1ddX_answer[ j ][ Xi_1.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1ddX, d2OverlapdXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX, d2OverlapddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXddX, d3OverlapddXddXddX_answer ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( R_nl ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl + delta, F, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl - delta, F, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdR_nl_answer[ j ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl + delta, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl - delta, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dR_nl_answer[ j ][ k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdR_nl_answer[ j ][ k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldR_nl_answer[ j + k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl + delta, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl - delta, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dR_nl_answer[ j ][ Xi_1.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldR_nl_answer[ j ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdR_nl, d3OverlapddXddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldR_nl_answer ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){

        floatVector delta( F.size( ), 0 );

        delta[ i ] += eps * std::fabs( F[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F + delta, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F - delta, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdF_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F + delta, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F - delta, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F + delta, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F - delta, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dF_answer[ j ][ Xi_1.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdF, d3OverlapddXddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldF, d3OverlapddXdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdF, d3OverlapddXdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdF, d3OverlapdR_nldFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdF, d3OverlapdFdFdF_answer ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi + delta, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi - delta, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi + delta, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi - delta, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi + delta, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi - delta, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi, d3OverlapddXdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi, d3OverlapddXdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi, d3OverlapdFdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi, d3OverlapdFdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi, d3OverlapdChidChidChi_answer ) );

    for ( unsigned int i = 0; i < chi_nl.size( ); i++ ){

        floatVector delta( chi_nl.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi_nl[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl + delta, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl - delta, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_nl_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl + delta, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl - delta, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d2OverlapdChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdChi_nlp[ j ][ k ] - dOverlapdChi_nlm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_nl, dOverlapdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi_nl, d2OverlapdXi_1dChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_nl, d2OverlapddXdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_nl, d2OverlapdR_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_nl, d2OverlapdFdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_nl, d2OverlapdChidChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_nldChi_nl, d2OverlapdChi_nldChi_nl_answer ) );

    // Tests of the gradients of the overlapped case. Gradients will be non-zero

    R_nl = 2.5;

    overlap_answer = { -0.0149698, -0.00571608, 0.00344577 };

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap,
                                                                    dOverlapdXi_1, dOverlapddX, dOverlapdR_nl, dOverlapdF, dOverlapdChi, dOverlapdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap_2,
                                                                    dOverlapdXi_1_2, dOverlapddX_2, dOverlapdR_nl_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdChi_nl_2,
                                                                    d2OverlapdXi_1dXi_1,      d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF,       d2OverlapdXi_1dChi,   d2OverlapdXi_1dChi_nl,
                                                                    d2OverlapddXddX,          d2OverlapddXdR_nl, d2OverlapddXdF,      d2OverlapddXdChi,       d2OverlapddXdChi_nl,
                                                                    d2OverlapdR_nldR_nl,      d2OverlapdR_nldF,  d2OverlapdR_nldChi,  d2OverlapdR_nldChi_nl,
                                                                    d2OverlapdFdF,            d2OverlapdFdChi,   d2OverlapdFdChi_nl,
                                                                    d2OverlapdChidChi,        d2OverlapdChidChi_nl,
                                                                    d2OverlapdChi_nldChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_2, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_2, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_2, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_2, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_2, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_2, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_nl_2, dOverlapdChi_nl ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl, overlap_3,
                                                                    dOverlapdXi_1_3, dOverlapddX_3, dOverlapdR_nl_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdChi_nl_3,
                                                                    d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dChi_nl_3,
                                                                    d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdChi_nl_3,
                                                                    d2OverlapdR_nldR_nl_3, d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldChi_nl_3,
                                                                    d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdChi_nl_3,
                                                                    d2OverlapdChidChi_3, d2OverlapdChidChi_nl_3,
                                                                    d2OverlapdChi_nldChi_nl_3,
                                                                    d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_nl,
                                                                    d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_nl,
                                                                    d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_nl,
                                                                    d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_nl,
                                                                    d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_nl,
                                                                    d3OverlapdXi_1dChi_nldChi_nl,
                                                                    d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdChi_nl,
                                                                    d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_nl,
                                                                    d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdChi_nl,
                                                                    d3OverlapddXdChidChi, d3OverlapddXdChidChi_nl,
                                                                    d3OverlapddXdChi_nldChi_nl,
                                                                    d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_nl,
                                                                    d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_nl,
                                                                    d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_nl,
                                                                    d3OverlapdR_nldChi_nldChi_nl,
                                                                    d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdChi_nl,
                                                                    d3OverlapdFdChidChi, d3OverlapdFdChidChi_nl,
                                                                    d3OverlapdFdChi_nldChi_nl,
                                                                    d3OverlapdChidChidChi, d3OverlapdChidChidChi_nl,
                                                                    d3OverlapdChidChi_nldChi_nl,
                                                                    d3OverlapdChi_nldChi_nldChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_3, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_3, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_3, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_3, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_3, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_3, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_nl_3, dOverlapdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX_3, d2OverlapddXddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl_3, d2OverlapddXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF_3, d2OverlapddXdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_3, d2OverlapddXdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_nl_3, d2OverlapddXdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl_3, d2OverlapdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF_3, d2OverlapdR_nldF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_3, d2OverlapdR_nldChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_nl_3, d2OverlapdR_nldChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF_3, d2OverlapdFdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_3, d2OverlapdFdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_nl_3, d2OverlapdFdChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_3, d2OverlapdChidChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_nl_3, d2OverlapdChidChi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_nldChi_nl_3, d2OverlapdChi_nldChi_nl ) );

    dOverlapdXi_1_answer    = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    dOverlapddX_answer      = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    dOverlapdR_nl_answer    = floatVector( overlap_answer.size( ), 0 );

    dOverlapdF_answer       = floatMatrix( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    dOverlapdChi_answer     = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    dOverlapdChi_nl_answer = floatMatrix( overlap_answer.size( ), floatVector( chi_nl.size( ), 0 ) );

    d2OverlapdXi_1dXi_1_answer    = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d2OverlapdXi_1ddX_answer      = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    d2OverlapdXi_1dR_nl_answer    = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    d2OverlapdXi_1dF_answer       = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    d2OverlapdXi_1dChi_answer     = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    d2OverlapdXi_1dChi_nl_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi_nl.size( ), 0 ) );

    d2OverlapddXddX_answer        = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    d2OverlapddXdR_nl_answer      = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    d2OverlapddXdF_answer         = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    d2OverlapddXdChi_answer       = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    d2OverlapddXdChi_nl_answer   = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * chi_nl.size( ), 0 ) );

    d2OverlapdR_nldR_nl_answer    = floatVector( overlap_answer.size( ), 0 );

    d2OverlapdR_nldF_answer       = floatMatrix( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    d2OverlapdR_nldChi_answer     = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    d2OverlapdR_nldChi_nl_answer = floatMatrix( overlap_answer.size( ), floatVector( chi_nl.size( ), 0 ) );

    d2OverlapdFdF_answer          = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    d2OverlapdFdChi_answer        = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    d2OverlapdFdChi_nl_answer    = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * chi_nl.size( ), 0 ) );

    d2OverlapdChidChi_answer      = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    d2OverlapdChidChi_nl_answer  = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ) * chi_nl.size( ), 0 ) );

    d2OverlapdChi_nldChi_nl_answer = floatMatrix( overlap_answer.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dXi_1_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dXi_1ddX_answer            = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dR_nl_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dF_answer             = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dChi_answer           = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dChi_nl_answer       = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdXi_1ddXddX_answer              = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1ddXdR_nl_answer            = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1ddXdF_answer               = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1ddXdChi_answer             = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1ddXdChi_nl_answer         = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdXi_1dR_nldR_nl_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dR_nldF_answer             = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dR_nldChi_answer           = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dR_nldChi_nl_answer       = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdXi_1dFdF_answer                = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dFdChi_answer              = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dFdChi_nl_answer          = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdXi_1dChidChi_answer            = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dChidChi_nl_answer        = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdXi_1dChi_nldChi_nl_answer    = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapddXddXddX_answer                = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * dX.size( ), 0 ) );

    d3OverlapddXddXdR_nl_answer              = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    d3OverlapddXddXdF_answer                 = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * F.size( ), 0 ) );

    d3OverlapddXddXdChi_answer               = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi.size( ), 0 ) );

    d3OverlapddXddXdChi_nl_answer           = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapddXdR_nldR_nl_answer            = floatMatrix( Xi_1.size( ), floatVector( dX.size( ), 0 ) );

    d3OverlapddXdR_nldF_answer               = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    d3OverlapddXdR_nldChi_answer             = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdR_nldChi_nl_answer         = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapddXdFdF_answer                  = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapddXdFdChi_answer                = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdFdChi_nl_answer            = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapddXdChidChi_answer              = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdChidChi_nl_answer          = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapddXdChi_nldChi_nl_answer      = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdR_nldR_nldR_nl_answer          = floatVector( Xi_1.size( ), 0 );

    d3OverlapdR_nldR_nldF_answer             = floatMatrix( Xi_1.size( ), floatVector( F.size( ), 0 ) );

    d3OverlapdR_nldR_nldChi_answer           = floatMatrix( Xi_1.size( ), floatVector( chi.size( ), 0 ) );

    d3OverlapdR_nldR_nldChi_nl_answer       = floatMatrix( Xi_1.size( ), floatVector( chi_nl.size( ), 0 ) );

    d3OverlapdR_nldFdF_answer                = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    d3OverlapdR_nldFdChi_answer              = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    d3OverlapdR_nldFdChi_nl_answer          = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdR_nldChidChi_answer            = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdR_nldChidChi_nl_answer        = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdR_nldChi_nldChi_nl_answer    = floatMatrix( Xi_1.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdFdFdF_answer                   = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapdFdFdChi_answer                 = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapdFdFdChi_nl_answer             = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdFdChidChi_answer               = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdFdChidChi_nl_answer           = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdFdChi_nldChi_nl_answer       = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdChidChidChi_answer             = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdChidChidChi_nl_answer         = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdChidChi_nldChi_nl_answer     = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    d3OverlapdChi_nldChi_nldChi_nl_answer = floatMatrix( Xi_1.size( ), floatVector( chi_nl.size( ) * chi_nl.size( ) * chi_nl.size( ), 0 ) );

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){

        floatVector delta( Xi_1.size( ), 0 );

        delta[ i ] += eps * std::fabs( Xi_1[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 + delta, dX, R_nl, F, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 - delta, dX, R_nl, F, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdXi_1_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 + delta, dX, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 - delta, dX, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 + delta, dX, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1 - delta, dX, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dXi_1_answer[ j ][ Xi_1.size( ) * Xi_1.size( ) * k + Xi_1.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1dXi_1_answer ) );

    for ( unsigned int i = 0; i < dX.size( ); i++ ){

        floatVector delta( dX.size( ), 0 );

        delta[ i ] += eps * std::fabs( dX[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX + delta, R_nl, F, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX - delta, R_nl, F, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapddX_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX + delta, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX - delta, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1ddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX + delta, R_nl, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX - delta, R_nl, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1ddX_answer[ j ][ Xi_1.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1ddX, d2OverlapdXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX, d2OverlapddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXddX, d3OverlapddXddXddX_answer ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( R_nl ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl + delta, F, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl - delta, F, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdR_nl_answer[ j ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl + delta, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl - delta, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dR_nl_answer[ j ][ k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdR_nl_answer[ j ][ k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldR_nl_answer[ j + k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl + delta, F, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl - delta, F, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dR_nl_answer[ j ][ Xi_1.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldR_nl_answer[ j ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdR_nl, d3OverlapddXddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldR_nl_answer ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){

        floatVector delta( F.size( ), 0 );

        delta[ i ] += eps * std::fabs( F[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F + delta, chi, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F - delta, chi, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdF_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F + delta, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F - delta, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F + delta, chi, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                        d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                        d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                        d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                        d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                        d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                        d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F - delta, chi, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                        d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                        d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                        d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                        d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                        d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                        d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dF_answer[ j ][ Xi_1.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdF, d3OverlapddXddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldF, d3OverlapddXdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdF, d3OverlapddXdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdF, d3OverlapdR_nldFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdF, d3OverlapdFdFdF_answer ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi + delta, chi_nl, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi - delta, chi_nl, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi + delta, chi_nl, overlapp,
                                                                        dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi - delta, chi_nl, overlapm,
                                                                        dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi + delta, chi_nl, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                  d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                  d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi - delta, chi_nl, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                  d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                  d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi, d3OverlapddXdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi, d3OverlapddXdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi, d3OverlapdFdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi, d3OverlapdFdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi, d3OverlapdChidChidChi_answer ) );

    for ( unsigned int i = 0; i < chi_nl.size( ); i++ ){

        floatVector delta( chi_nl.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi_nl[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl + delta, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl - delta, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_nl_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_nlp, dOverlapdChi_nlm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                d2OverlapdChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * k + i ] = ( dOverlapdChi_nlp[ j ][ k ] - dOverlapdChi_nlm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_nlp,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_nlp,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_nlp,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_nlp,
                    d2OverlapdChidChip, d2OverlapdChidChi_nlp,
                    d2OverlapdChi_nldChi_nlp;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_nlm,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_nlm,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_nlm,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_nlm,
                    d2OverlapdChidChim, d2OverlapdChidChi_nlm,
                    d2OverlapdChi_nldChi_nlm;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_nlp,
                                                                  d2OverlapdXi_1dXi_1p,      d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp,       d2OverlapdXi_1dChip,   d2OverlapdXi_1dChi_nlp,
                                                                  d2OverlapddXddXp,          d2OverlapddXdR_nlp, d2OverlapddXdFp,      d2OverlapddXdChip,       d2OverlapddXdChi_nlp,
                                                                  d2OverlapdR_nldR_nlp,      d2OverlapdR_nldFp,  d2OverlapdR_nldChip,  d2OverlapdR_nldChi_nlp,
                                                                  d2OverlapdFdFp,            d2OverlapdFdChip,   d2OverlapdFdChi_nlp,
                                                                  d2OverlapdChidChip,        d2OverlapdChidChi_nlp,
                                                                  d2OverlapdChi_nldChi_nlp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlapChi_nl( Xi_1, dX, R_nl, F, chi, chi_nl - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_nlm,
                                                                  d2OverlapdXi_1dXi_1m,      d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm,       d2OverlapdXi_1dChim,   d2OverlapdXi_1dChi_nlm,
                                                                  d2OverlapddXddXm,          d2OverlapddXdR_nlm, d2OverlapddXdFm,      d2OverlapddXdChim,       d2OverlapddXdChi_nlm,
                                                                  d2OverlapdR_nldR_nlm,      d2OverlapdR_nldFm,  d2OverlapdR_nldChim,  d2OverlapdR_nldChi_nlm,
                                                                  d2OverlapdFdFm,            d2OverlapdFdChim,   d2OverlapdFdChi_nlm,
                                                                  d2OverlapdChidChim,        d2OverlapdChidChi_nlm,
                                                                  d2OverlapdChi_nldChi_nlm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_nl_answer[ j ][ Xi_1.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_nl_answer[ j ][ dX.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_nl_answer[ j ][ chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_nl_answer[ j ][ F.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_nl_answer[ j ][ chi.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3OverlapdXi_1dChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2OverlapdXi_1dChi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_nl_answer[ j ][ dX.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_nl_answer[ j ][ chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_nl_answer[ j ][ F.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_nl_answer[ j ][ chi.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3OverlapddXdChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapddXdChi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2OverlapddXdChi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_nl_answer[ j ][ chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_nl_answer[ j ][ F.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_nl_answer[ j ][ chi.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ chi.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3OverlapdR_nldChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdR_nldChi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2OverlapdR_nldChi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_nl_answer[ j ][ F.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_nl_answer[ j ][ chi.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3OverlapdFdChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdFdChi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2OverlapdFdChi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_nl_answer[ j ][ chi.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3OverlapdChidChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdChidChi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2OverlapdChidChi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    d3OverlapdChi_nldChi_nldChi_nl_answer[ j ][ chi_nl.size( ) * chi_nl.size( ) * k + chi_nl.size( ) * l + i ]
                        = ( d2OverlapdChi_nldChi_nlp[ j ][ chi_nl.size( ) * k + l ] - d2OverlapdChi_nldChi_nlm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_nl, dOverlapdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi_nl, d2OverlapdXi_1dChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_nl, d2OverlapddXdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_nl, d2OverlapdR_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_nl, d2OverlapdFdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_nl, d2OverlapdChidChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_nldChi_nl, d2OverlapdChi_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi_nl, d3OverlapdXi_1dXi_1dChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi_nl, d3OverlapdXi_1ddXdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi_nl, d3OverlapdXi_1dR_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi_nl, d3OverlapdXi_1dFdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi_nl, d3OverlapdXi_1dChidChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChi_nldChi_nl, d3OverlapdXi_1dChi_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi_nl, d3OverlapddXddXdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi_nl, d3OverlapddXdR_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi_nl, d3OverlapddXdFdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi_nl, d3OverlapddXdChidChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChi_nldChi_nl, d3OverlapddXdChi_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi_nl, d3OverlapdR_nldR_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi_nl, d3OverlapdR_nldFdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi_nl, d3OverlapdR_nldChidChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChi_nldChi_nl, d3OverlapdR_nldChi_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi_nl, d3OverlapdFdFdChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi_nl, d3OverlapdFdChidChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChi_nldChi_nl, d3OverlapdFdChi_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi_nl, d3OverlapdChidChidChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChi_nldChi_nl, d3OverlapdChidChi_nldChi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChi_nldChi_nldChi_nl, d3OverlapdChi_nldChi_nldChi_nl_answer ) );

}

BOOST_AUTO_TEST_CASE( test_computeParticleOverlap_2 ){

    floatVector Xi_1 = { 1, 0, 0 };

    floatVector dX = { 2, 0, 0 };

    floatType R_nl = 1;

    floatVector F = { 0.75, 0.0, 0.0,
                      0.00, 1.0, 0.0,
                      0.00, 0.0, 1.0 };

    floatVector chi_nl_basis = { 1.50, 0.0, 0.0,
                                 0.00, 1.0, 0.0,
                                 0.00, 0.0, 1.0 };

    floatVector chi = { 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0 };

    floatVector gradChi( chi.size( ) * dX.size( ) );

    floatVector overlap_answer = { -1.0, 0.0, 0.0 };

    floatVector overlap;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    Xi_1 = { 0.39293837, -0.42772133, -0.54629709 };

    dX = { 0.10262954,  0.43893794, -0.15378708 };

    R_nl = 1.961528396769231;

    F = { 1.36965948, -0.0381362 , -0.21576496,
         -0.31364397,  1.45809941, -0.12285551,
         -0.88064421, -0.20391149,  1.47599081 };

    chi = { 0.36498346, -0.64909649,  0.06310275,
            0.06365517,  1.26880192,  0.69886359,
            0.44891065,  0.22204702,  1.44488677 };

    chi_nl_basis = chi;

    gradChi = { -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225,
                -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942,
                -0.37547755, -0.14729739,  0.78677833,  0.88832004,  0.00367335,
                 0.2479059 , -0.76876321, -0.36542904, -0.17034758,  0.73261832,
                -0.49908927, -0.03393147,  0.97111957,  0.03897024,  0.22578905,
                -0.75874267,  0.6526816 };

    overlap_answer = floatVector( Xi_1.size( ), 0 );

    overlap.clear( );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    floatMatrix dOverlapdXi_1, dOverlapddX, dOverlapdF, dOverlapdChi, dOverlapdChi_NL_B, dOverlapdGradChi;

    floatVector dOverlapdR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap,
                                                              dOverlapdXi_1, dOverlapddX, dOverlapdR_nl, dOverlapdF, dOverlapdChi, dOverlapdChi_NL_B, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    floatVector overlap_2;

    floatMatrix dOverlapdXi_1_2, dOverlapddX_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdChi_NL_B_2, dOverlapdGradChi_2;

    floatVector dOverlapdR_nl_2;

    floatMatrix d2OverlapdXi_1dXi_1, d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF, d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_NL_B, d2OverlapdXi_1dGradChi,
                d2OverlapddXddX, d2OverlapddXdR_nl, d2OverlapddXdF, d2OverlapddXdChi, d2OverlapddXdChi_NL_B, d2OverlapddXdGradChi,
                d2OverlapdR_nldF, d2OverlapdR_nldChi, d2OverlapdR_nldChi_NL_B, d2OverlapdR_nldGradChi,
                d2OverlapdFdF, d2OverlapdFdChi, d2OverlapdFdChi_NL_B, d2OverlapdFdGradChi,
                d2OverlapdChidChi, d2OverlapdChidChi_NL_B, d2OverlapdChidGradChi,
                d2OverlapdChi_NL_BdChi_NL_B, d2OverlapdChi_NL_BdGradChi,
                d2OverlapdGradChidGradChi;

    floatVector d2OverlapdR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap_2,
                                                              dOverlapdXi_1_2, dOverlapddX_2, dOverlapdR_nl_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdChi_NL_B_2, dOverlapdGradChi_2,
                                                              d2OverlapdXi_1dXi_1, d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF, d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_NL_B, d2OverlapdXi_1dGradChi,
                                                              d2OverlapddXddX, d2OverlapddXdR_nl, d2OverlapddXdF, d2OverlapddXdChi, d2OverlapddXdChi_NL_B, d2OverlapddXdGradChi,
                                                              d2OverlapdR_nldR_nl, d2OverlapdR_nldF, d2OverlapdR_nldChi, d2OverlapdR_nldChi_NL_B, d2OverlapdR_nldGradChi,
                                                              d2OverlapdFdF, d2OverlapdFdChi, d2OverlapdFdChi_NL_B, d2OverlapdFdGradChi,
                                                              d2OverlapdChidChi, d2OverlapdChidChi_NL_B, d2OverlapdChidGradChi,
                                                              d2OverlapdChi_NL_BdChi_NL_B, d2OverlapdChi_NL_BdGradChi,
                                                              d2OverlapdGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_2, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_2, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_2, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_2, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_2, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_2, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_NL_B_2, dOverlapdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_2, dOverlapdGradChi ) );

    floatMatrix dOverlapdXi_1_answer( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix dOverlapddX_answer( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    floatVector dOverlapdR_nl_answer( overlap_answer.size( ), 0 );

    floatMatrix dOverlapdF_answer( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix dOverlapdChi_answer( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix dOverlapdChi_NL_B_answer( overlap_answer.size( ), floatVector( chi_nl_basis.size( ), 0 ) );

    floatMatrix dOverlapdGradChi_answer( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dXi_1_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1ddX_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dR_nl_answer( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dF_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dChi_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dChi_NL_B_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d2OverlapdXi_1dGradChi_answer( overlap_answer.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d2OverlapddXddX_answer( overlap_answer.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d2OverlapddXdR_nl_answer( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    floatMatrix d2OverlapddXdF_answer( overlap_answer.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapddXdChi_answer( overlap_answer.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapddXdChi_NL_B_answer( overlap_answer.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d2OverlapddXdGradChi_answer( overlap_answer.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    floatVector d2OverlapdR_nldR_nl_answer( overlap_answer.size( ), 0 );

    floatMatrix d2OverlapdR_nldF_answer( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix d2OverlapdR_nldChi_answer( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix d2OverlapdR_nldChi_NL_B_answer( overlap_answer.size( ), floatVector( chi_nl_basis.size( ), 0 ) );

    floatMatrix d2OverlapdR_nldGradChi_answer( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdFdF_answer( overlap_answer.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    floatMatrix d2OverlapdFdChi_answer( overlap_answer.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdFdChi_NL_B_answer( overlap_answer.size( ), floatVector( F.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d2OverlapdFdGradChi_answer( overlap_answer.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdChidChi_answer( overlap_answer.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d2OverlapdChidChi_NL_B_answer( overlap_answer.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d2OverlapdChidGradChi_answer( overlap_answer.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdChi_NL_BdChi_NL_B_answer( overlap_answer.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d2OverlapdChi_NL_BdGradChi_answer( overlap_answer.size( ), floatVector( chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d2OverlapdGradChidGradChi_answer( overlap_answer.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    floatVector overlap_3;

    floatMatrix dOverlapdXi_1_3, dOverlapddX_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdChi_NL_B_3, dOverlapdGradChi_3;

    floatVector dOverlapdR_nl_3;

    floatMatrix d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dChi_NL_B_3, d2OverlapdXi_1dGradChi_3,
                d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdChi_NL_B_3, d2OverlapddXdGradChi_3,
                d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldChi_NL_B_3, d2OverlapdR_nldGradChi_3,
                d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdChi_NL_B_3, d2OverlapdFdGradChi_3,
                d2OverlapdChidChi_3, d2OverlapdChidChi_NL_B_3, d2OverlapdChidGradChi_3,
                d2OverlapdChi_NL_BdChi_NL_B_3, d2OverlapdChi_NL_BdGradChi_3,
                d2OverlapdGradChidGradChi_3;

    floatVector d2OverlapdR_nldR_nl_3;

    floatMatrix d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_NL_B, d3OverlapdXi_1dXi_1dGradChi,
                d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_NL_B, d3OverlapdXi_1ddXdGradChi,
                d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_NL_B, d3OverlapdXi_1dR_nldGradChi,
                d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_NL_B, d3OverlapdXi_1dFdGradChi,
                d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_NL_B, d3OverlapdXi_1dChidGradChi,
                d3OverlapdXi_1dChi_NL_BdChi_NL_B, d3OverlapdXi_1dChi_NL_BdGradChi,
                d3OverlapdXi_1dGradChidGradChi,
                d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdChi_NL_B, d3OverlapddXddXdGradChi,
                d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_NL_B, d3OverlapddXdR_nldGradChi,
                d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdChi_NL_B, d3OverlapddXdFdGradChi,
                d3OverlapddXdChidChi, d3OverlapddXdChidChi_NL_B, d3OverlapddXdChidGradChi,
                d3OverlapddXdChi_NL_BdChi_NL_B, d3OverlapddXdChi_NL_BdGradChi,
                d3OverlapddXdGradChidGradChi,
                d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_NL_B, d3OverlapdR_nldR_nldGradChi,
                d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_NL_B, d3OverlapdR_nldFdGradChi,
                d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_NL_B, d3OverlapdR_nldChidGradChi,
                d3OverlapdR_nldChi_NL_BdChi_NL_B, d3OverlapdR_nldChi_NL_BdGradChi,
                d3OverlapdR_nldGradChidGradChi,
                d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdChi_NL_B, d3OverlapdFdFdGradChi,
                d3OverlapdFdChidChi, d3OverlapdFdChidChi_NL_B, d3OverlapdFdChidGradChi,
                d3OverlapdFdChi_NL_BdChi_NL_B, d3OverlapdFdChi_NL_BdGradChi,
                d3OverlapdFdGradChidGradChi,
                d3OverlapdChidChidChi, d3OverlapdChidChidChi_NL_B, d3OverlapdChidChidGradChi,
                d3OverlapdChidChi_NL_BdChi_NL_B, d3OverlapdChidChi_NL_BdGradChi,
                d3OverlapdChidGradChidGradChi,
                d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B, d3OverlapdChi_NL_BdChi_NL_BdGradChi,
                d3OverlapdChi_NL_BdGradChidGradChi,
                d3OverlapdGradChidGradChidGradChi;

    floatVector d3OverlapdR_nldR_nldR_nl;

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap_3,
                                                              dOverlapdXi_1_3, dOverlapddX_3, dOverlapdR_nl_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdChi_NL_B_3, dOverlapdGradChi_3,
                                                              d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dChi_NL_B_3, d2OverlapdXi_1dGradChi_3,
                                                              d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdChi_NL_B_3, d2OverlapddXdGradChi_3,
                                                              d2OverlapdR_nldR_nl_3, d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldChi_NL_B_3, d2OverlapdR_nldGradChi_3,
                                                              d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdChi_NL_B_3, d2OverlapdFdGradChi_3,
                                                              d2OverlapdChidChi_3, d2OverlapdChidChi_NL_B_3, d2OverlapdChidGradChi_3,
                                                              d2OverlapdChi_NL_BdChi_NL_B_3, d2OverlapdChi_NL_BdGradChi_3,
                                                              d2OverlapdGradChidGradChi_3,
                                                              d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_NL_B, d3OverlapdXi_1dXi_1dGradChi,
                                                              d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_NL_B, d3OverlapdXi_1ddXdGradChi,
                                                              d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_NL_B, d3OverlapdXi_1dR_nldGradChi,
                                                              d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_NL_B, d3OverlapdXi_1dFdGradChi,
                                                              d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_NL_B, d3OverlapdXi_1dChidGradChi,
                                                              d3OverlapdXi_1dChi_NL_BdChi_NL_B, d3OverlapdXi_1dChi_NL_BdGradChi,
                                                              d3OverlapdXi_1dGradChidGradChi,
                                                              d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdChi_NL_B, d3OverlapddXddXdGradChi,
                                                              d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_NL_B, d3OverlapddXdR_nldGradChi,
                                                              d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdChi_NL_B, d3OverlapddXdFdGradChi,
                                                              d3OverlapddXdChidChi, d3OverlapddXdChidChi_NL_B, d3OverlapddXdChidGradChi,
                                                              d3OverlapddXdChi_NL_BdChi_NL_B, d3OverlapddXdChi_NL_BdGradChi,
                                                              d3OverlapddXdGradChidGradChi,
                                                              d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_NL_B, d3OverlapdR_nldR_nldGradChi,
                                                              d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_NL_B, d3OverlapdR_nldFdGradChi,
                                                              d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_NL_B, d3OverlapdR_nldChidGradChi,
                                                              d3OverlapdR_nldChi_NL_BdChi_NL_B, d3OverlapdR_nldChi_NL_BdGradChi,
                                                              d3OverlapdR_nldGradChidGradChi,
                                                              d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdChi_NL_B, d3OverlapdFdFdGradChi,
                                                              d3OverlapdFdChidChi, d3OverlapdFdChidChi_NL_B, d3OverlapdFdChidGradChi,
                                                              d3OverlapdFdChi_NL_BdChi_NL_B, d3OverlapdFdChi_NL_BdGradChi,
                                                              d3OverlapdFdGradChidGradChi,
                                                              d3OverlapdChidChidChi, d3OverlapdChidChidChi_NL_B, d3OverlapdChidChidGradChi,
                                                              d3OverlapdChidChi_NL_BdChi_NL_B, d3OverlapdChidChi_NL_BdGradChi,
                                                              d3OverlapdChidGradChidGradChi,
                                                              d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B, d3OverlapdChi_NL_BdChi_NL_BdGradChi,
                                                              d3OverlapdChi_NL_BdGradChidGradChi,
                                                              d3OverlapdGradChidGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_3, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_3, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_3, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_3, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_3, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_3, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_NL_B_3, dOverlapdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_3, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX_3, d2OverlapddXddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl_3, d2OverlapddXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF_3, d2OverlapddXdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_3, d2OverlapddXdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_NL_B_3, d2OverlapddXdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi_3, d2OverlapddXdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl_3, d2OverlapdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF_3, d2OverlapdR_nldF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_3, d2OverlapdR_nldChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_NL_B_3, d2OverlapdR_nldChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi_3, d2OverlapdR_nldGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF_3, d2OverlapdFdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_3, d2OverlapdFdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_NL_B_3, d2OverlapdFdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi_3, d2OverlapdFdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_3, d2OverlapdChidChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_NL_B_3, d2OverlapdChidChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi_3, d2OverlapdChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdChi_NL_B_3, d2OverlapdChi_NL_BdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdGradChi_3, d2OverlapdChi_NL_BdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi_3, d2OverlapdGradChidGradChi ) );

    floatMatrix d3OverlapdXi_1dXi_1dXi_1_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1ddX_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dChi_NL_B_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dXi_1dGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXddX_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdChi_NL_B_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1ddXdGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldR_nl_answer( Xi_1.size( ), floatVector( Xi_1.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldChi_NL_B_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dR_nldGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdF_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdChi_NL_B_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dFdGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChidChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChidChi_NL_B_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChidGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChi_NL_BdChi_NL_B_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dChi_NL_BdGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdXi_1dGradChidGradChi_answer( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXddXddX_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdR_nl_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdF_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdChi_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdChi_NL_B_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapddXddXdGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldR_nl_answer( Xi_1.size( ), floatVector( dX.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldF_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldChi_NL_B_answer( Xi_1.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapddXdR_nldGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdF_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdChi_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdChi_NL_B_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapddXdFdGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdChidChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapddXdChidChi_NL_B_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapddXdChidGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdChi_NL_BdChi_NL_B_answer( Xi_1.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapddXdChi_NL_BdGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapddXdGradChidGradChi_answer( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatVector d3OverlapdR_nldR_nldR_nl_answer( Xi_1.size( ), 0 );

    floatMatrix d3OverlapdR_nldR_nldF_answer( Xi_1.size( ), floatVector( F.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldR_nldChi_answer( Xi_1.size( ), floatVector( chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldR_nldChi_NL_B_answer( Xi_1.size( ), floatVector( chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldR_nldGradChi_answer( Xi_1.size( ), floatVector( gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdF_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdChi_NL_B_answer( Xi_1.size( ), floatVector( F.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldFdGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChidChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChidChi_NL_B_answer( Xi_1.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChidGradChi_answer( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChi_NL_BdChi_NL_B_answer( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldChi_NL_BdGradChi_answer( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdR_nldGradChidGradChi_answer( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdF_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * F.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdChi_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdChi_NL_B_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdFdFdGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdFdChidChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdFdChidChi_NL_B_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdFdChidGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdFdChi_NL_BdChi_NL_B_answer( Xi_1.size( ), floatVector( F.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdFdChi_NL_BdGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdFdGradChidGradChi_answer( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdChidChidChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi.size( ), 0 ) );

    floatMatrix d3OverlapdChidChidChi_NL_B_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdChidChidGradChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdChidChi_NL_BdChi_NL_B_answer( Xi_1.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdChidChi_NL_BdGradChi_answer( Xi_1.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdChidGradChidGradChi_answer( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B_answer( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    floatMatrix d3OverlapdChi_NL_BdChi_NL_BdGradChi_answer( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdChi_NL_BdGradChidGradChi_answer( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    floatMatrix d3OverlapdGradChidGradChidGradChi_answer( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    // Tests of the gradients for the non-overlapped case. Everything should be zero!

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){

        floatVector delta( Xi_1.size( ), 0 );

        delta[ i ] += eps * std::fabs( Xi_1[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdXi_1_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,     d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp,  d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dXi_1_answer[ j ][ Xi_1.size( ) * Xi_1.size( ) * k + Xi_1.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1dXi_1_answer ) );

    for ( unsigned int i = 0; i < dX.size( ); i++ ){

        floatVector delta( dX.size( ), 0 );

        delta[ i ] += eps * std::fabs( dX[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapddX_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1ddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1ddX_answer[ j ][ Xi_1.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1ddX, d2OverlapdXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX, d2OverlapddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXddX, d3OverlapddXddXddX_answer ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( R_nl ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdR_nl_answer[ j ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dR_nl_answer[ j ][ k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdR_nl_answer[ j ][ k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldR_nl_answer[ j + k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dR_nl_answer[ j ][ Xi_1.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldR_nl_answer[ j ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdR_nl, d3OverlapddXddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldR_nl_answer ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){

        floatVector delta( F.size( ), 0 );

        delta[ i ] += eps * std::fabs( F[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdF_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dF_answer[ j ][ Xi_1.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdF, d3OverlapddXddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldF, d3OverlapddXdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdF, d3OverlapddXdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdF, d3OverlapdR_nldFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdF, d3OverlapdFdFdF_answer ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi, d3OverlapddXdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi, d3OverlapddXdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi, d3OverlapdFdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi, d3OverlapdFdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi, d3OverlapdChidChidChi_answer ) );

    for ( unsigned int i = 0; i < chi_nl_basis.size( ); i++ ){

        floatVector delta( chi_nl_basis.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi_nl_basis[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis + delta, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis - delta, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_NL_B_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                d2OverlapdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * k + i ] = ( dOverlapdChi_NL_Bp[ j ][ k ] - dOverlapdChi_NL_Bm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );


        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_NL_B_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_NL_B_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_NL_B_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_NL_B_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_NL_B_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdXi_1dChi_NL_BdChi_NL_B_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChi_NL_Bp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChi_NL_Bm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_NL_B_answer[ j ][ dX.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_NL_B_answer[ j ][ F.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapddXdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdChi_NL_Bp[ j ][ F.size( ) * k + l ] - d2OverlapddXdChi_NL_Bm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_NL_B_answer[ j ][ F.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdR_nldChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldChi_NL_Bp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChi_NL_Bm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_NL_B_answer[ j ][ F.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ chi.size( ) * k + l ] - d2OverlapdFdChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdFdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdFdChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdFdChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdChidChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdChidChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdChi_NL_BdChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdChi_NL_BdChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_NL_B, dOverlapdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi_NL_B, d2OverlapdXi_1dChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_NL_B, d2OverlapddXdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_NL_B, d2OverlapdR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_NL_B, d2OverlapdFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_NL_B, d2OverlapdChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdChi_NL_B, d2OverlapdChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi_NL_B, d3OverlapdXi_1dXi_1dChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi_NL_B, d3OverlapdXi_1ddXdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi_NL_B, d3OverlapdXi_1dR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi_NL_B, d3OverlapdXi_1dFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi_NL_B, d3OverlapdXi_1dChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChi_NL_BdChi_NL_B, d3OverlapdXi_1dChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi_NL_B, d3OverlapddXddXdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi_NL_B, d3OverlapddXdR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi_NL_B, d3OverlapddXdFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi_NL_B, d3OverlapddXdChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChi_NL_BdChi_NL_B, d3OverlapddXdChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi_NL_B, d3OverlapdR_nldR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi_NL_B, d3OverlapdR_nldFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi_NL_B, d3OverlapdR_nldChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChi_NL_BdChi_NL_B, d3OverlapdR_nldChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi_NL_B, d3OverlapdFdFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi_NL_B, d3OverlapdFdChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChi_NL_BdChi_NL_B, d3OverlapdFdChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi_NL_B, d3OverlapdChidChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChi_NL_BdChi_NL_B, d3OverlapdChidChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B, d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B_answer ) );

    for ( unsigned int i = 0; i < gradChi.size( ); i++ ){

        floatVector delta( gradChi.size( ), 0 );

        delta[ i ] += eps * std::fabs( gradChi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi + delta, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi - delta, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdGradChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                d2OverlapdChi_NL_BdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdChi_NL_Bp[ j ][ k ] - dOverlapdChi_NL_Bm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < gradChi.size( ); k++ ){

                d2OverlapdGradChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdGradChip[ j ][ k ] - dOverlapdGradChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );


        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dGradChi_answer[ j ][ Xi_1.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdGradChi_answer[ j ][ dX.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ chi.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdXi_1dChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdXi_1dChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdXi_1dGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdXi_1dGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdGradChi_answer[ j ][ dX.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ chi.size( ) * k + l ] - d2OverlapddXdChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapddXdChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapddXdChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapddXdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapddXdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ chi.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdR_nldChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdR_nldChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdR_nldGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdR_nldGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ chi.size( ) * k + l ] - d2OverlapdFdChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdFdChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdFdChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdFdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdFdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdChidChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdChidGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdChidGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdChi_NL_BdChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChi_NL_BdChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdChi_NL_BdChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdChi_NL_BdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChi_NL_BdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdChi_NL_BdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < gradChi.size( ); k++ ){

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdGradChidGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdGradChidGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdGradChidGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi, dOverlapdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dGradChi, d2OverlapdXi_1dGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi, d2OverlapddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi, d2OverlapdR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi, d2OverlapdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi, d2OverlapdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdGradChi, d2OverlapdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi, d2OverlapdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dGradChi,        d3OverlapdXi_1dXi_1dGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdGradChi,          d3OverlapdXi_1ddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldGradChi,        d3OverlapdXi_1dR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdGradChi,           d3OverlapdXi_1dFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidGradChi,         d3OverlapdXi_1dChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChi_NL_BdGradChi,    d3OverlapdXi_1dChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dGradChidGradChi,     d3OverlapdXi_1dGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdGradChi,            d3OverlapdXi_1ddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldGradChi,          d3OverlapddXdR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdGradChi,             d3OverlapddXdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidGradChi,           d3OverlapddXdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChi_NL_BdGradChi,      d3OverlapddXdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdGradChidGradChi,       d3OverlapddXdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldGradChi,        d3OverlapdR_nldR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdGradChi,           d3OverlapdR_nldFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidGradChi,         d3OverlapdR_nldChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChi_NL_BdGradChi,    d3OverlapdR_nldChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldGradChidGradChi,     d3OverlapdR_nldGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdGradChi,              d3OverlapdFdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidGradChi,            d3OverlapdFdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChi_NL_BdGradChi,       d3OverlapdFdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdGradChidGradChi,        d3OverlapdFdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidGradChi,           d3OverlapdChidChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChi_NL_BdGradChi,      d3OverlapdChidChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidGradChidGradChi,       d3OverlapdChidGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChi_NL_BdChi_NL_BdGradChi, d3OverlapdChi_NL_BdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChi_NL_BdGradChidGradChi,  d3OverlapdChi_NL_BdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdGradChidGradChidGradChi,   d3OverlapdGradChidGradChidGradChi_answer ) );

    // Tests of the gradients of the overlapped case. Gradients will be non-zero

    R_nl = 2.5;

    overlap_answer = { -0.0149698, -0.00571608, 0.00344577 };

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap,
                                                              dOverlapdXi_1, dOverlapddX, dOverlapdR_nl, dOverlapdF, dOverlapdChi, dOverlapdChi_NL_B, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap, overlap_answer ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap_2,
                                                              dOverlapdXi_1_2, dOverlapddX_2, dOverlapdR_nl_2, dOverlapdF_2, dOverlapdChi_2, dOverlapdChi_NL_B_2, dOverlapdGradChi_2,
                                                              d2OverlapdXi_1dXi_1, d2OverlapdXi_1ddX, d2OverlapdXi_1dR_nl, d2OverlapdXi_1dF, d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_NL_B, d2OverlapdXi_1dGradChi,
                                                              d2OverlapddXddX, d2OverlapddXdR_nl, d2OverlapddXdF, d2OverlapddXdChi, d2OverlapddXdChi_NL_B, d2OverlapddXdGradChi,
                                                              d2OverlapdR_nldR_nl, d2OverlapdR_nldF, d2OverlapdR_nldChi, d2OverlapdR_nldChi_NL_B, d2OverlapdR_nldGradChi,
                                                              d2OverlapdFdF, d2OverlapdFdChi, d2OverlapdFdChi_NL_B, d2OverlapdFdGradChi,
                                                              d2OverlapdChidChi, d2OverlapdChidChi_NL_B, d2OverlapdChidGradChi,
                                                              d2OverlapdChi_NL_BdChi_NL_B, d2OverlapdChi_NL_BdGradChi,
                                                              d2OverlapdGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_2, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_2, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_2, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_2, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_2, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_2, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_NL_B_2, dOverlapdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_2, dOverlapdGradChi ) );

    dOverlapdXi_1_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    dOverlapddX_answer = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    dOverlapdR_nl_answer = floatVector( overlap_answer.size( ), 0 );

    dOverlapdF_answer = floatMatrix( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    dOverlapdChi_answer = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    dOverlapdChi_NL_B_answer = floatMatrix( overlap_answer.size( ), floatVector( chi_nl_basis.size( ), 0 ) );

    dOverlapdGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    d2OverlapdXi_1dXi_1_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d2OverlapdXi_1ddX_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    d2OverlapdXi_1dR_nl_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ), 0 ) );

    d2OverlapdXi_1dF_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    d2OverlapdXi_1dChi_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    d2OverlapdXi_1dChi_NL_B_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ), 0 ) );

    d2OverlapdXi_1dGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    d2OverlapddXddX_answer = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    d2OverlapddXdR_nl_answer = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ), 0 ) );

    d2OverlapddXdF_answer = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    d2OverlapddXdChi_answer = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    d2OverlapddXdChi_NL_B_answer = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ), 0 ) );

    d2OverlapddXdGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    d2OverlapdR_nldR_nl_answer = floatVector( overlap_answer.size( ), 0 );

    d2OverlapdR_nldF_answer = floatMatrix( overlap_answer.size( ), floatVector( F.size( ), 0 ) );

    d2OverlapdR_nldChi_answer = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ), 0 ) );

    d2OverlapdR_nldChi_NL_B_answer = floatMatrix( overlap_answer.size( ), floatVector( chi_nl_basis.size( ), 0 ) );

    d2OverlapdR_nldGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( gradChi.size( ), 0 ) );

    d2OverlapdFdF_answer = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    d2OverlapdFdChi_answer = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    d2OverlapdFdChi_NL_B_answer = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * chi_nl_basis.size( ), 0 ) );

    d2OverlapdFdGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    d2OverlapdChidChi_answer = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    d2OverlapdChidChi_NL_B_answer = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ), 0 ) );

    d2OverlapdChidGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    d2OverlapdChi_NL_BdChi_NL_B_answer = floatMatrix( overlap_answer.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    d2OverlapdChi_NL_BdGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    d2OverlapdGradChidGradChi_answer = floatMatrix( overlap_answer.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlap_3,
                                                              dOverlapdXi_1_3, dOverlapddX_3, dOverlapdR_nl_3, dOverlapdF_3, dOverlapdChi_3, dOverlapdChi_NL_B_3, dOverlapdGradChi_3,
                                                              d2OverlapdXi_1dXi_1_3, d2OverlapdXi_1ddX_3, d2OverlapdXi_1dR_nl_3, d2OverlapdXi_1dF_3, d2OverlapdXi_1dChi_3, d2OverlapdXi_1dChi_NL_B_3, d2OverlapdXi_1dGradChi_3,
                                                              d2OverlapddXddX_3, d2OverlapddXdR_nl_3, d2OverlapddXdF_3, d2OverlapddXdChi_3, d2OverlapddXdChi_NL_B_3, d2OverlapddXdGradChi_3,
                                                              d2OverlapdR_nldR_nl_3, d2OverlapdR_nldF_3, d2OverlapdR_nldChi_3, d2OverlapdR_nldChi_NL_B_3, d2OverlapdR_nldGradChi_3,
                                                              d2OverlapdFdF_3, d2OverlapdFdChi_3, d2OverlapdFdChi_NL_B_3, d2OverlapdFdGradChi_3,
                                                              d2OverlapdChidChi_3, d2OverlapdChidChi_NL_B_3, d2OverlapdChidGradChi_3,
                                                              d2OverlapdChi_NL_BdChi_NL_B_3, d2OverlapdChi_NL_BdGradChi_3,
                                                              d2OverlapdGradChidGradChi_3,
                                                              d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_NL_B, d3OverlapdXi_1dXi_1dGradChi,
                                                              d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_NL_B, d3OverlapdXi_1ddXdGradChi,
                                                              d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_NL_B, d3OverlapdXi_1dR_nldGradChi,
                                                              d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_NL_B, d3OverlapdXi_1dFdGradChi,
                                                              d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_NL_B, d3OverlapdXi_1dChidGradChi,
                                                              d3OverlapdXi_1dChi_NL_BdChi_NL_B, d3OverlapdXi_1dChi_NL_BdGradChi,
                                                              d3OverlapdXi_1dGradChidGradChi,
                                                              d3OverlapddXddXddX, d3OverlapddXddXdR_nl, d3OverlapddXddXdF, d3OverlapddXddXdChi, d3OverlapddXddXdChi_NL_B, d3OverlapddXddXdGradChi,
                                                              d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldF, d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_NL_B, d3OverlapddXdR_nldGradChi,
                                                              d3OverlapddXdFdF, d3OverlapddXdFdChi, d3OverlapddXdFdChi_NL_B, d3OverlapddXdFdGradChi,
                                                              d3OverlapddXdChidChi, d3OverlapddXdChidChi_NL_B, d3OverlapddXdChidGradChi,
                                                              d3OverlapddXdChi_NL_BdChi_NL_B, d3OverlapddXdChi_NL_BdGradChi,
                                                              d3OverlapddXdGradChidGradChi,
                                                              d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_NL_B, d3OverlapdR_nldR_nldGradChi,
                                                              d3OverlapdR_nldFdF, d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_NL_B, d3OverlapdR_nldFdGradChi,
                                                              d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_NL_B, d3OverlapdR_nldChidGradChi,
                                                              d3OverlapdR_nldChi_NL_BdChi_NL_B, d3OverlapdR_nldChi_NL_BdGradChi,
                                                              d3OverlapdR_nldGradChidGradChi,
                                                              d3OverlapdFdFdF, d3OverlapdFdFdChi, d3OverlapdFdFdChi_NL_B, d3OverlapdFdFdGradChi,
                                                              d3OverlapdFdChidChi, d3OverlapdFdChidChi_NL_B, d3OverlapdFdChidGradChi,
                                                              d3OverlapdFdChi_NL_BdChi_NL_B, d3OverlapdFdChi_NL_BdGradChi,
                                                              d3OverlapdFdGradChidGradChi,
                                                              d3OverlapdChidChidChi, d3OverlapdChidChidChi_NL_B, d3OverlapdChidChidGradChi,
                                                              d3OverlapdChidChi_NL_BdChi_NL_B, d3OverlapdChidChi_NL_BdGradChi,
                                                              d3OverlapdChidGradChidGradChi,
                                                              d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B, d3OverlapdChi_NL_BdChi_NL_BdGradChi,
                                                              d3OverlapdChi_NL_BdGradChidGradChi,
                                                              d3OverlapdGradChidGradChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( overlap_3, overlap_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1_3, dOverlapdXi_1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX_3, dOverlapddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl_3, dOverlapdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF_3, dOverlapdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_3, dOverlapdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_NL_B_3, dOverlapdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi_3, dOverlapdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX_3, d2OverlapddXddX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl_3, d2OverlapddXdR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF_3, d2OverlapddXdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_3, d2OverlapddXdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_NL_B_3, d2OverlapddXdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi_3, d2OverlapddXdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl_3, d2OverlapdR_nldR_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF_3, d2OverlapdR_nldF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_3, d2OverlapdR_nldChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_NL_B_3, d2OverlapdR_nldChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi_3, d2OverlapdR_nldGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF_3, d2OverlapdFdF ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_3, d2OverlapdFdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_NL_B_3, d2OverlapdFdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi_3, d2OverlapdFdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_3, d2OverlapdChidChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_NL_B_3, d2OverlapdChidChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi_3, d2OverlapdChidGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdChi_NL_B_3, d2OverlapdChi_NL_BdChi_NL_B ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdGradChi_3, d2OverlapdChi_NL_BdGradChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi_3, d2OverlapdGradChidGradChi ) );

    d3OverlapdXi_1dXi_1dXi_1_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dXi_1ddX_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dR_nl_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dF_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdXi_1dXi_1dGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * Xi_1.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1ddXddX_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1ddXdR_nl_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ), 0 ) );

    d3OverlapdXi_1ddXdF_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1ddXdChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1ddXdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdXi_1ddXdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dR_nldR_nl_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ), 0 ) );

    d3OverlapdXi_1dR_nldF_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dR_nldChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dR_nldChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdXi_1dR_nldGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dFdF_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapdXi_1dFdChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dFdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdXi_1dFdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dChidChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdXi_1dChidChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdXi_1dChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dChi_NL_BdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdXi_1dChi_NL_BdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdXi_1dGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( Xi_1.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXddXddX_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * dX.size( ), 0 ) );

    d3OverlapddXddXdR_nl_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ), 0 ) );

    d3OverlapddXddXdF_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * F.size( ), 0 ) );

    d3OverlapddXddXdChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi.size( ), 0 ) );

    d3OverlapddXddXdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapddXddXdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * dX.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdR_nldR_nl_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ), 0 ) );

    d3OverlapddXdR_nldF_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ), 0 ) );

    d3OverlapddXdR_nldChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdR_nldChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapddXdR_nldGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdFdF_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapddXdFdChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdFdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapddXdFdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdChidChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapddXdChidChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapddXdChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdChi_NL_BdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapddXdChi_NL_BdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    d3OverlapddXdGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( dX.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdR_nldR_nldR_nl_answer = floatVector( Xi_1.size( ), 0 );

    d3OverlapdR_nldR_nldF_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ), 0 ) );

    d3OverlapdR_nldR_nldChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ), 0 ) );

    d3OverlapdR_nldR_nldChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( chi_nl_basis.size( ), 0 ) );

    d3OverlapdR_nldR_nldGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( gradChi.size( ), 0 ) );

    d3OverlapdR_nldFdF_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ), 0 ) );

    d3OverlapdR_nldFdChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ), 0 ) );

    d3OverlapdR_nldFdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdR_nldFdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdR_nldChidChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdR_nldChidChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdR_nldChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdR_nldChi_NL_BdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdR_nldChi_NL_BdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdR_nldGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdFdFdF_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * F.size( ), 0 ) );

    d3OverlapdFdFdChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi.size( ), 0 ) );

    d3OverlapdFdFdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdFdFdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * F.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdFdChidChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdFdChidChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdFdChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdFdChi_NL_BdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdFdChi_NL_BdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdFdGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( F.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdChidChidChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi.size( ), 0 ) );

    d3OverlapdChidChidChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdChidChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdChidChi_NL_BdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdChidChi_NL_BdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdChidGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B_answer = floatMatrix( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ) * chi_nl_basis.size( ), 0 ) );

    d3OverlapdChi_NL_BdChi_NL_BdGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * chi_nl_basis.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdChi_NL_BdGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( chi_nl_basis.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    d3OverlapdGradChidGradChidGradChi_answer = floatMatrix( Xi_1.size( ), floatVector( gradChi.size( ) * gradChi.size( ) * gradChi.size( ), 0 ) );

    // Tests of the gradients for the non-overlapped case. Everything should be zero!

    for ( unsigned int i = 0; i < Xi_1.size( ); i++ ){

        floatVector delta( Xi_1.size( ), 0 );

        delta[ i ] += eps * std::fabs( Xi_1[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdXi_1_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dXi_1_answer[ j ][ Xi_1.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 + delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,     d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp,  d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1 - delta, dX, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dXi_1_answer[ j ][ Xi_1.size( ) * Xi_1.size( ) * k + Xi_1.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dXi_1, d3OverlapdXi_1dXi_1dXi_1_answer ) );

    for ( unsigned int i = 0; i < dX.size( ); i++ ){

        floatVector delta( dX.size( ), 0 );

        delta[ i ] += eps * std::fabs( dX[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapddX_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1ddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXddX_answer[ j ][ dX.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX + delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX - delta, R_nl, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1ddX_answer[ j ][ Xi_1.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXddX_answer[ j ][ dX.size( ) * dX.size( ) * k + dX.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1ddX, d2OverlapdXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX, d2OverlapddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1ddX, d3OverlapdXi_1dXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXddX, d3OverlapdXi_1ddXddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXddX, d3OverlapddXddXddX_answer ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( R_nl ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdR_nl_answer[ j ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dR_nl_answer[ j ][ k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdR_nl_answer[ j ][ k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldR_nl_answer[ j + k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl + delta, F, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl - delta, F, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dR_nl_answer[ j ][ Xi_1.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdR_nl_answer[ j ][ dX.size( ) * k + l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldR_nl_answer[ j ][ k + l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldR_nl_answer[ j ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dR_nl, d3OverlapdXi_1dXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdR_nl, d3OverlapdXi_1ddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldR_nl, d3OverlapdXi_1dR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdR_nl, d3OverlapddXddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldR_nl, d3OverlapddXdR_nldR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldR_nl, d3OverlapdR_nldR_nldR_nl_answer ) );

    for ( unsigned int i = 0; i < F.size( ); i++ ){

        floatVector delta( F.size( ), 0 );

        delta[ i ] += eps * std::fabs( F[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdF_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdF_answer[ j ][ F.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F + delta, chi, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F - delta, chi, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dF_answer[ j ][ Xi_1.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdF_answer[ j ][ dX.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldF_answer[ j ][ F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdF_answer[ j ][ F.size( ) * F.size( ) * k + F.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dF, d3OverlapdXi_1dXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdF, d3OverlapdXi_1ddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldF, d3OverlapdXi_1dR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdF, d3OverlapdXi_1dFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdF, d3OverlapddXddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldF, d3OverlapddXdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdF, d3OverlapddXdFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldF, d3OverlapdR_nldR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdF, d3OverlapdR_nldFdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdF, d3OverlapdFdFdF_answer ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, chi_nl_basis, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, chi_nl_basis, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi + delta, chi_nl_basis, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi - delta, chi_nl_basis, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ F.size( ) * k + l ] - d2OverlapdFdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi, d3OverlapdXi_1dXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi, d3OverlapdXi_1ddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi, d3OverlapdXi_1dR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi, d3OverlapdXi_1dFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi, d3OverlapdXi_1dChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi, d3OverlapddXdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi, d3OverlapddXdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi, d3OverlapddXdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi, d3OverlapdR_nldR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi, d3OverlapdR_nldFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi, d3OverlapdR_nldChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi, d3OverlapdFdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi, d3OverlapdFdChidChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi, d3OverlapdChidChidChi_answer ) );

    for ( unsigned int i = 0; i < chi_nl_basis.size( ); i++ ){

        floatVector delta( chi_nl_basis.size( ), 0 );

        delta[ i ] += eps * std::fabs( chi_nl_basis[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis + delta, gradChi, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis - delta, gradChi, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdChi_NL_B_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidChi_NL_B_answer[ j ][ chi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                d2OverlapdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * k + i ] = ( dOverlapdChi_NL_Bp[ j ][ k ] - dOverlapdChi_NL_Bm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis + delta, gradChi, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis - delta, gradChi, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );


        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dChi_NL_B_answer[ j ][ Xi_1.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdChi_NL_B_answer[ j ][ dX.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldChi_NL_B_answer[ j ][ chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdChi_NL_B_answer[ j ][ F.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidChi_NL_B_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdXi_1dChi_NL_BdChi_NL_B_answer[ j ][ chi.size( ) * chi.size( ) * k + chi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChi_NL_Bp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dChi_NL_Bm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdChi_NL_B_answer[ j ][ dX.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdChi_NL_B_answer[ j ][ F.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ F.size( ) * k + l ] - d2OverlapddXdChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapddXdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapddXdChi_NL_Bp[ j ][ F.size( ) * k + l ] - d2OverlapddXdChi_NL_Bm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdChi_NL_B_answer[ j ][ F.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdR_nldChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdR_nldChi_NL_Bp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldChi_NL_Bm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdChi_NL_B_answer[ j ][ F.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ chi.size( ) * k + l ] - d2OverlapdFdChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdFdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdFdChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdFdChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidChi_NL_B_answer[ j ][ chi.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdChidChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdChidChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B_answer[ j ][ chi_nl_basis.size( ) * chi_nl_basis.size( ) * k + chi_nl_basis.size( ) * l + i ]
                        = ( d2OverlapdChi_NL_BdChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdChi_NL_BdChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi_NL_B, dOverlapdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi_NL_B, d2OverlapdXi_1dChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi_NL_B, d2OverlapddXdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi_NL_B, d2OverlapdR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi_NL_B, d2OverlapdFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi_NL_B, d2OverlapdChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdChi_NL_B, d2OverlapdChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dChi_NL_B, d3OverlapdXi_1dXi_1dChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdChi_NL_B, d3OverlapdXi_1ddXdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldChi_NL_B, d3OverlapdXi_1dR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdChi_NL_B, d3OverlapdXi_1dFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidChi_NL_B, d3OverlapdXi_1dChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChi_NL_BdChi_NL_B, d3OverlapdXi_1dChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi, d3OverlapddXddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdChi_NL_B, d3OverlapddXddXdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldChi_NL_B, d3OverlapddXdR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdChi_NL_B, d3OverlapddXdFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidChi_NL_B, d3OverlapddXdChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChi_NL_BdChi_NL_B, d3OverlapddXdChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldChi_NL_B, d3OverlapdR_nldR_nldChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdChi_NL_B, d3OverlapdR_nldFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidChi_NL_B, d3OverlapdR_nldChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChi_NL_BdChi_NL_B, d3OverlapdR_nldChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdChi_NL_B, d3OverlapdFdFdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidChi_NL_B, d3OverlapdFdChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChi_NL_BdChi_NL_B, d3OverlapdFdChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidChi_NL_B, d3OverlapdChidChidChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChi_NL_BdChi_NL_B, d3OverlapdChidChi_NL_BdChi_NL_B_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B, d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B_answer ) );

    for ( unsigned int i = 0; i < gradChi.size( ); i++ ){

        floatVector delta( gradChi.size( ), 0 );

        delta[ i ] += eps * std::fabs( gradChi[ i ] ) + eps;

        floatVector overlapp, overlapm;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi + delta, overlapp ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi - delta, overlapm ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            dOverlapdGradChi_answer[ j ][ i ] += ( overlapp[ j ] - overlapm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dOverlapdXi_1p, dOverlapdXi_1m;

        floatMatrix dOverlapddXp, dOverlapddXm;

        floatVector dOverlapdR_nlp, dOverlapdR_nlm;

        floatMatrix dOverlapdFp, dOverlapdFm;

        floatMatrix dOverlapdChip, dOverlapdChim;

        floatMatrix dOverlapdChi_NL_Bp, dOverlapdChi_NL_Bm;

        floatMatrix dOverlapdGradChip, dOverlapdGradChim;

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp, dOverlapdGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim ) );

        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                d2OverlapdXi_1dGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdXi_1p[ j ][ k ] - dOverlapdXi_1m[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                d2OverlapddXdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapddXp[ j ][ k ] - dOverlapddXm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                d2OverlapdR_nldGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdR_nlp[ j + k ] - dOverlapdR_nlm[ j + k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                d2OverlapdFdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdFp[ j ][ k ] - dOverlapdFm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                d2OverlapdChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdChip[ j ][ k ] - dOverlapdChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                d2OverlapdChi_NL_BdGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdChi_NL_Bp[ j ][ k ] - dOverlapdChi_NL_Bm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < gradChi.size( ); k++ ){

                d2OverlapdGradChidGradChi_answer[ j ][ gradChi.size( ) * k + i ] = ( dOverlapdGradChip[ j ][ k ] - dOverlapdGradChim[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

        }

        floatMatrix d2OverlapdXi_1dXi_1p, d2OverlapdXi_1ddXp, d2OverlapdXi_1dR_nlp, d2OverlapdXi_1dFp, d2OverlapdXi_1dChip, d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                    d2OverlapddXddXp, d2OverlapddXdR_nlp, d2OverlapddXdFp, d2OverlapddXdChip, d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                    d2OverlapdR_nldFp, d2OverlapdR_nldChip, d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                    d2OverlapdFdFp, d2OverlapdFdChip, d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                    d2OverlapdChidChip, d2OverlapdChidChi_NL_Bp, d2OverlapdChidGradChip,
                    d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                    d2OverlapdGradChidGradChip;
    
        floatMatrix d2OverlapdXi_1dXi_1m, d2OverlapdXi_1ddXm, d2OverlapdXi_1dR_nlm, d2OverlapdXi_1dFm, d2OverlapdXi_1dChim, d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                    d2OverlapddXddXm, d2OverlapddXdR_nlm, d2OverlapddXdFm, d2OverlapddXdChim, d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                    d2OverlapdR_nldFm, d2OverlapdR_nldChim, d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                    d2OverlapdFdFm, d2OverlapdFdChim, d2OverlapdFdChi_NL_Bm, d2OverlapdFdGradChim,
                    d2OverlapdChidChim, d2OverlapdChidChi_NL_Bm, d2OverlapdChidGradChim,
                    d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                    d2OverlapdGradChidGradChim;
    
        floatVector d2OverlapdR_nldR_nlp, d2OverlapdR_nldR_nlm;
    
        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi + delta, overlapp,
                                                                  dOverlapdXi_1p, dOverlapddXp, dOverlapdR_nlp, dOverlapdFp, dOverlapdChip, dOverlapdChi_NL_Bp,  dOverlapdGradChip,
                                                                  d2OverlapdXi_1dXi_1p,         d2OverlapdXi_1ddXp,         d2OverlapdXi_1dR_nlp,  d2OverlapdXi_1dFp,        d2OverlapdXi_1dChip,    d2OverlapdXi_1dChi_NL_Bp, d2OverlapdXi_1dGradChip,
                                                                  d2OverlapddXddXp,             d2OverlapddXdR_nlp,         d2OverlapddXdFp,       d2OverlapddXdChip,        d2OverlapddXdChi_NL_Bp, d2OverlapddXdGradChip,
                                                                  d2OverlapdR_nldR_nlp,         d2OverlapdR_nldFp,          d2OverlapdR_nldChip,   d2OverlapdR_nldChi_NL_Bp, d2OverlapdR_nldGradChip,
                                                                  d2OverlapdFdFp,               d2OverlapdFdChip,           d2OverlapdFdChi_NL_Bp, d2OverlapdFdGradChip,
                                                                  d2OverlapdChidChip,           d2OverlapdChidChi_NL_Bp,    d2OverlapdChidGradChip,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bp, d2OverlapdChi_NL_BdGradChip,
                                                                  d2OverlapdGradChidGradChip ) );

        BOOST_CHECK( !tractionSeparation::computeParticleOverlap( Xi_1, dX, R_nl, F, chi, chi_nl_basis, gradChi - delta, overlapm,
                                                                  dOverlapdXi_1m, dOverlapddXm, dOverlapdR_nlm, dOverlapdFm, dOverlapdChim, dOverlapdChi_NL_Bm, dOverlapdGradChim,
                                                                  d2OverlapdXi_1dXi_1m,         d2OverlapdXi_1ddXm,         d2OverlapdXi_1dR_nlm,   d2OverlapdXi_1dFm,        d2OverlapdXi_1dChim,    d2OverlapdXi_1dChi_NL_Bm, d2OverlapdXi_1dGradChim,
                                                                  d2OverlapddXddXm,             d2OverlapddXdR_nlm,         d2OverlapddXdFm,        d2OverlapddXdChim,        d2OverlapddXdChi_NL_Bm, d2OverlapddXdGradChim,
                                                                  d2OverlapdR_nldR_nlm,         d2OverlapdR_nldFm,          d2OverlapdR_nldChim,    d2OverlapdR_nldChi_NL_Bm, d2OverlapdR_nldGradChim,
                                                                  d2OverlapdFdFm,               d2OverlapdFdChim,           d2OverlapdFdChi_NL_Bm,  d2OverlapdFdGradChim,
                                                                  d2OverlapdChidChim,           d2OverlapdChidChi_NL_Bm,    d2OverlapdChidGradChim,
                                                                  d2OverlapdChi_NL_BdChi_NL_Bm, d2OverlapdChi_NL_BdGradChim,
                                                                  d2OverlapdGradChidGradChim ) );


        for ( unsigned int j = 0; j < overlap_answer.size( ); j++ ){

            for ( unsigned int k = 0; k < Xi_1.size( ); k++ ){

                for ( unsigned int l = 0; l < Xi_1.size( ); l++ ){

                    d3OverlapdXi_1dXi_1dGradChi_answer[ j ][ Xi_1.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dXi_1p[ j ][ Xi_1.size( ) * k + l ] - d2OverlapdXi_1dXi_1m[ j ][ Xi_1.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapdXi_1ddXdGradChi_answer[ j ][ dX.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1ddXp[ j ][ dX.size( ) * k + l ] - d2OverlapdXi_1ddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdXi_1dR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dR_nlp[ j ][ k + l ] - d2OverlapdXi_1dR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdXi_1dFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dFp[ j ][ F.size( ) * k + l ] - d2OverlapdXi_1dFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdXi_1dChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChip[ j ][ chi.size( ) * k + l ] - d2OverlapdXi_1dChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdXi_1dChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdXi_1dChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdXi_1dGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdXi_1dGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdXi_1dGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < dX.size( ); k++ ){

                for ( unsigned int l = 0; l < dX.size( ); l++ ){

                    d3OverlapddXddXdGradChi_answer[ j ][ dX.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXddXp[ j ][ dX.size( ) * k + l ] - d2OverlapddXddXm[ j ][ dX.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapddXdR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdR_nlp[ j ][ k + l ] - d2OverlapddXdR_nlm[ j ][ k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapddXdFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdFp[ j ][ F.size( ) * k + l ] - d2OverlapddXdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapddXdChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdChip[ j ][ chi.size( ) * k + l ] - d2OverlapddXdChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapddXdChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapddXdChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapddXdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapddXdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapddXdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < 1; k++ ){

                for ( unsigned int l = 0; l < 1; l++ ){

                    d3OverlapdR_nldR_nldGradChi_answer[ j ][ gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldR_nlp[ j ] - d2OverlapdR_nldR_nlm[ j ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdR_nldFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldFp[ j ][ F.size( ) * k + l ] - d2OverlapdR_nldFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdR_nldChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChip[ j ][ chi.size( ) * k + l ] - d2OverlapdR_nldChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdR_nldChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdR_nldChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdR_nldGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdR_nldGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdR_nldGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < F.size( ); k++ ){

                for ( unsigned int l = 0; l < F.size( ); l++ ){

                    d3OverlapdFdFdGradChi_answer[ j ][ F.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdFp[ j ][ F.size( ) * k + l ] - d2OverlapdFdFm[ j ][ F.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdFdChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdChip[ j ][ chi.size( ) * k + l ] - d2OverlapdFdChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdFdChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdFdChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdFdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdFdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdFdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi.size( ); k++ ){

                for ( unsigned int l = 0; l < chi.size( ); l++ ){

                    d3OverlapdChidChidGradChi_answer[ j ][ chi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidChip[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChim[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdChidChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidChi_NL_Bp[ j ][ chi.size( ) * k + l ] - d2OverlapdChidChi_NL_Bm[ j ][ chi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdChidGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChidGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdChidGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < chi_nl_basis.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl_basis.size( ); l++ ){

                    d3OverlapdChi_NL_BdChi_NL_BdGradChi_answer[ j ][ chi_nl_basis.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChi_NL_BdChi_NL_Bp[ j ][ chi_nl_basis.size( ) * k + l ] - d2OverlapdChi_NL_BdChi_NL_Bm[ j ][ chi_nl_basis.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdChi_NL_BdGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdChi_NL_BdGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdChi_NL_BdGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

            for ( unsigned int k = 0; k < gradChi.size( ); k++ ){

                for ( unsigned int l = 0; l < gradChi.size( ); l++ ){

                    d3OverlapdGradChidGradChidGradChi_answer[ j ][ gradChi.size( ) * gradChi.size( ) * k + gradChi.size( ) * l + i ]
                        = ( d2OverlapdGradChidGradChip[ j ][ gradChi.size( ) * k + l ] - d2OverlapdGradChidGradChim[ j ][ gradChi.size( ) * k + l ] ) / ( 2 * delta[ i ] );

                }

            }

        }
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdGradChi, dOverlapdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dGradChi, d2OverlapdXi_1dGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdGradChi, d2OverlapddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldGradChi, d2OverlapdR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdGradChi, d2OverlapdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidGradChi, d2OverlapdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChi_NL_BdGradChi, d2OverlapdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdGradChidGradChi, d2OverlapdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dXi_1dGradChi,        d3OverlapdXi_1dXi_1dGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1ddXdGradChi,          d3OverlapdXi_1ddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dR_nldGradChi,        d3OverlapdXi_1dR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dFdGradChi,           d3OverlapdXi_1dFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChidGradChi,         d3OverlapdXi_1dChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dChi_NL_BdGradChi,    d3OverlapdXi_1dChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdXi_1dGradChidGradChi,     d3OverlapdXi_1dGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXddXdGradChi,            d3OverlapddXddXdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdR_nldGradChi,          d3OverlapddXdR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdFdGradChi,             d3OverlapddXdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChidGradChi,           d3OverlapddXdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdChi_NL_BdGradChi,      d3OverlapddXdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapddXdGradChidGradChi,       d3OverlapddXdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldR_nldGradChi,        d3OverlapdR_nldR_nldGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldFdGradChi,           d3OverlapdR_nldFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChidGradChi,         d3OverlapdR_nldChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldChi_NL_BdGradChi,    d3OverlapdR_nldChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdR_nldGradChidGradChi,     d3OverlapdR_nldGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdFdGradChi,              d3OverlapdFdFdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChidGradChi,            d3OverlapdFdChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdChi_NL_BdGradChi,       d3OverlapdFdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdFdGradChidGradChi,        d3OverlapdFdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChidGradChi,           d3OverlapdChidChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidChi_NL_BdGradChi,      d3OverlapdChidChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChidGradChidGradChi,       d3OverlapdChidGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChi_NL_BdChi_NL_BdGradChi, d3OverlapdChi_NL_BdChi_NL_BdGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdChi_NL_BdGradChidGradChi,  d3OverlapdChi_NL_BdGradChidGradChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d3OverlapdGradChidGradChidGradChi,   d3OverlapdGradChidGradChidGradChi_answer ) );

}

BOOST_AUTO_TEST_CASE( test_computeLinearTraction ){

    floatVector dn = { 1, 2, 3 };

    floatVector dt = { -0.2, 4.5, .7 };

    floatVector parameters = { 3.6, 1.8 };

    floatVector answer = parameters[ 0 ] * dn + parameters[ 1 ] * dt;

    floatVector result;

    BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt, parameters, result ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    floatMatrix dtractionddn, dtractionddn_answer( answer.size( ), floatVector( dn.size( ), 0 ) );

    floatMatrix dtractionddt, dtractionddt_answer( answer.size( ), floatVector( dt.size( ), 0 ) );

    floatMatrix dtractiondp, dtractiondp_answer( answer.size( ), floatVector( parameters.size( ), 0 ) );

    floatMatrix d2tractionddndp, d2tractionddndp_answer( answer.size( ), floatVector( dn.size( ) * parameters.size( ), 0 ) );

    floatMatrix d2tractionddtdp, d2tractionddtdp_answer( answer.size( ), floatVector( dt.size( ) * parameters.size( ), 0 ) );

    floatVector result_2;

    BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt, parameters, result_2, dtractionddn, dtractionddt, dtractiondp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, result_2 ) );

    floatVector result_3;

    floatMatrix dtractionddn_3, dtractionddt_3, dtractiondp_3;

    BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt, parameters, result_3, dtractionddn_3, dtractionddt_3, dtractiondp_3, d2tractionddndp, d2tractionddtdp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, result_3 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dtractionddn, dtractionddn_3 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dtractionddt, dtractionddt_3 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dtractiondp, dtractiondp_3 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < dn.size( ); i++ ){

        floatVector delta( dn.size( ), 0 );

        delta[ i ] = eps * std::fabs( dn[ i ] ) + eps;

        floatVector rp, rm;

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn + delta, dt, parameters, rp ) );

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn - delta, dt, parameters, rm ) );

        for ( unsigned int j = 0; j < answer.size( ); j++ ){

            dtractionddn_answer[ j ][ i ] = ( rp[ j ] - rm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dtractionddnp, dtractionddnm;

        floatMatrix dtractionddtp, dtractionddtm;

        floatMatrix dtractiondpp, dtractiondpm;

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt + delta, parameters, rp, dtractionddnp, dtractionddtp, dtractiondpp ) );

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt - delta, parameters, rm, dtractionddnm, dtractionddtm, dtractiondpm ) );

        for ( unsigned int j = 0; j < answer.size( ); j++ ){

            for ( unsigned int k = 0; k < dn.size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( ( dtractionddnp[ j ][ k ] - dtractionddnm[ j ][ k ] ) / ( 2 * delta[ i ] ), 0. ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dtractionddn, dtractionddn_answer ) );

    for ( unsigned int i = 0; i < dt.size( ); i++ ){

        floatVector delta( dt.size( ), 0 );

        delta[ i ] = eps * std::fabs( dt[ i ] ) + eps;

        floatVector rp, rm;

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt + delta, parameters, rp ) );

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt - delta, parameters, rm ) );

        for ( unsigned int j = 0; j < answer.size( ); j++ ){

            dtractionddt_answer[ j ][ i ] = ( rp[ j ] - rm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dtractionddnp, dtractionddnm;

        floatMatrix dtractionddtp, dtractionddtm;

        floatMatrix dtractiondpp, dtractiondpm;

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt + delta, parameters, rp, dtractionddnp, dtractionddtp, dtractiondpp ) );

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt - delta, parameters, rm, dtractionddnm, dtractionddtm, dtractiondpm ) );

        for ( unsigned int j = 0; j < answer.size( ); j++ ){

            for ( unsigned int k = 0; k < dn.size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( ( dtractionddnp[ j ][ k ] - dtractionddnm[ j ][ k ] ) / ( 2 * delta[ i ] ), 0. ) );

            }

            for ( unsigned int k = 0; k < dt.size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( ( dtractionddtp[ j ][ k ] - dtractionddtm[ j ][ k ] ) / ( 2 * delta[ i ] ), 0. ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dtractionddt, dtractionddt_answer ) );

    for ( unsigned int i = 0; i < parameters.size( ); i++ ){

        floatVector delta( parameters.size( ), 0 );

        delta[ i ] = eps * std::fabs( parameters[ i ] ) + eps;

        floatVector rp, rm;

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt, parameters + delta, rp ) );

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt, parameters - delta, rm ) );

        for ( unsigned int j = 0; j < answer.size( ); j++ ){

            dtractiondp_answer[ j ][ i ] = ( rp[ j ] - rm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatMatrix dtractionddnp, dtractionddnm;

        floatMatrix dtractionddtp, dtractionddtm;

        floatMatrix dtractiondpp, dtractiondpm;

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt, parameters + delta, rp, dtractionddnp, dtractionddtp, dtractiondpp ) );

        BOOST_CHECK_NO_THROW( tractionSeparation::computeLinearTraction( dn, dt, parameters - delta, rm, dtractionddnm, dtractionddtm, dtractiondpm ) );

        for ( unsigned int j = 0; j < answer.size( ); j++ ){

            for ( unsigned int k = 0; k < dn.size( ); k++ ){

                d2tractionddndp_answer[ j ][ parameters.size( ) * k + i ] = ( dtractionddnp[ j ][ k ] - dtractionddnm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < dt.size( ); k++ ){

                d2tractionddtdp_answer[ j ][ parameters.size( ) * k + i ] = ( dtractionddtp[ j ][ k ] - dtractionddtm[ j ][ k ] ) / ( 2 * delta[ i ] );

            }

            for ( unsigned int k = 0; k < parameters.size( ); k++ ){

                BOOST_CHECK( vectorTools::fuzzyEquals( ( dtractiondpp[ j ][ k ] - dtractiondpm[ j ][ k ] ) / ( 2 * delta[ i ] ), 0. ) );

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dtractiondp, dtractiondp_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2tractionddndp, d2tractionddndp_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2tractionddtdp, d2tractionddtdp_answer ) );

}
