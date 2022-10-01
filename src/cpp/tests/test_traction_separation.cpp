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

    floatType eps = 1e-6;

    floatVector denergydn_answer( dn.size( ), 0 );

    floatVector denergydt_answer( dt.size( ), 0 );

    floatVector d2energydndn_answer( dn.size( ) * dn.size( ), 0 );

    floatVector d2energydndt_answer( dn.size( ) * dt.size( ), 0 );

    floatVector d2energydtdt_answer( dt.size( ) * dt.size( ), 0 );

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

            d2energydndn_answer[ dn.size( ) * i + j ] += ( dednp[ j ] - dednm[ j ] ) / ( 2 * delta[ i ] );

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

            d2energydndt_answer[ dt.size( ) * i + j ] += ( dednp[ j ] - dednm[ j ] ) / ( 2 * delta[ i ] );

            d2energydtdt_answer[ dt.size( ) * i + j ] += ( dedtp[ j ] - dedtm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( denergydt, denergydt_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydndt, d2energydndt_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2energydtdt, d2energydtdt_answer ) );

}
