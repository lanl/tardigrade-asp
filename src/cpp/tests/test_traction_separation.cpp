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

                    d4LdXdxi_tdR_nldR_nl_answer[ xi_t.size( ) * xi_t.size( ) * j + xi_t.size( ) * k + l + i ]
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

    floatMatrix RHSTERM;

    BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl, distance_2, dDistancedchi_nl_2, dDistancedxi_t_2, dDistancedR_nl_2,
                                                            d2Distancedchi_nldchi_nl, d2Distancedchi_nldxi_t, d2Distancedchi_nldR_nl,
                                                            d2Distancedxi_tdxi_t, d2Distancedxi_tdR_nl,
                                                            d2DistancedR_nldR_nl, RHSTERM ) );

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

    floatMatrix RHSTERM_GRAD;

    BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl, distance_3, dDistancedchi_nl_3, dDistancedxi_t_3, dDistancedR_nl_3,
                                                            d2Distancedchi_nldchi_nl_3, d2Distancedchi_nldxi_t_3, d2Distancedchi_nldR_nl_3,
                                                            d2Distancedxi_tdxi_t_3, d2Distancedxi_tdR_nl_3,
                                                            d2DistancedR_nldR_nl_3,
                                                            d3Distancedchi_nldchi_nldchi_nl, d3Distancedchi_nldchi_nldxi_t, d3Distancedchi_nldchi_nldR_nl,
                                                            d3Distancedchi_nldxi_tdxi_t, d3Distancedchi_nldxi_tdR_nl,
                                                            d3Distancedchi_nldR_nldR_nl,
                                                            d3Distancedxi_tdxi_tdxi_t, d3Distancedxi_tdxi_tdR_nl,
                                                            d3Distancedxi_tdR_nldR_nl,
                                                            d3DistancedR_nldR_nldR_nl, RHSTERM_GRAD ) );

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

    floatMatrix RHSTERM_GRAD_ANSWER( 4, floatVector( chi_nl.size( ) * chi_nl.size( ), 0 ) );

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

        floatMatrix RHSTERMp, RHSTERMm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl + delta, xi_t, R_nl, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp,
                                                                d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp,
                                                                d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp,
                                                                d2DistancedR_nldR_nlp, RHSTERMp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl - delta, xi_t, R_nl, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm,
                                                                d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm,
                                                                d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm,
                                                                d2DistancedR_nldR_nlm, RHSTERMm ) );

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

        floatMatrix RHSTERMp, RHSTERMm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t + delta, R_nl, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp,
                                                                d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp,
                                                                d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp,
                                                                d2DistancedR_nldR_nlp, RHSTERMp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t - delta, R_nl, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm,
                                                                d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm,
                                                                d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm,
                                                                d2DistancedR_nldR_nlm, RHSTERMm ) );

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

        floatMatrix RHSTERMp, RHSTERMm;

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl + delta, dp, dDistancedchi_nlp, dDistancedxi_tp, dDistancedR_nlp,
                                                                d2Distancedchi_nldchi_nlp, d2Distancedchi_nldxi_tp, d2Distancedchi_nldR_nlp,
                                                                d2Distancedxi_tdxi_tp, d2Distancedxi_tdR_nlp,
                                                                d2DistancedR_nldR_nlp, RHSTERMp ) );

        BOOST_CHECK( !tractionSeparation::solveOverlapDistance( chi_nl, xi_t, R_nl - delta, dm, dDistancedchi_nlm, dDistancedxi_tm, dDistancedR_nlm,
                                                                d2Distancedchi_nldchi_nlm, d2Distancedchi_nldxi_tm, d2Distancedchi_nldR_nlm,
                                                                d2Distancedxi_tdxi_tm, d2Distancedxi_tdR_nlm,
                                                                d2DistancedR_nldR_nlm, RHSTERMm ) );

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

        for ( unsigned int j = 0; j < 4; j++ ){

            for ( unsigned int k = 0; k < chi_nl.size( ); k++ ){

                for ( unsigned int l = 0; l < chi_nl.size( ); l++ ){

                    RHSTERM_GRAD_ANSWER[ j ][ chi_nl.size( ) * k + l + i ] = ( RHSTERMp[ j ][ chi_nl.size( ) * k + l ] - RHSTERMm[ j ][ chi_nl.size( ) * k + l ] ) / ( 2 * delta );

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
    std::cout << "d3DistancedR_nldR_nldR_nl:\n"; vectorTools::print( d3DistancedR_nldR_nldR_nl );
    std::cout << "d3DistancedR_nldR_nldR_nl_answer:\n"; vectorTools::print( d3DistancedR_nldR_nldR_nl_answer );

    BOOST_CHECK( vectorTools::fuzzyEquals( RHSTERM_GRAD, RHSTERM_GRAD_ANSWER ) );
    std::cout << "RHSTERM_GRAD:\n"; vectorTools::print( RHSTERM_GRAD );
    std::cout << "RHSTERM_GRAD_ANSWER:\n"; vectorTools::print( RHSTERM_GRAD_ANSWER );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdXi_1, dOverlapdXi_1_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dXi_1, d2OverlapdXi_1dXi_1_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapddX, dOverlapddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1ddX, d2OverlapdXi_1ddX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXddX, d2OverlapddXddX_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdR_nl, dOverlapdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dR_nl, d2OverlapdXi_1dR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdR_nl, d2OverlapddXdR_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldR_nl, d2OverlapdR_nldR_nl_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdF, dOverlapdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dF, d2OverlapdXi_1dF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdF, d2OverlapddXdF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldF, d2OverlapdR_nldF_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdF, d2OverlapdFdF_answer ) );

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

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dOverlapdChi, dOverlapdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdXi_1dChi, d2OverlapdXi_1dChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapddXdChi, d2OverlapddXdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdR_nldChi, d2OverlapdR_nldChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdFdChi, d2OverlapdFdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2OverlapdChidChi, d2OverlapdChidChi_answer ) );

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

}
