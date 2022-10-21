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

    floatType L_answer = 0.595416169431248;

    floatType L;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, L ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L, L_answer ) );

    floatVector dLdX, dLdchi_nl, dLdxi_t;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, L, dLdX, dLdchi_nl, dLdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L, L_answer ) );

    floatType L_2;

    floatVector dLdX_2, dLdchi_nl_2, dLdxi_t_2;

    floatVector d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t, d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t, d2Ldxi_tdxi_t;

    BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t, L_2, dLdX_2, dLdchi_nl_2, dLdxi_t_2,
                                                                        d2LdXdX, d2LdXdchi_nl, d2LdXdxi_t,
                                                                        d2Ldchi_nldchi_nl, d2Ldchi_nldxi_t,
                                                                        d2Ldxi_tdxi_t ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( L_2, L_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdX_2, dLdX ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdchi_nl_2, dLdchi_nl ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdxi_t_2, dLdxi_t ) );

    floatType eps = 1e-6;

    floatVector dLdX_answer( X.size( ), 0 );

    floatVector dLdchi_nl_answer( chi_nl.size( ), 0 );

    floatVector dLdxi_t_answer( xi_t.size( ), 0 );

    floatVector d2LdXdX_answer( X.size( ) * X.size( ), 0 );

    floatVector d2LdXdchi_nl_answer( X.size( ) * chi_nl.size( ), 0 );

    floatVector d2LdXdxi_t_answer( X.size( ) * xi_t.size( ), 0 );

    floatVector d2Ldchi_nldchi_nl_answer( chi_nl.size( ) * chi_nl.size( ), 0 );

    floatVector d2Ldchi_nldxi_t_answer( chi_nl.size( ) * xi_t.size( ), 0 );

    floatVector d2Ldxi_tdxi_t_answer( xi_t.size( ) * xi_t.size( ), 0 );

    for ( unsigned int i = 0; i < X.size( ); i++ ){

        floatVector delta( X.size( ), 0 );

        delta[ i ] += eps * std::abs( X[ i ] ) + eps;

        floatType Lp, Lm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X + delta, chi_nl, xi_t, Lp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X - delta, chi_nl, xi_t, Lm ) );

        dLdX_answer[ i ] += ( Lp - Lm ) / ( 2 * delta[ i ] );

        floatVector dLdXp, dLdXm;

        floatVector dLdchi_nlp, dLdchi_nlm;

        floatVector dLdxi_tp, dLdxi_tm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X + delta, chi_nl, xi_t, Lp, dLdXp, dLdchi_nlp, dLdxi_tp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X - delta, chi_nl, xi_t, Lm, dLdXm, dLdchi_nlm, dLdxi_tm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            d2LdXdX_answer[ X.size( ) * j + i ] = ( dLdXp[ j ] - dLdXm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdX, dLdX_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdX, d2LdXdX_answer ) );

    for ( unsigned int i = 0; i < chi_nl.size( ); i++ ){

        floatVector delta( chi_nl.size( ), 0 );

        delta[ i ] += eps * std::abs( chi_nl[ i ] ) + eps;

        floatType Lp, Lm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl + delta, xi_t, Lp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl - delta, xi_t, Lm ) );

        dLdchi_nl_answer[ i ] += ( Lp - Lm ) / ( 2 * delta[ i ] );

        floatVector dLdXp, dLdXm;

        floatVector dLdchi_nlp, dLdchi_nlm;

        floatVector dLdxi_tp, dLdxi_tm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl + delta, xi_t, Lp, dLdXp, dLdchi_nlp, dLdxi_tp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl - delta, xi_t, Lm, dLdXm, dLdchi_nlm, dLdxi_tm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            d2LdXdchi_nl_answer[ chi_nl.size( ) * j + i ] = ( dLdXp[ j ] - dLdXm[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < chi_nl.size( ); j++ ){

            d2Ldchi_nldchi_nl_answer[ chi_nl.size( ) * j + i ] = ( dLdchi_nlp[ j ] - dLdchi_nlm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdchi_nl, dLdchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdchi_nl, d2LdXdchi_nl_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldchi_nl, d2Ldchi_nldchi_nl_answer ) );

    floatMatrix _d2Ldchi_nldxi_t_answer( chi_nl.size( ), floatVector( xi_t.size( ), 0 ) );

    for ( unsigned int i = 0; i < xi_t.size( ); i++ ){

        floatVector delta( xi_t.size( ), 0 );

        delta[ i ] += eps * std::abs( xi_t[ i ] ) + eps;

        floatType Lp, Lm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t + delta, Lp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t - delta, Lm ) );

        dLdxi_t_answer[ i ] += ( Lp - Lm ) / ( 2 * delta[ i ] );

        floatVector dLdXp, dLdXm;

        floatVector dLdchi_nlp, dLdchi_nlm;

        floatVector dLdxi_tp, dLdxi_tm;

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t + delta, Lp, dLdXp, dLdchi_nlp, dLdxi_tp ) );

        BOOST_CHECK( !tractionSeparation::computeOverlapDistanceLagrangian( X, chi_nl, xi_t - delta, Lm, dLdXm, dLdchi_nlm, dLdxi_tm ) );

        for ( unsigned int j = 0; j < X.size( ); j++ ){

            d2LdXdxi_t_answer[ xi_t.size( ) * j + i ] = ( dLdXp[ j ] - dLdXm[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < chi_nl.size( ); j++ ){

            d2Ldchi_nldxi_t_answer[ xi_t.size( ) * j + i ] = ( dLdchi_nlp[ j ] - dLdchi_nlm[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < xi_t.size( ); j++ ){

            d2Ldxi_tdxi_t_answer[ xi_t.size( ) * j + i ] = ( dLdxi_tp[ j ] - dLdxi_tm[ j ] ) / ( 2 * delta[ i ] );
        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dLdxi_t, dLdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2LdXdxi_t, d2LdXdxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldchi_nldxi_t, d2Ldchi_nldxi_t_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2Ldxi_tdxi_t, d2Ldxi_tdxi_t_answer ) );

}
