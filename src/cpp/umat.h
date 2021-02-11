/**
  ******************************************************************************
  * \file umat.h
  ******************************************************************************
  * The declarations and definitions required for an Abaqus UMAT c++ interface.
  ******************************************************************************
  */

#include<iostream>
#include<vector>
#include<string.h>
#include<stdio.h>

#include<vector_tools.h>

#include<cpp_stub.h>

extern "C" void UMAT( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                      double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                      const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                      const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                      const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                      const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                      const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                      const int *JSTEP,     const int &KINC );

template< typename T >
std::vector< std::vector< T > > columnToRowMajor( const T *column_major,  const int &height, const int &width ){
    /*!
     * Convert column major two dimensional arrays to row major.
     *
     * Specifically, convert pointers to Fortran column major arrays to c++ row major vector of vectors.
     *
     * \param *column_major: The pointer to the start of a column major array
     * \param &height: The height of the array, e.g. number of rows
     * \param &width: The width of the array, e.g. number of columns
     * \return row_major: A c++ row major vector of vectors
     */
    std::vector< std::vector< T > > row_major;
    int column_major_index;
    for ( int row = 0; row < height; row++ ){
        std::vector< T > row_vector;
        for ( int col = 0; col < width; col++ ){
            column_major_index = col*height + row;
            row_vector.push_back( *( column_major + column_major_index ) );
        }
        row_major.push_back( row_vector );
    }
    return row_major;
}

template< typename T >
void rowToColumnMajor( T *column_major, const std::vector< std::vector< T > > &row_major_array,
                       const int &height, const int &width ){
    /*!
     * Convert row major two dimensional arrays to column major
     *
     * Specifically, c++ row major vector of vectors or arrays to Fortran column major arrays using the column major
     * pointer.
     *
     * \param *column_major: The pointer to the start of a column major array
     * \param &row_major_array: A c++ two dimensional, row major vector of vectors
     * \param &height: The height of the array, e.g. number of rows
     * \param &width: The width of the array, e.g. number of columns
     */
    int column_major_index;
    for ( int row = 0; row < height; row++ ){
        for ( int col = 0; col < width; col++ ){
            column_major_index = col*height + row;
            column_major[column_major_index] = row_major_array[row][col];
        }
    }

    return;
}

template< typename T >
void rowToColumnMajor( T *column_major, const std::vector< T > &row_major, const int &height, const int &width ){
    /*!
     * Convert row major two dimensional arrays stored as vector to column major array
     *
     * Specifically, c++ row major vector to Fortran column major arrays using the column major pointer.
     *
     * \param *column_major: The pointer to the start of a column major array
     * \param &row_major_array: A c++ two dimensional array stored as row major vector
     * \param &height: The height of the array, e.g. number of rows. The c++ row count (1) for 1D arrays.
     * \param &width: The width of the array, e.g. number of columns. The c++ column count (size) for 1D arrays.
     */
    int row_major_index;
    int column_major_index;
    for ( int row = 0; row < height; row++ ){
        for ( int col = 0; col < width; col++ ){
            row_major_index = row*width + col;
            column_major_index = col*height + row;
            column_major[column_major_index] = row_major[row_major_index];
        }
    }
}

template< typename T >
std::vector< T > expandAbaqusStandardStressVector( const std::vector< T > &abaqus_vector,
                                                   const int &NDI, const int &NSHR ){
    /*!
     * Expand stress and strain type components to full Abaqus vectors
     *
     * See the Abaqus documentation > Introduction & Spatial Modeling > Conventions chapter > Convention used for stress
     * and strain components.
     *
     * The stress vector components for Abaqus/Standard (UMAT) are
     *
     * \f$ \left { \sigma_{11}, \sigma_{22}, \sigma_{33}, \tau_{12}, \tau_{13}, \tau_{23} } \f$
     *
     * and the strain vector components match as
     *
     * \f$ \left { \epsilon_{11}, \epsilon_{22}, \epsilon_{33}, \gamma_{12}, \gamma_{13}, \gamma_{23} } \f$
     *
     * where components that are zero-valued by definition, e.g. plane stress, are omitted. The shear strain is the
     * engineering shear strain where
     *
     * \f$ \gamma_{ij} = \epsilon_{ij} + \epsilon_{ji} \f$
     *
     * \param &abaqus_vector: an abaqus stress-type vector with no by-definition-zero components. Length NDI + NSHR.
     * \param &NDI: The number of direct components.
     * \param &NSHR: The number of shear components.
     * \returns vector_expansion: c++ type vector of length 6.
     */

    //Initialize expanded vector with zero values
    std::vector< T > vector_expansion( 6, 0 );

    //Unpack direct components of Abaqus/Standard stress-type vector
    for ( int index = 0; index < NDI; index++ ){
        vector_expansion[ index ] = abaqus_vector[ index ];
    }

    //Unpack shear components of Abaqus/Standard stress-type vector
    for ( int index = 0; index < NSHR; index++ ){
        vector_expansion[ 3 + index ] = abaqus_vector[ NDI + index ];
    }

    return vector_expansion;
}

template< typename T >
std::vector< T > contractAbaqusStandardStressVector( const std::vector< T > &full_abaqus_vector,
                                                     const int &NDI, const int &NSHR ){
    /*!
     * Contract stress and strain type components from full Abaqus vectors
     *
     * See the Abaqus documentation > Introduction & Spatial Modeling > Conventions chapter > Convention used for stress
     * and strain components.
     *
     * The stress vector components for Abaqus/Standard (UMAT) are
     *
     * \f$ \left { \sigma_{11}, \sigma_{22}, \sigma_{33}, \tau_{12}, \tau_{13}, \tau_{23} } \f$
     *
     * and the strain vector components match as
     *
     * \f$ \left { \epsilon_{11}, \epsilon_{22}, \epsilon_{33}, \gamma_{12}, \gamma_{13}, \gamma_{23} } \f$
     *
     * where components that are zero-valued by definition, e.g. plane stress, are omitted. The shear strain is the
     * engineering shear strain where
     *
     * \f$ \gamma_{ij} = \epsilon_{ij} + \epsilon_{ji} \f$
     *
     * \param &full_abaqus_vector: a previously expanded abaqus stress-type vector. Length 6.
     * \param &NDI: The number of direct components.
     * \param &NSHR: The number of shear components.
     * \returns vector_contraction: c++ type vector of length NDI + NSHR.
     */

    //Initialize contracted vector
    std::vector< T > vector_contraction( NDI + NSHR );

    //Pack non-zero direct components of Abaqus/Standard stress-type vector
    for ( int index = 0; index < NDI; index++ ){
        vector_contraction[ index ] = full_abaqus_vector[ index ];
    }

    //Pack non-zero shear components of Abaqus/Standard stress-type vector
    for ( int index = 0; index < NSHR; index++ ){
        vector_contraction[ NDI + index ] = full_abaqus_vector[ 3 + index ];
    }

    return vector_contraction;
}

template< typename T >
std::vector< std::vector < T > > contractAbaqusStandardNTENSMatrix( const std::vector< std::vector< T > > &full_abaqus_matrix, 
                                                                    const int &NDI, const int &NSHR ){
    std::vector< std::vector< T > > matrix_contraction( 6, std::vector< double >( 6, -666. ) );
    return matrix_contraction;
}

char *FtoCString( int stringLength, const char* fString );
