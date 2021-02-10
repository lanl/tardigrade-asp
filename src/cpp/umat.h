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

char *FtoCString( int stringLength, const char* fString );
