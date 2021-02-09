/**
  ******************************************************************************
  * \file umat.cpp
  ******************************************************************************
  * The Abaqus UMAT c++ interface template
  ******************************************************************************
  */

#include<umat.h>

extern "C" void umat_( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                       double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                       const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                       const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                       const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                       const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                       const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                       const int *JSTEP,     const int &KINC ){
    /*!
     * A template Abaqus UMAT c++ interface.
     *
     * The variables defined in this interface are described more completely in the Abaqus User Subroutines
     * manual entry for UMAT user subroutines.
     *
     * The 80 character assumption for string conversion is based on the data line limits on character strings as
     * described in the Abaqus User's Manual > Introduction & Spatial Modeling > Input Syntax Rules and the sample
     * Fortran UMAT in the Abaqus User's Manual > User Subroutines > Abaqus/Standard User Subroutines > UMAT.
     *
     * \param *STRESS: Cauchy stress tensor at beginning of time increment stored in vector form, \f$ \sigma \f$.
     * \param *STATEV: State variable vector.
     * \param *DDSDDE: Jacobian matrix \f$ \frac{\delta \Delta \sigma}{\delta \Delta \epsilon} \f$.
     * \param &SSE: Specific elastic strain energy.
     * \param &SPD: Specific plastic dissipation.
     * \param &SCD: Specific "creep" dissipation.
     * \param &RPL: Volumetric heat generation per unit time at the end of the time increment.
     * \param *DDSDDT: Variation of stress increment with respect to temperature.
     * \param *DRPLDE: Variation of RPL with respect to strain increment.
     * \param &DRPLDT: Variation of RPL with respect to temperature .
     * \param *STRAN: Strain tensor at the beginning of the time increment stored in vector form, \f$ \epsilon \f$.
     * \param *DSTRAN: Strain increment tensor stored in vector form.
     * \param *TIME: Time vector at beginning of time increment = {Step time, Total time}.
     * \param &DTIME: Time increment.
     * \param &TEMP: Temperature at the start of the time increment.
     * \param &DTEMP: Temperature increment.
     * \param *PREDEF: Interpolated predefined field variables at current integration point.
     * \param *DPRED: Change in predefined field variables.
     * \param *CMNAME: User defined material name. Left justified as passed by FORTRAN.
     * \param &NDI: Number of direct stress components at this integration point.
     * \param &NSHR: Number of engineering shear stress components at this integration point.
     * \param &NTENS: Size of the stress and strain component array. NTENS = NDI + NSHR
     * \param &NSTATEV: Number of state variables for this material, CMNAME.
     * \param *PROPS: Material model constants defined as part of the *MATERIAL keyword in the input file.
     * \param &NPROPS: Number of user defined material constants.
     * \param *COORDS: Coordinates of the current Gauss point.
     * \param *DROT: Rigid body rotation increment matrix.
     * \param &PNEWDT: Ratio of suggested new time increment to current time increment, e.g. the next increment's DTIME.
     * \param &CELENT: Characteristic element length.
     * \param *DFGRD0: Deformation gradient matrix at the beginning of the time increment.
     * \param *DFGRD1: Deformation gradient at the end of the current time increment.
     * \param &NOEL: Current element number.
     * \param &NPT: Integration point number.
     * \param &LAYER: Layer number for composite shells and layered solids.
     * \param &KSPT: Section point number within the current layer.
     * \param *JSTEP: Meta data vector = {Step number, procedure type key, NLGEOM switch, linear perturbation switch}
     * \param &KINC: Increment number.
     */

    //Define the size of geometry related tensor dimensions
    const int geometricSize = 3;

    //Map FORTRAN UMAT variables to C++ types as necessary. Use case sensitivity to distinguish.
    //TODO: Decide if case sensitive variable names is a terrible idea or not
    //Vectors can be created directly with pointer arithmetic
    std::vector< double > stress( STRESS, STRESS + NTENS );
    std::vector< double > statev( STATEV, STATEV + NSTATV );
    std::vector< double > ddsddt( DDSDDT, DDSDDT + NTENS );
    std::vector< double > drplde( DRPLDE, DRPLDE + NTENS );
    const std::vector< double > strain( STRAN, STRAN + NTENS );
    const std::vector< double > dstrain( DSTRAN, DSTRAN + NTENS );
    const std::vector< double > time( TIME, TIME + 2 );
    const std::vector< double > predef( PREDEF, PREDEF + 1 );
    const std::vector< double > dpred( DPRED, DPRED + 1 );
    const std::string cmname( FtoCString( 80, CMNAME ) );
    const std::vector< double > props( PROPS, PROPS + NPROPS );
    const std::vector< double > coords( COORDS, COORDS + geometricSize );
    const std::vector< int > jstep( JSTEP, JSTEP + 4 );
    //Fortran two-dimensional arrays require careful column to row major conversions to c++ types
    std::vector< std::vector< double > > ddsdde = columnToRowMajor( DDSDDE, NTENS, NTENS );
    const std::vector< std::vector< double > > drot = columnToRowMajor( DROT, geometricSize, geometricSize );
    const std::vector< std::vector< double > > dfgrd0 = columnToRowMajor( DFGRD0, geometricSize, geometricSize );
    const std::vector< std::vector< double > > dfgrd1 = columnToRowMajor( DFGRD1, geometricSize, geometricSize );

    //Call the appropriate subroutine interface
    //Show example use of c++ library in UMAT
    if ( KINC == 1 && NOEL == 1 && NPT == 1 ){
        cppStub::sayHello( "Abaqus" );
    }

    //Re-pack C++ objects into FORTRAN memory to return values to Abaqus
    //Scalars were passed by reference and will update correctly
    //Vectors don't require row/column major considerations
    for ( int row = 0; row < NTENS; row++ ){
        STRESS[row] = stress[row];
        DDSDDT[row] = ddsddt[row];
        DRPLDE[row] = drplde[row];
    }
    for ( int row = 0; row < NSTATV; row++ ){
        STATEV[row] = statev[row];
    }
    //Arrays require vector of vector to column major conversion
    int column_major_index;
    for ( int row = 0; row < geometricSize; row++ ){
        for ( int col = 0; col < geometricSize; col++ ){
            column_major_index = col*geometricSize + row;
            DDSDDE[column_major_index] = ddsdde[row][col];
        }
    }

    return;
}

template< typename T >
std::vector< std::vector< double > > columnToRowMajor( T *column_major, const int &width, const int &height ){
    /*!
     * Convert column major two dimensional arrays to row major.
     *
     * Specifically, convert pointers to Fortran column major arrays to c++ row major vector of vectors.
     *
     * \param *column_major: The pointer to the start of a column major array
     * \param &width: The width of the array, e.g. number of columns
     * \param &height: The height of the array, e.g. number of rows
     */
    std::vector< std::vector< double > > row_major;
    int column_major_index;
    for ( int row = 0; row < height; row++ ){
        std::vector< double > row_vector;
        for ( int col = 0; col < width; col++ ){
            column_major_index = col*height + row;
            row_vector.push_back( *( column_major + column_major_index ) );
        }
        row_major.push_back( row_vector );
    }
    return row_major;
}

char *FtoCString( int stringLength, const char* fString ){
    /*!
     * Converts a Fortran string to C-string. Trims trailing white space during processing.
     *
     * Code excerpt from a c++ Abaqus FILM subroutine in the Abaqus Knowledge Base:
     * https://kb.dsxclient.3ds.com/mashup-ui/page/resultqa?from=search%3fq%3dwriting%2bsubroutine%2bc%252B%252B&id=QA00000008005e&q=writing%20subroutine%20c%2B%2B 
     *
     * TODO: update coding style to match project.
     *
     * \param stringLength: The length of the Fortran string.
     * \param *fString: The pointer to the start of the Fortran string.
     */
    int stringLen = stringLength;
    for ( int k1 = stringLength - 1; k1 >= 0; k1-- )
	{
	    if ( fString[ k1 ] != ' ' ) break;
	    stringLen = k1;
	}
    char* cString  = new char [ stringLen + 1 ];
    memcpy ( cString, fString, stringLen );
    cString[ stringLen ] = '\0';
    return cString;
}
