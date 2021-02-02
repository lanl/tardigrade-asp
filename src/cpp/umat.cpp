/**
  ******************************************************************************
  * \file umat.cpp
  ******************************************************************************
  * The Abaqus UMAT c++ interface template
  ******************************************************************************
  */

#include<umat.h>

extern "C" void umat_(double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                      double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                      const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                      const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                      const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                      const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                      const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                      const int *JSTEP,     const int &KINC){
    /*!
     * The Abaqus UMAT c++ interface for the (V)iscous (I)sotropic (Por)osity based model.
     *
     * The variables defined in this interface are described more completely in the Abaqus User Subroutines
     * manual entry for UMAT user subroutines.
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

    if (KINC ==1 && NOEL == 1 && NPT == 1){
        std::cout << std::endl << "Sign of life" << std::endl;
        std::cout << std::endl;
    }

    //Map FORTRAN UMAT variables to C++ types as necessary. Use case sensitivity to distinguish.
    //TODO: Decide if case sensitive variable names is a terrible idea or not
    //Vectors are straight forward
    std::vector<double> stress( STRESS, STRESS + NTENS );
    std::vector<double> statev( STATEV, STATEV + NSTATV );
    std::vector<double> ddsddt( DDSDDT, DDSDDT + NTENS );
    const std::vector<double> strain( STRAN, STRAN + NTENS );
    const std::vector<double> dstrain( DSTRAN, DSTRAN + NTENS );
    const std::vector<double> time( TIME, TIME + 2 );
    const std::vector<double> predef( PREDEF, PREDEF + 1 );
    const std::vector<double> dpred( DPRED, DPRED + 1 );
    //TODO: cmname?
    const std::vector<double> props( PROPS, PROPS + NPROPS );
    const std::vector<double> coords( COORDS, COORDS + 3 );
    const std::vector<int> jstep( JSTEP, JSTEP + 4 );
    //Matrices require careful row/column major conversions
    std::vector< std::vector< double > > ddsdde;
          std::vector< std::vector< double > > drot;
          std::vector< std::vector< double > > dfgrd0;
          std::vector< std::vector< double > > dfgrd1;

    std::vector< std::vector< double > > dummy = columnToRowMajor(DDSDDE, NTENS, NTENS);

    // Perform a column to row major conversion to populate c++ ddsdde
    int height = NTENS;
    int width = NTENS;
    int column_major_index;
    for (int row = 0; row < height; row++){
        std::vector< double > ddsdde_row;
        for (int col = 0; col < width; col++){
            column_major_index = row*width + col;
            ddsdde_row.push_back(*(DDSDDE + column_major_index));
        }
        ddsdde.push_back(ddsdde_row);
    }

    // Perform a column to row major conversion to populate c++ two dimensional arrays of 3x3
    height = 3;
    width = 3;
    column_major_index;
    for (int row = 0; row < height; row++){
        std::vector< double > drot_row;
        std::vector< double > dfgrd0_row;
        std::vector< double > dfgrd1_row;
        for (int col = 0; col < width; col++){
            column_major_index = row*width + col;
            drot_row.push_back(*(DROT + column_major_index));
            dfgrd0_row.push_back(*(DFGRD0 + column_major_index));
            dfgrd1_row.push_back(*(DFGRD1 + column_major_index));
        }
        drot.push_back(drot_row);
        dfgrd0.push_back(dfgrd0_row);
        dfgrd1.push_back(dfgrd1_row);
    }

    // Print two dimensional array(s) to check work
    if (NOEL == 1 && NPT == 1){
        std::cout << "KINC: " << KINC << "; DFGRD0" << std::endl;
        height = 3;
        width = 3;
        column_major_index;
        for (int row = 0; row < height; row++){
            std::vector< double > vector_row;
            for (int col = 0; col < width; col++){
                column_major_index = row*width + col;
                std::cout << DFGRD0[column_major_index] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "KINC: " << KINC << "; dfgrd0" << std::endl;
        for (int row = 0; row < height; row++){
            std::vector< double > vector_row;
            for (int col = 0; col < width; col++){
                std::cout << dfgrd0[row][col] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    return;
}

std::vector< std::vector< double > > columnToRowMajor(double *myPointer, const int &width, const int &height){
    std::vector< std::vector< double > > row_major;
    return row_major;
}
