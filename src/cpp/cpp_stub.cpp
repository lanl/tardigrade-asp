/**
  ******************************************************************************
  * \file cpp_stub.cpp
  ******************************************************************************
  * A C++ library for printing messages to stdout. Used as a stub repo example.
  ******************************************************************************
  */

#include<cpp_stub.h>

namespace cppStub{

    /// Say hello
    /// @param message The message to print
    void sayHello(std::string message) {
        std::cout << "Hello " << message << std::endl;
    }

    void abaqusInterface( floatVector &stress,             floatVector &statev,        floatMatrix &ddsdde,       floatType &SSE,            floatType &SPD,
                          floatType &SCD,                  floatType &RPL,             floatVector &ddsddt,       floatVector &drplde,       floatType &DRPLDT,
                          const floatVector &strain,       const floatVector &dstrain, const floatVector &time,   const floatType &DTIME,    const floatType &TEMP,
                          const floatType &DTEMP,          const floatVector &predef,  const floatVector &dpred,  const std::string &cmname, const int &NDI,
                          const int &NSHR,                 const int &NTENS,           const int &NSTATV,         const floatVector &props,  const int &NPROPS,
                          const floatVector &coords,       const floatMatrix &drot,    floatType &PNEWDT,         const floatType &CELENT,   const floatMatrix &dfgrd0,
                          const floatMatrix &dfgrd1,       const int &NOEL,            const int &NPT,            const int &LAYER,          const int &KSPT,
                          const std::vector< int > &jstep, const int &KINC ){
        /*!
         * A template Abaqus UMAT c++ interface using c++ STL types. Variables in all caps reference ABAQUS FORTRAN
         * memory directly. Variables in lower case are native c++ type conversions stored separately from the original
         * ABAQUS FORTRAN memory.
         */

        //Perform type conversions and variable calculations for constitutive model interface
        //Common example conversions is to expand NTENS vectors (stress, strain, etc) and arrays (ddsdde) to the Abaqus
        //Voigt notation in 1x6 or 6x6 shape.

        //Call constitutive model
        sayHello( "Abaqus" );

        //Re-pack original c++ UMAT variables (lower case) from constitutive model interface

        return;
    }

}
