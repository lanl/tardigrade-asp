/**
  ******************************************************************************
  * \file cpp_stub.h
  ******************************************************************************
  * A C++ library for printing messages to stdout. Used as a stub repo example.
  ******************************************************************************
  */

#ifndef CPP_STUB_H
    #define CPP_STUB_H
#endif

#include<string>
#include<iostream>
#include<vector>

namespace cppStub{

    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    /// Say hello
    /// @param message The message to print
    void sayHello(std::string message);

    void abaqusInterface( floatVector &STRESS,       floatVector &STATEV,       floatVector &DDSDDE,       floatType &SSE,           floatType &SPD,
                          floatType &SCD,            floatType &RPL,            floatVector &DDSDDT,       floatVector &DRPLDE,      floatType &DRPLDT,
                          const floatVector &STRAN,  const floatVector &DSTRAN, const floatVector &TIME,   const floatType &DTIME,   const floatType &TEMP,
                          const floatType &DTEMP,    const floatVector &PREDEF, const floatVector &DPRED,  const char *CMNAME,       const int &NDI,
                          const int &NSHR,           const int &NTENS,          const int &NSTATV,         const floatVector &PROPS, const int &NPROPS,
                          const floatVector &COORDS, const floatVector &DROT,   floatType &PNEWDT,         const floatType &CELENT,  const floatVector &DFGRD0,
                          const floatVector &DFGRD1, const int &NOEL,           const int &NPT,            const int &LAYER,         const int &KSPT,
                          const int *JSTEP,          const int &KINC );

}
