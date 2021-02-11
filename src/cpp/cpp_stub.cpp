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

    void abaqusInterface( floatVector &STRESS,       floatVector &STATEV,       floatVector &DDSDDE,       floatType &SSE,           floatType &SPD,
                          floatType &SCD,            floatType &RPL,            floatVector &DDSDDT,       floatVector &DRPLDE,      floatType &DRPLDT,
                          const floatVector &STRAN,  const floatVector &DSTRAN, const floatVector &TIME,   const floatType &DTIME,   const floatType &TEMP,
                          const floatType &DTEMP,    const floatVector &PREDEF, const floatVector &DPRED,  const char *CMNAME,       const int &NDI,
                          const int &NSHR,           const int &NTENS,          const int &NSTATV,         const floatVector &PROPS, const int &NPROPS,
                          const floatVector &COORDS, const floatVector &DROT,   floatType &PNEWDT,         const floatType &CELENT,  const floatVector &DFGRD0,
                          const floatVector &DFGRD1, const int &NOEL,           const int &NPT,            const int &LAYER,         const int &KSPT,
                          const int *JSTEP,          const int &KINC ){

        return;
    }

}
