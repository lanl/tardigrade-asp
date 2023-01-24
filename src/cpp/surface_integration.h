/**
  ******************************************************************************
  * \file surfaceIntegration.h
  ******************************************************************************
  * The header file for an implementation of various functions required for
  * surface integrations for use in ASP.
  ******************************************************************************
  */

#ifndef SURFACEINTEGRATION_H
#define SURFACEINTEGRATION_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<constitutive_tools.h>

namespace surfaceIntegration{

    errorOut decomposeSphere( const floatType &radius, unsigned int &elementCount,
                              floatMatrix &points, std:vector<unsigned int> &connectivity );

    template <class T>
    errorOut integrateProperty( const floatMatrix &points, const std::vector<unsigned int> &connectivity,
                                const std::vector< T > &nodalData, T &answer );

}
