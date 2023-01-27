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

    typedef constitutiveTools::errorNode errorNode; //!< Redefinition for the error node
    typedef constitutiveTools::errorOut errorOut; //!< Redefinition for a pointer to the error node
    typedef constitutiveTools::floatType floatType; //!< Define the float values type.
    typedef constitutiveTools::floatVector floatVector; //!< Define a vector of floats
    typedef constitutiveTools::floatMatrix floatMatrix; //!< Define a matrix of floats

    errorOut decomposeSphere( const floatType &radius, unsigned int &elementCount,
                              floatMatrix &points, std::vector<unsigned int> &connectivity );

    errorOut buildSurfacePoints( const floatType &x0, const floatType &y0, const floatType &z0,
                                 const floatType &dx, const floatType &dy,
                                 const unsigned int n_points_x, const unsigned int n_points_y,
                                 floatVector &points );

    errorOut rotatePoints( const floatVector &points,
                           const floatType &thetaX, const floatType &thetaY, const floatType &thetaZ,
                           floatVector &rotatedPoints );

    errorOut formSurfaceConnectivity( const std::vector< unsigned int > &surfaceIDs,
                                      const unsigned int &n_elements_x, const unsigned int &n_elements_y,
                                      unsigned int &index, std::vector< unsigned int > &connectivity );

    template <class T>
    errorOut integrateProperty( const floatVector &points, const std::vector<unsigned int> &connectivity,
                                const std::vector< T > &nodalData, T &answer );

}

#endif
