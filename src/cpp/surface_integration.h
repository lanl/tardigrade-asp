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

    errorOut decomposeSphere( const floatType &radius, const unsigned int &elementCount,
                              floatVector &points, std::vector<unsigned int> &connectivity );

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

    errorOut formBaseCubePoints( const unsigned int &elementCount, floatVector &points );

    errorOut formCubeConnectivity( const unsigned int &elementCount, std::vector< unsigned int > &connectivity );

    errorOut evaluateQuadraticShapeFunctions( const floatType &xi, const floatType &eta, floatVector &shapeFunctions );

    errorOut evaluateGradQuadraticShapeFunctions( const floatType &xi, const floatType &eta, floatMatrix &gradShapeFunctions );

    errorOut interpolateFunction( const floatType &xi, const floatType &eta, floatMatrix &nodalValues, floatVector &answer );

    errorOut localGradientFunction( const floatType &xi, const floatType &eta, floatMatrix &nodalValues, floatMatrix &answer );

    errorOut localJacobian( const floatType &xi, const floatType &eta, floatMatrix &nodalPositions, floatType &jacobian );

}

#endif
