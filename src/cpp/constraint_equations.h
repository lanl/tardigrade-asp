/**
  ******************************************************************************
  * \file constraint_equations.h
  ******************************************************************************
  * The header file for an implementation of various constraint equations.
  ******************************************************************************
  */

#ifndef CONSTRAINTEQUATIONS_H
#define CONSTRAINTEQUATIONS_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<constitutive_tools.h>

namespace constraintEquations{

    typedef constitutiveTools::errorNode errorNode; //!< Redefinition for the error node
    typedef constitutiveTools::errorOut errorOut; //!< Redefinition for a pointer to the error node
    typedef constitutiveTools::floatType floatType; //!< Define the float values type.
    typedef constitutiveTools::floatVector floatVector; //!< Define a vector of floats
    typedef constitutiveTools::floatMatrix floatMatrix; //!< Define a matrix of floats

    errorOut tractionConstraint( const floatVector &cauchyStress, const floatVector &n, const floatVector &traction, const floatType &P, floatType &C );

    errorOut tractionConstraint( const floatVector &cauchyStress, const floatVector &n, const floatVector &traction, const floatType &P, floatType &C,
                                 floatVector &dCdCauchyStress, floatVector &dCdNormal, floatVector &dCdTraction, floatType &dCdP );

    errorOut tractionConstraint( const floatVector &cauchyStress, const floatVector &n, const floatVector &traction, const floatType &P, floatType &C,
                                 floatVector &dCdCauchyStress, floatVector &dCdNormal, floatVector &dCdTraction, floatType &dCdP,
                                 floatVector &d2CdCauchyStressdNormal, floatVector &d2CdCauchyStressdP,
                                 floatVector &d2CdNormaldP,            floatVector &d2CdTractiondP );

}

#endif
