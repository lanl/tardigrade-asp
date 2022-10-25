/**
  ******************************************************************************
  * \file traction_separation.h
  ******************************************************************************
  * The header file for an implementation of various traction separation laws
  * for use in ASP. We will start with linear laws and, if required, progress
  * from there.
  ******************************************************************************
  */

#ifndef TRACTIONSEPARATION_H
#define TRACTIONSEPARATION_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<constitutive_tools.h>

namespace tractionSeparation{

    typedef constitutiveTools::errorNode errorNode; //!< Redefinition for the error node
    typedef constitutiveTools::errorOut errorOut; //!< Redefinition for a pointer to the error node
    typedef constitutiveTools::floatType floatType; //!< Define the float values type.
    typedef constitutiveTools::floatVector floatVector; //!< Define a vector of floats
    typedef constitutiveTools::floatMatrix floatMatrix; //!< Define a matrix of floats

    errorOut computeCurrentDistance( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &d );

    errorOut computeCurrentDistance( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &d,
                                     floatMatrix &dddXi_1, floatMatrix &dddXi_2, floatMatrix &dddD,
                                     floatMatrix &dddF, floatMatrix &dddChi, floatMatrix &dddGradChi );

    errorOut computeCurrentDistance( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &d,
                                     floatMatrix &dddXi_1, floatMatrix &dddXi_2, floatMatrix &dddD,
                                     floatMatrix &dddF, floatMatrix &dddChi, floatMatrix &dddGradChi,
                                     floatMatrix &d2ddFdXi_1,       floatMatrix &d2ddFdXi_2,       floatMatrix &d2ddFdD,
                                     floatMatrix &d2ddChidXi_1,     floatMatrix &d2ddChidXi_2,     floatMatrix &d2ddChidD,
                                     floatMatrix &d2ddGradChidXi_1, floatMatrix &d2ddGradChidXi_2, floatMatrix &d2ddGradChidD );

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt );

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt,
                              floatMatrix &ddndd, floatMatrix &ddndn,
                              floatMatrix &ddtdd, floatMatrix &ddtdn );

    errorOut decomposeVector( const floatVector &d, const floatVector &n,
                              floatVector &dn, floatVector &dt,
                              floatMatrix &ddndd, floatMatrix &ddndn,
                              floatMatrix &ddtdd, floatMatrix &ddtdn,
                              floatMatrix &d2dndddd, floatMatrix &d2dndddn,
                              floatMatrix &d2dndndn,
                              floatMatrix &d2dtdddd, floatMatrix &d2dtdddn,
                              floatMatrix &d2dtdndn );

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy );

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt );

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt,
                                          floatVector &d2energyddnddn, floatVector &d2energyddnddt,
                                          floatVector &d2energyddtddt );

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt, floatVector &denergydParameters );

    errorOut computeLinearTractionEnergy( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                          const floatVector &parameters, floatType &energy,
                                          floatVector &denergyddn, floatVector &denergyddt, floatVector &denergydParameters,
                                          floatVector &d2energyddnddn, floatVector &d2energyddnddt, floatVector &d2energyddndParameters,
                                          floatVector &d2energyddtddt, floatVector &d2energyddtdParameters,
                                          floatVector &d2energydParametersdParameters );

    errorOut computeNansonsRelation( const floatVector &deformationGradient, const floatVector &dAN, floatVector &dan );

    errorOut computeNansonsRelation( const floatVector &deformationGradient, const floatVector &dAN, floatVector &dan,
                                     floatMatrix &ddandF, floatMatrix &ddanddAN );

    errorOut computeNansonsRelation( const floatVector &deformationGradient, const floatVector &dAN, floatVector &dan,
                                     floatMatrix &ddandF, floatMatrix &ddanddAN,
                                     floatMatrix &d2dandFdF, floatMatrix &d2dandFddAN );

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatType &R_2, const floatVector &D,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &overlap );

    errorOut computeOverlapDistanceLagrangian( const floatVector &X, const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatType &L );

    errorOut computeOverlapDistanceLagrangian( const floatVector &X, const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatType &L,
                                               floatVector &dLdX, floatVector &dLdchi_nl, floatVector &dLdxi_t, floatType &dLdR_nl );

    errorOut computeOverlapDistanceLagrangian( const floatVector &X, const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatType &L,
                                               floatVector &dLdX, floatVector &dLdchi_nl, floatVector &dLdxi_t, floatType &dLdR_nl,
                                               floatVector &d2LdXdX, floatVector &d2LdXdchi_nl, floatVector &d2LdXdxi_t, floatVector &d2LdXdR_nl,
                                               floatVector &d2Ldchi_nldchi_nl, floatVector &d2Ldchi_nldxi_t, floatVector &d2Ldchi_nldR_nl,
                                               floatVector &d2Ldxi_tdxi_t, floatVector &d2Ldxi_tdR_nl,
                                               floatType &d2LdR_nldR_nl );

    errorOut computeOverlapDistanceLagrangian( const floatVector &X, const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatType &L,
                                               floatVector &dLdX, floatVector &dLdchi_nl, floatVector &dLdxi_t, floatType &dLdR_nl,
                                               floatVector &d2LdXdX, floatVector &d2LdXdchi_nl, floatVector &d2LdXdxi_t, floatVector &d2LdXdR_nl,
                                               floatVector &d2Ldchi_nldchi_nl, floatVector &d2Ldchi_nldxi_t, floatVector &d2Ldchi_nldR_nl,
                                               floatVector &d2Ldxi_tdxi_t, floatVector &d2Ldxi_tdR_nl,
                                               floatType &d2LdR_nldR_nl,
                                               floatVector &d3LdXdXdX, floatVector &d3LdXdXdchi_nl, floatVector &d3LdXdchi_nldchi_nl,
                                               floatVector &d3LdXdXdxi_t, floatVector &d3LdXdchi_nldxi_t,
                                               floatVector &d3LdXdXdR_nl, floatVector &d3LdXdchi_nldR_nl,
                                               floatVector &d3LdXdxi_tdxi_t, floatVector &d3LdXdxi_tdR_nl,
                                               floatVector &d3LdXdR_nldR_nl );

    errorOut solveOverlapDistance( const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatVector &d,
                                   const floatType tolr = 1e-9, const floatType tola = 1e-9, const unsigned int max_iteration = 20,
                                   const unsigned int max_ls = 5, const floatType alpha_ls = 1e-4 );

    errorOut solveOverlapDistance( const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatVector &d,
                                   floatMatrix &dddchi_nl, floatMatrix &dddxi_t, floatVector &dddR_nl,
                                   const floatType tolr = 1e-9, const floatType tola = 1e-9, const unsigned int max_iteration = 20,
                                   const unsigned int max_ls = 5, const floatType alpha_ls = 1e-4 );

    errorOut solveOverlapDistance( const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatVector &d,
                                   floatMatrix &dddchi_nl, floatMatrix &dddxi_t, floatVector &dddR_nl,
                                   floatMatrix &d2ddchi_nldchi_nl, floatMatrix &d2ddchi_nldxi_t, floatMatrix &d2ddchi_nldR_nl,
                                   floatMatrix &d2ddxi_tdxi_t, floatMatrix &d2ddxi_tdR_nl,
                                   floatVector &d2ddR_nldR_nl,
                                   const floatType tolr = 1e-9, const floatType tola = 1e-9, const unsigned int max_iteration = 20,
                                   const unsigned int max_ls = 5, const floatType alpha_ls = 1e-4 );

}

#endif
