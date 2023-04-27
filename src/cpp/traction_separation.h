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

    errorOut computeCurrentDistanceGeneral( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                            const floatVector &F,    const floatVector &chi,  const floatVector &chiNL,
                                            floatVector &d );

    errorOut computeCurrentDistanceGeneral( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                            const floatVector &F,    const floatVector &chi,  const floatVector &chiNL,
                                            floatVector &d,
                                            floatMatrix &dddXi_1, floatMatrix &dddXi_2, floatMatrix &dddD,
                                            floatMatrix &dddF, floatMatrix &dddchi, floatMatrix &dddchiNL );

    errorOut computeCurrentDistanceGeneral( const floatVector &Xi_1, const floatVector &Xi_2, const floatVector &D,
                                            const floatVector &F,    const floatVector &chi,  const floatVector &chiNL,
                                            floatVector &d,
                                            floatMatrix &dddXi_1, floatMatrix &dddXi_2, floatMatrix &dddD,
                                            floatMatrix &dddF, floatMatrix &dddchi, floatMatrix &dddchiNL,
                                            floatMatrix &d2ddFdXi_1, floatMatrix &d2ddchidXi_1,
                                            floatMatrix &d2ddFdXi_2, floatMatrix &d2ddchiNLdXi_2,
                                            floatMatrix &d2ddFdD );

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

    errorOut computeLinearTraction( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                    const floatVector &parameters, floatVector &traction );

    errorOut computeLinearTraction( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                    const floatVector &parameters, floatVector &traction,
                                    floatMatrix &dtractionddn, floatMatrix &dtractionddt, floatMatrix &dtractiondp );

    errorOut computeLinearTraction( const floatVector &normalDeformationMeasure, const floatVector &tangentialDeformationMeasure,
                                    const floatVector &parameters, floatVector &traction,
                                    floatMatrix &dtractionddn, floatMatrix &dtractionddt, floatMatrix &dtractiondp,
                                    floatMatrix &d2tractionddndp, floatMatrix &d2tractionddtdp );

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

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &overlap );

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &overlap,
                                     floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                     floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdGradChi );

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &overlap,
                                     floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                     floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdGradChi,
                                     floatMatrix &d2OverlapdXi_1dXi_1, floatMatrix &d2OverlapdXi_1ddX, floatMatrix &d2OverlapdXi_1dR_nl, floatMatrix &d2OverlapdXi_1dF, floatMatrix &d2OverlapdXi_1dChi, floatMatrix &d2OverlapdXi_1dGradChi,
                                     floatMatrix &d2OverlapddXddX, floatMatrix &d2OverlapddXdR_nl, floatMatrix &d2OverlapddXdF, floatMatrix &d2OverlapddXdChi, floatMatrix &d2OverlapddXdGradChi,
                                     floatVector &d2OverlapdR_nldR_nl, floatMatrix &d2OverlapdR_nldF, floatMatrix &d2OverlapdR_nldChi, floatMatrix &d2OverlapdR_nldGradChi,
                                     floatMatrix &d2OverlapdFdF, floatMatrix &d2OverlapdFdChi, floatMatrix &d2OverlapdFdGradChi,
                                     floatMatrix &d2OverlapdChidChi, floatMatrix &d2OverlapdChidGradChi,
                                     floatMatrix &d2OverlapdGradChidGradChi );

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &gradChi,
                                     floatVector &overlap,
                                     floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                     floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdGradChi,
                                     floatMatrix &d2OverlapdXi_1dXi_1, floatMatrix &d2OverlapdXi_1ddX, floatMatrix &d2OverlapdXi_1dR_nl, floatMatrix &d2OverlapdXi_1dF, floatMatrix &d2OverlapdXi_1dChi, floatMatrix &d2OverlapdXi_1dGradChi,
                                     floatMatrix &d2OverlapddXddX, floatMatrix &d2OverlapddXdR_nl, floatMatrix &d2OverlapddXdF, floatMatrix &d2OverlapddXdChi, floatMatrix &d2OverlapddXdGradChi,
                                     floatVector &d2OverlapdR_nldR_nl, floatMatrix &d2OverlapdR_nldF, floatMatrix &d2OverlapdR_nldChi, floatMatrix &d2OverlapdR_nldGradChi,
                                     floatMatrix &d2OverlapdFdF, floatMatrix &d2OverlapdFdChi, floatMatrix &d2OverlapdFdGradChi,
                                     floatMatrix &d2OverlapdChidChi, floatMatrix &d2OverlapdChidGradChi,
                                     floatMatrix &d2OverlapdGradChidGradChi,
                                     floatMatrix &d3OverlapdXi_1dXi_1dXi_1, floatMatrix &d3OverlapdXi_1dXi_1ddX, floatMatrix &d3OverlapdXi_1dXi_1dR_nl, floatMatrix &d3OverlapdXi_1dXi_1dF, floatMatrix &d3OverlapdXi_1dXi_1dChi, floatMatrix &d3OverlapdXi_1dXi_1dGradChi,
                                     floatMatrix &d3OverlapdXi_1ddXddX, floatMatrix &d3OverlapdXi_1ddXdR_nl, floatMatrix &d3OverlapdXi_1ddXdF, floatMatrix &d3OverlapdXi_1ddXdChi, floatMatrix &d3OverlapdXi_1ddXdGradChi,
                                     floatMatrix &d3OverlapdXi_1dR_nldR_nl, floatMatrix &d3OverlapdXi_1dR_nldF, floatMatrix &d3OverlapdXi_1dR_nldChi, floatMatrix &d3OverlapdXi_1dR_nldGradChi,
                                     floatMatrix &d3OverlapdXi_1dFdF, floatMatrix &d3OverlapdXi_1dFdChi, floatMatrix &d3OverlapdXi_1dFdGradChi,
                                     floatMatrix &d3OverlapdXi_1dChidChi, floatMatrix &d3OverlapdXi_1dChidGradChi,
                                     floatMatrix &d3OverlapdXi_1dGradChidGradChi,
                                     floatMatrix &d3OverlapddXddXddX, floatMatrix &d3OverlapddXddXdR_nl, floatMatrix &d3OverlapddXddXdF, floatMatrix &d3OverlapddXddXdChi, floatMatrix &d3OverlapddXddXdGradChi,
                                     floatMatrix &d3OverlapddXdR_nldR_nl, floatMatrix &d3OverlapddXdR_nldF, floatMatrix &d3OverlapddXdR_nldChi, floatMatrix &d3OverlapddXdR_nldGradChi,
                                     floatMatrix &d3OverlapddXdFdF, floatMatrix &d3OverlapddXdFdChi, floatMatrix &d3OverlapddXdFdGradChi,
                                     floatMatrix &d3OverlapddXdChidChi, floatMatrix &d3OverlapddXdChidGradChi,
                                     floatMatrix &d3OverlapddXdGradChidGradChi,
                                     floatVector &d3OverlapdR_nldR_nldR_nl, floatMatrix &d3OverlapdR_nldR_nldF, floatMatrix &d3OverlapdR_nldR_nldChi, floatMatrix &d3OverlapdR_nldR_nldGradChi,
                                     floatMatrix &d3OverlapdR_nldFdF, floatMatrix &d3OverlapdR_nldFdChi, floatMatrix &d3OverlapdR_nldFdGradChi,
                                     floatMatrix &d3OverlapdR_nldChidChi, floatMatrix &d3OverlapdR_nldChidGradChi,
                                     floatMatrix &d3OverlapdR_nldGradChidGradChi,
                                     floatMatrix &d3OverlapdFdFdF, floatMatrix &d3OverlapdFdFdChi, floatMatrix &d3OverlapdFdFdGradChi,
                                     floatMatrix &d3OverlapdFdChidChi, floatMatrix &d3OverlapdFdChidGradChi,
                                     floatMatrix &d3OverlapdFdGradChidGradChi,
                                     floatMatrix &d3OverlapdChidChidChi, floatMatrix &d3OverlapdChidChidGradChi,
                                     floatMatrix &d3OverlapdChidGradChidGradChi,
                                     floatMatrix &d3OverlapdGradChidGradChidGradChi );

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                     const floatVector &F,    const floatVector &chi,  const floatVector &chi_nl_basis, const floatVector &gradChi,
                                     floatVector &overlap );

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                     const floatVector &F,    const floatVector &chi, const floatVector &chi_nl_basis, const floatVector &gradChi,
                                     floatVector &overlap,
                                     floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                     floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdChi_NL_B, floatMatrix &dOverlapdGradChi );

    errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                     const floatVector &F,    const floatVector &chi, const floatVector &chi_nl_basis, const floatVector &gradChi,
                                     floatVector &overlap,
                                     floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                     floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdChi_NL_B, floatMatrix &dOverlapdGradChi,
                                     floatMatrix &d2OverlapdXi_1dXi_1, floatMatrix &d2OverlapdXi_1ddX, floatMatrix &d2OverlapdXi_1dR_nl, floatMatrix &d2OverlapdXi_1dF, floatMatrix &d2OverlapdXi_1dChi, floatMatrix &d2OverlapdXi_1dChi_NL_B, floatMatrix &d2OverlapdXi_1dGradChi,
                                     floatMatrix &d2OverlapddXddX, floatMatrix &d2OverlapddXdR_nl, floatMatrix &d2OverlapddXdF, floatMatrix &d2OverlapddXdChi, floatMatrix &d2OerlapddXdChi_NL_B, floatMatrix &d2OverlapddXdGradChi,
                                     floatVector &d2OverlapdR_nldR_nl, floatMatrix &d2OverlapdR_nldF, floatMatrix &d2OverlapdR_nldChi, floatMatrix &d2OverlapdR_nldChi_NL_B, floatMatrix &d2OverlapdR_nldGradChi,
                                     floatMatrix &d2OverlapdFdF, floatMatrix &d2OverlapdFdChi, floatMatrix &d2OverlapdFdChi_NL_B, floatMatrix &d2OverlapdFdGradChi,
                                     floatMatrix &d2OverlapdChidChi, floatMatrix &d2OverlapdChidChi_NL_B, floatMatrix &d2OverlapdChidGradChi,
                                     floatMatrix &d2OverlapdChi_NL_BdChi_NL_B, floatMatrix &d2OverlapdChi_NL_BdGradChi,
                                     floatMatrix &d2OverlapdGradChidGradChi );

 errorOut computeParticleOverlap( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                  const floatVector &F,    const floatVector &chi, const floatVector &chi_nl_basis, const floatVector &gradChi,
                                  floatVector &overlap,
                                  floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                  floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdChi_NL_B, floatMatrix &dOverlapdGradChi,
                                  floatMatrix &d2OverlapdXi_1dXi_1, floatMatrix &d2OverlapdXi_1ddX, floatMatrix &d2OverlapdXi_1dR_nl, floatMatrix &d2OverlapdXi_1dF, floatMatrix &d2OverlapdXi_1dChi, floatMatrix &d2OverlapdXi_1dChi_NL_B, floatMatrix &d2OverlapdXi_1dGradChi,
                                  floatMatrix &d2OverlapddXddX, floatMatrix &d2OverlapddXdR_nl, floatMatrix &d2OverlapddXdF, floatMatrix &d2OverlapddXdChi, floatMatrix &d2OverlapddXdChi_NL_B, floatMatrix &d2OverlapddXdGradChi,
                                  floatVector &d2OverlapdR_nldR_nl, floatMatrix &d2OverlapdR_nldF, floatMatrix &d2OverlapdR_nldChi, floatMatrix &d2OverlapdR_nldChi_NL_B, floatMatrix &d2OverlapdR_nldGradChi,
                                  floatMatrix &d2OverlapdFdF, floatMatrix &d2OverlapdFdChi, floatMatrix &d2OverlapdFdChi_NL_B, floatMatrix &d2OverlapdFdGradChi,
                                  floatMatrix &d2OverlapdChidChi, floatMatrix &d2OverlapdChidChi_NL_B, floatMatrix &d2OverlapdChidGradChi,
                                  floatMatrix &d2OverlapdChi_NL_BdChi_NL_B, floatMatrix &d2OverlapdChi_NL_BdGradChi,
                                  floatMatrix &d2OverlapdGradChidGradChi,
                                  floatMatrix &d3OverlapdXi_1dXi_1dXi_1, floatMatrix &d3OverlapdXi_1dXi_1ddX, floatMatrix &d3OverlapdXi_1dXi_1dR_nl, floatMatrix &d3OverlapdXi_1dXi_1dF, floatMatrix &d3OverlapdXi_1dXi_1dChi, floatMatrix &d3OverlapdXi_1dXi_1dChi_NL_B, floatMatrix &d3OverlapdXi_1dXi_1dGradChi,
                                  floatMatrix &d3OverlapdXi_1ddXddX, floatMatrix &d3OverlapdXi_1ddXdR_nl, floatMatrix &d3OverlapdXi_1ddXdF, floatMatrix &d3OverlapdXi_1ddXdChi, floatMatrix &d3OverlapdXi_1ddXdChi_NL_B, floatMatrix &d3OverlapdXi_1ddXdGradChi,
                                  floatMatrix &d3OverlapdXi_1dR_nldR_nl, floatMatrix &d3OverlapdXi_1dR_nldF, floatMatrix &d3OverlapdXi_1dR_nldChi, floatMatrix &d3OverlapdXi_1dR_nldChi_NL_B, floatMatrix &d3OverlapdXi_1dR_nldGradChi,
                                  floatMatrix &d3OverlapdXi_1dFdF, floatMatrix &d3OverlapdXi_1dFdChi, floatMatrix &d3OverlapdXi_1dFdChi_NL_B, floatMatrix &d3OverlapdXi_1dFdGradChi,
                                  floatMatrix &d3OverlapdXi_1dChidChi, floatMatrix &d3OverlapdXi_1dChidChi_NL_B, floatMatrix &d3OverlapdXi_1dChidGradChi,
                                  floatMatrix &d3OverlapdXi_1dChi_NL_BdChi_NL_B, floatMatrix &d3OverlapdXi_1dChi_NL_BdGradChi,
                                  floatMatrix &d3OverlapdXi_1dGradChidGradChi,
                                  floatMatrix &d3OverlapddXddXddX, floatMatrix &d3OverlapddXddXdR_nl, floatMatrix &d3OverlapddXddXdF, floatMatrix &d3OverlapddXddXdChi, floatMatrix &d3OverlapddXddXdChi_NL_B, floatMatrix &d3OverlapddXddXdGradChi,
                                  floatMatrix &d3OverlapddXdR_nldR_nl, floatMatrix &d3OverlapddXdR_nldF, floatMatrix &d3OverlapddXdR_nldChi, floatMatrix &d3OverlapddXdR_nldChi_NL_B, floatMatrix &d3OverlapddXdR_nldGradChi,
                                  floatMatrix &d3OverlapddXdFdF, floatMatrix &d3OverlapddXdFdChi, floatMatrix &d3OverlapddXdFdChi_NL_B, floatMatrix &d3OverlapddXdFdGradChi,
                                  floatMatrix &d3OverlapddXdChidChi, floatMatrix &d3OverlapddXdChidChi_NL_B, floatMatrix &d3OverlapddXdChidGradChi,
                                  floatMatrix &d3OverlapddXdChi_NL_BdChi_NL_B, floatMatrix &d3OverlapddXdChi_NL_BdGradChi,
                                  floatMatrix &d3OverlapddXdGradChidGradChi,
                                  floatVector &d3OverlapdR_nldR_nldR_nl, floatMatrix &d3OverlapdR_nldR_nldF, floatMatrix &d3OverlapdR_nldR_nldChi, floatMatrix &d3OverlapdR_nldR_nldChi_NL_B, floatMatrix &d3OverlapdR_nldR_nldGradChi,
                                  floatMatrix &d3OverlapdR_nldFdF, floatMatrix &d3OverlapdR_nldFdChi, floatMatrix &d3OverlapdR_nldFdChi_NL_B, floatMatrix &d3OverlapdR_nldFdGradChi,
                                  floatMatrix &d3OverlapdR_nldChidChi, floatMatrix &d3OverlapdR_nldChidChi_NL_B, floatMatrix &d3OverlapdR_nldChidGradChi,
                                  floatMatrix &d3OverlapdR_nldChi_NL_BdChi_NL_B, floatMatrix &d3OverlapdR_nldChi_NL_BdGradChi,
                                  floatMatrix &d3OverlapdR_nldGradChidGradChi,
                                  floatMatrix &d3OverlapdFdFdF, floatMatrix &d3OverlapdFdFdChi, floatMatrix &d3OverlapdFdFdChi_NL_B, floatMatrix &d3OverlapdFdFdGradChi,
                                  floatMatrix &d3OverlapdFdChidChi, floatMatrix &d3OverlapdFdChidChi_NL_B, floatMatrix &d3OverlapdFdChidGradChi,
                                  floatMatrix &d3OverlapdFdChi_NL_BdChi_NL_B, floatMatrix &d3OverlapdFdChi_NL_BdGradChi,
                                  floatMatrix &d3OverlapdFdGradChidGradChi,
                                  floatMatrix &d3OverlapdChidChidChi, floatMatrix &d3OverlapdChidChidChi_NL_B, floatMatrix &d3OverlapdChidChidGradChi,
                                  floatMatrix &d3OverlapdChidChi_NL_BdChi_NL_B, floatMatrix &d3OverlapdChidChi_NL_BdGradChi,
                                  floatMatrix &d3OverlapdChidGradChidGradChi,
                                  floatMatrix &d3OverlapdChi_NL_BdChi_NL_BdChi_NL_B, floatMatrix &d3OverlapdChi_NL_BdChi_NL_BdGradChi,
                                  floatMatrix &d3OverlapdChi_NL_BdGradChidGradChi,
                                  floatMatrix &d3OverlapdGradChidGradChidGradChi );

    errorOut computeParticleOverlapChi_nl( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                           const floatVector &F,    const floatVector &chi,  const floatVector &chi_nl,
                                           floatVector &overlap );

    errorOut computeParticleOverlapChi_nl( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                           const floatVector &F,    const floatVector &chi,  const floatVector &chi_nl,
                                           floatVector &overlap,
                                           floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                           floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdChi_nl );

    errorOut computeParticleOverlapChi_nl( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                           const floatVector &F,    const floatVector &chi,  const floatVector &chi_nl,
                                           floatVector &overlap,
                                           floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                           floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdChi_nl,
                                           floatMatrix &d2OverlapdXi_1dXi_1, floatMatrix &d2OverlapdXi_1ddX, floatMatrix &d2OverlapdXi_1dR_nl, floatMatrix &d2OverlapdXi_1dF, floatMatrix &d2OverlapdXi_1dChi, floatMatrix &d2OverlapdXi_1dChi_nl,
                                           floatMatrix &d2OverlapddXddX, floatMatrix &d2OverlapddXdR_nl, floatMatrix &d2OverlapddXdF, floatMatrix &d2OverlapddXdChi, floatMatrix &d2OverlapddXdChi_nl,
                                           floatVector &d2OverlapdR_nldR_nl, floatMatrix &d2OverlapdR_nldF, floatMatrix &d2OverlapdR_nldChi, floatMatrix &d2OverlapdR_nldChi_nl,
                                           floatMatrix &d2OverlapdFdF, floatMatrix &d2OverlapdFdChi, floatMatrix &d2OverlapdFdChi_nl,
                                           floatMatrix &d2OverlapdChidChi, floatMatrix &d2OverlapdChidChi_nl,
                                           floatMatrix &d2OverlapdChi_nldChi_nl );
 
     errorOut computeParticleOverlapChi_nl( const floatVector &Xi_1, const floatVector &dX, const floatType &R_nl,
                                           const floatVector &F,    const floatVector &chi,  const floatVector &chi_nl,
                                           floatVector &overlap,
                                           floatMatrix &dOverlapdXi_1, floatMatrix &dOverlapddX, floatVector &dOverlapdR_nl,
                                           floatMatrix &dOverlapdF, floatMatrix &dOverlapdChi, floatMatrix &dOverlapdChi_nl,
                                           floatMatrix &d2OverlapdXi_1dXi_1, floatMatrix &d2OverlapdXi_1ddX, floatMatrix &d2OverlapdXi_1dR_nl, floatMatrix &d2OverlapdXi_1dF, floatMatrix &d2OverlapdXi_1dChi, floatMatrix &d2OverlapdXi_1dChi_nl,
                                           floatMatrix &d2OverlapddXddX, floatMatrix &d2OverlapddXdR_nl, floatMatrix &d2OverlapddXdF, floatMatrix &d2OverlapddXdChi, floatMatrix &d2OverlapddXdChi_nl,
                                           floatVector &d2OverlapdR_nldR_nl, floatMatrix &d2OverlapdR_nldF, floatMatrix &d2OverlapdR_nldChi, floatMatrix &d2OverlapdR_nldChi_nl,
                                           floatMatrix &d2OverlapdFdF, floatMatrix &d2OverlapdFdChi, floatMatrix &d2OverlapdFdChi_nl,
                                           floatMatrix &d2OverlapdChidChi, floatMatrix &d2OverlapdChidChi_nl,
                                           floatMatrix &d2OverlapdChi_nldChi_nl,
                                           floatMatrix &d3OverlapdXi_1dXi_1dXi_1, floatMatrix &d3OverlapdXi_1dXi_1ddX, floatMatrix &d3OverlapdXi_1dXi_1dR_nl, floatMatrix &d3OverlapdXi_1dXi_1dF, floatMatrix &d3OverlapdXi_1dXi_1dChi, floatMatrix &d3OverlapdXi_1dXi_1dChi_nl,
                                           floatMatrix &d3OverlapdXi_1ddXddX, floatMatrix &d3OverlapdXi_1ddXdR_nl, floatMatrix &d3OverlapdXi_1ddXdF, floatMatrix &d3OverlapdXi_1ddXdChi, floatMatrix &d3OverlapdXi_1ddXdChi_nl,
                                           floatMatrix &d3OverlapdXi_1dR_nldR_nl, floatMatrix &d3OverlapdXi_1dR_nldF, floatMatrix &d3OverlapdXi_1dR_nldChi, floatMatrix &d3OverlapdXi_1dR_nldChi_nl,
                                           floatMatrix &d3OverlapdXi_1dFdF, floatMatrix &d3OverlapdXi_1dFdChi, floatMatrix &d3OverlapdXi_1dFdChi_nl,
                                           floatMatrix &d3OverlapdXi_1dChidChi, floatMatrix &d3OverlapdXi_1dChidChi_nl,
                                           floatMatrix &d3OverlapdXi_1dChi_nldChi_nl,
                                           floatMatrix &d3OverlapddXddXddX, floatMatrix &d3OverlapddXddXdR_nl, floatMatrix &d3OverlapddXddXdF, floatMatrix &d3OverlapddXddXdChi, floatMatrix &d3OverlapddXddXdChi_nl,
                                           floatMatrix &d3OverlapddXdR_nldR_nl, floatMatrix &d3OverlapddXdR_nldF, floatMatrix &d3OverlapddXdR_nldChi, floatMatrix &d3OverlapddXdR_nldChi_nl,
                                           floatMatrix &d3OverlapddXdFdF, floatMatrix &d3OverlapddXdFdChi, floatMatrix &d3OverlapddXdFdChi_nl,
                                           floatMatrix &d3OverlapddXdChidChi, floatMatrix &d3OverlapddXdChidChi_nl,
                                           floatMatrix &d3OverlapddXdChi_nldChi_nl,
                                           floatVector &d3OverlapdR_nldR_nldR_nl, floatMatrix &d3OverlapdR_nldR_nldF, floatMatrix &d3OverlapdR_nldR_nldChi, floatMatrix &d3OverlapdR_nldR_nldChi_nl,
                                           floatMatrix &d3OverlapdR_nldFdF, floatMatrix &d3OverlapdR_nldFdChi, floatMatrix &d3OverlapdR_nldFdChi_nl,
                                           floatMatrix &d3OverlapdR_nldChidChi, floatMatrix &d3OverlapdR_nldChidChi_nl,
                                           floatMatrix &d3OverlapdR_nldChi_nldChi_nl,
                                           floatMatrix &d3OverlapdFdFdF, floatMatrix &d3OverlapdFdFdChi, floatMatrix &d3OverlapdFdFdChi_nl,
                                           floatMatrix &d3OverlapdFdChidChi, floatMatrix &d3OverlapdFdChidChi_nl,
                                           floatMatrix &d3OverlapdFdChi_nldChi_nl,
                                           floatMatrix &d3OverlapdChidChidChi, floatMatrix &d3OverlapdChidChidChi_nl,
                                           floatMatrix &d3OverlapdChidChi_nldChi_nl,
                                           floatMatrix &d3OverlapdChi_nldChi_nldChi_nl );

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
                                               floatVector &d3LdXdR_nldR_nl,
                                               floatVector &d4LdXdXdchi_nldchi_nl );

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

    errorOut solveOverlapDistance( const floatVector &chi_nl, const floatVector &xi_t, const floatType &R_nl, floatVector &d,
                                   floatMatrix &dddchi_nl, floatMatrix &dddxi_t, floatVector &dddR_nl,
                                   floatMatrix &d2ddchi_nldchi_nl, floatMatrix &d2ddchi_nldxi_t, floatMatrix &d2ddchi_nldR_nl,
                                   floatMatrix &d2ddxi_tdxi_t, floatMatrix &d2ddxi_tdR_nl,
                                   floatVector &d2ddR_nldR_nl,
                                   floatMatrix &d3ddchi_nldchi_nldchi_nl, floatMatrix &d3ddchi_nldchi_nldxi_t, floatMatrix &d3ddchi_nldchi_nldR_nl,
                                   floatMatrix &d3ddchi_nldxi_tdxi_t, floatMatrix &d3ddchi_nldxi_tdR_nl,
                                   floatMatrix &d3ddchi_nldR_nldR_nl,
                                   floatMatrix &d3ddxi_tdxi_tdxi_t, floatMatrix &d3ddxi_tdxi_tdR_nl,
                                   floatMatrix &d3ddxi_tdR_nldR_nl,
                                   floatVector &d3ddR_nldR_nldR_nl,
                                   const floatType tolr = 1e-9, const floatType tola = 1e-9, const unsigned int max_iteration = 20,
                                   const unsigned int max_ls = 5, const floatType alpha_ls = 1e-4 );

}

#endif
