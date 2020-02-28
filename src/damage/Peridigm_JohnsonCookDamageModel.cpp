/*! \file Peridigm_JohnsonCookDamageModel.cpp */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include "Peridigm_JohnsonCookDamageModel.hpp"
#include "Peridigm_Field.hpp"
#include <algorithm>    // std::max
#include "material_utilities.h" // to use Influence Function 
#include "correspondence.h" // to compute weighted volume

using namespace std;

PeridigmNS::JohnsonCookDamageModel::JohnsonCookDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params), 
    m_thresholdDamage(0.0), 
    m_criticalDamage(0.0), 
    m_damageParameter1(0.0),
    m_damageParameter2(0.0),
    m_damageParameter3(0.0),
    m_criticalDamageToNeglectMaterialPoint(0.95),
    m_damageFieldId(-1), 
    m_flyingPointFlagFieldId(-1),
    m_bondDamageFieldId(-1),
    m_bondLevelEquivalentPlasticStrainFieldId(-1),
    m_bondLevelStressTriaxialityFieldId(-1),
    m_bondLevelDamageParameterFieldId(-1),
    m_brokenBondVolumeAveragedFieldId(-1),
    m_horizonFieldId(-1), 
    m_volumeFieldId(-1),
    m_coordinatesFieldId(-1), 
    m_jacobianDeterminantFieldId(-1),
    m_undamagedWeightedVolumeFieldId(-1)
{

  // Johnson-Cook Correspondence Model
  m_damageParameter1 = params.get<double>("Damage Parameter 1");
  m_damageParameter2 = params.get<double>("Damage Parameter 2");
  m_damageParameter3 = params.get<double>("Damage Parameter 3");
  m_criticalDamage = params.get<double>("Critical Damage");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(1.0 <= m_criticalDamage, "**** Error: Critical Damage should be less than 1.0");
  m_thresholdDamage = params.get<double>("Threshold Damage");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(m_criticalDamage <= m_thresholdDamage, "**** Error: Threshold Damage should be less than Critical Damage");

  if(params.isParameter("Critical Damage To Neglect Material Point")) 
    m_criticalDamageToNeglectMaterialPoint = params.get<double>("Critical Damage To Neglect Material Point"); // default value is 95%

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_damageFieldId                                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_flyingPointFlagFieldId                         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Flying_Point_Flag");
  m_bondDamageFieldId                              = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_brokenBondVolumeAveragedFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Broken_Bond_Volume_Averaged");
  m_horizonFieldId                                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId                                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_coordinatesFieldId                             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_jacobianDeterminantFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Jacobian_Determinant");
  m_undamagedWeightedVolumeFieldId                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Undamaged_Weighted_Volume");

  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_flyingPointFlagFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_brokenBondVolumeAveragedFieldId);
  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_jacobianDeterminantFieldId);
  m_fieldIds.push_back(m_undamagedWeightedVolumeFieldId);

  m_bondLevelEquivalentPlasticStrainFieldId        = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Equivalent_Plastic_Strain");
  m_bondLevelStressTriaxialityFieldId              = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Stress_Triaxiality");
  m_bondLevelDamageParameterFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage_Parameter");
  m_fieldIds.push_back(m_bondLevelEquivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_bondLevelStressTriaxialityFieldId);
  m_fieldIds.push_back(m_bondLevelDamageParameterFieldId);
}

PeridigmNS::JohnsonCookDamageModel::~JohnsonCookDamageModel()
{
}

void
PeridigmNS::JohnsonCookDamageModel::initialize(const double dt,
                                               const int numOwnedPoints,
                                               const int* ownedIDs,
                                               const int* neighborhoodList,
                                               PeridigmNS::DataManager& dataManager) const
{
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  dataManager.getData(m_bondLevelDamageParameterFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelDamageParameterFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
}

void
PeridigmNS::JohnsonCookDamageModel::computeDamage(const double dt,
                                                  const int numOwnedPoints,
                                                  const int* ownedIDs,
                                                  const int* neighborhoodList,
                                                  PeridigmNS::DataManager& dataManager) const
{
  //NOTE: this fixes the inter-processor and inter-block double counting problem
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *damageN, *damageNP1, *brokenBondVolumeAveragedNP1;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_N)->ExtractView(&damageN);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damageNP1);
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_NP1)->ExtractView(&brokenBondVolumeAveragedNP1);

  double *horizon, *volume, *coordinates;
  double *undamagedWeightedVolume, *jacobianDeterminant;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_undamagedWeightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&undamagedWeightedVolume);
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_N)->ExtractView(&jacobianDeterminant);

  double* coord = coordinates;
  double* neighborCoord;
  double* delta = horizon;

  double deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  double neighborVolume, omega, scalarTemp;

  double deqps, epsF, deltaDamage;

  const int *neighborListPtr = neighborhoodList;
  int numNeighbors, neighborIndex;

  // Zero out the bond damage to compute intact weighted volume (to
  // eventually evaluate damage of material points)
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_NP1)->PutScalar(-1.0);

  double *bondDamage, *flyingPointFlag;
  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_NP1)->ExtractView(&flyingPointFlag);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  // Set the bond damage to the previous value
  *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N));

  // the mixed state-bond based formulation is considered else.
  // Update the bond damage directly based on its own properties
  double *equivalentPlasticStrainN, *equivalentPlasticStrainNP1, *stressTriaxiality, *damageParameterN, *damageParameterNP1, *brokenBondVolumeAveragedN;
  dataManager.getData(m_bondLevelEquivalentPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&equivalentPlasticStrainN);
  dataManager.getData(m_bondLevelEquivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&equivalentPlasticStrainNP1);
  dataManager.getData(m_bondLevelStressTriaxialityFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressTriaxiality);
  dataManager.getData(m_bondLevelDamageParameterFieldId, PeridigmField::STEP_N)->ExtractView(&damageParameterN);
  dataManager.getData(m_bondLevelDamageParameterFieldId, PeridigmField::STEP_NP1)->ExtractView(&damageParameterNP1);
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_N)->ExtractView(&brokenBondVolumeAveragedN);

  for(int iID=0 ; iID<numOwnedPoints ; ++iID, ++damageN, ++damageNP1,
        ++flyingPointFlag, ++delta, ++undamagedWeightedVolume, coord+=3, 
        ++brokenBondVolumeAveragedN, ++brokenBondVolumeAveragedNP1)
  {
    if(*flyingPointFlag < 0.0){
      // damage = volume-average of broken bonds
      *damageNP1 = 0.0;
      *brokenBondVolumeAveragedNP1 = 0.0;

      // Loop over the neighbors and compute bond damage
      numNeighbors = *neighborListPtr; neighborListPtr++;

      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
            bondDamage++, equivalentPlasticStrainN++, 
            equivalentPlasticStrainNP1++, stressTriaxiality++,
            damageParameterN++, damageParameterNP1++){

        neighborIndex = *neighborListPtr;

        // D is a non-decreasing function 
        // D_dot = eqps_dot / eps_f
        // eps_f = d_1 + d_2 exp(-d_3 * sigma_m/sigma_e)
        *damageParameterNP1 = *damageParameterN;

        // Damage model is evaluated before material model.
        // The updated values from the previous time step have been swapped
        // between N and NP1 states.
        deqps = *equivalentPlasticStrainN - *equivalentPlasticStrainNP1;

        epsF = m_damageParameter1 + m_damageParameter2 * exp(-m_damageParameter3 * *stressTriaxiality);

        deltaDamage = deqps / epsF;
        if(deltaDamage > 0.0)
          *damageParameterNP1 += deltaDamage; 

        if(*damageParameterNP1 >= 1.0){
          *damageParameterNP1 = 1.0;
        }

        if(*damageParameterNP1 > m_thresholdDamage){
          *bondDamage = (*damageParameterNP1 - m_thresholdDamage) / (1.0 - m_thresholdDamage);
        }
        else{
          *bondDamage = 0.0;
        }
        
        neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
        neighborCoord = coordinates + 3*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

        *damageNP1 += *damageParameterNP1 * omega / *undamagedWeightedVolume * neighborVolume;
        *brokenBondVolumeAveragedNP1 += *bondDamage * omega / *undamagedWeightedVolume * neighborVolume;
      }

      if(*brokenBondVolumeAveragedNP1 > m_criticalDamageToNeglectMaterialPoint) // to prevent cases that only a couple of bonds are left, causing instability issues. 
        *flyingPointFlag = 1.0;
    }
    else{
      // for the visualization purposes
      *damageNP1 = *damageN;
      *brokenBondVolumeAveragedNP1 = *brokenBondVolumeAveragedN;

      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamage += numNeighbors; 
      equivalentPlasticStrainN += numNeighbors; 
      equivalentPlasticStrainNP1 += numNeighbors; 
      stressTriaxiality += numNeighbors; 
      damageParameterN += numNeighbors; 
      damageParameterNP1 += numNeighbors; 
    }
  }
}
