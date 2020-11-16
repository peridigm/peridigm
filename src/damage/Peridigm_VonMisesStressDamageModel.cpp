/*! \file Peridigm_VonMisesStressDamageModel.cpp */

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

#include "Peridigm_VonMisesStressDamageModel.hpp"
#include "Peridigm_Field.hpp"
#include "material_utilities.h" // to use Influence Function 
#include "correspondence.h" // to compute weighted volume

using namespace std;

PeridigmNS::VonMisesStressDamageModel::VonMisesStressDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params), 
    m_thresholdDamage(0.99),
    m_criticalDamage(1.00), 
    m_criticalVonMisesStress(0.0),
    m_criticalDamageToNeglectMaterialPoint(0.95),
    m_flyingPointFlagFieldId(-1),
    m_bondDamageFieldId(-1),
    m_bondLevelVonMisesStressFieldId(-1),
    m_brokenBondVolumeAveragedFieldId(-1),
    m_horizonFieldId(-1), 
    m_volumeFieldId(-1),
    m_coordinatesFieldId(-1), 
    m_jacobianDeterminantFieldId(-1),
    m_undamagedWeightedVolumeFieldId(-1)
{

  // Johnson-Cook Correspondence Model
  m_criticalVonMisesStress = params.get<double>("Critical Von Mises Stress");

  if(params.isParameter("Critical Damage"))
    m_criticalDamage = params.get<double>("Critical Damage");
  if(params.isParameter("Threshold Damage"))
    m_thresholdDamage = params.get<double>("Threshold Damage");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(m_criticalDamage < m_thresholdDamage, "**** Error: Threshold Damage cannot be greater than Critical Damage");

  if(params.isParameter("Critical Damage To Neglect Material Point")) 
    m_criticalDamageToNeglectMaterialPoint = params.get<double>("Critical Damage To Neglect Material Point"); // default value is 95%

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_flyingPointFlagFieldId                         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Flying_Point_Flag");
  m_bondDamageFieldId                              = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_bondLevelVonMisesStressFieldId                 = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Von_Mises_Stress");
  m_brokenBondVolumeAveragedFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Broken_Bond_Volume_Averaged");
  m_horizonFieldId                                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId                                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_coordinatesFieldId                             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_jacobianDeterminantFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Jacobian_Determinant");
  m_undamagedWeightedVolumeFieldId                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Undamaged_Weighted_Volume");

  m_fieldIds.push_back(m_flyingPointFlagFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_bondLevelVonMisesStressFieldId);
  m_fieldIds.push_back(m_brokenBondVolumeAveragedFieldId);
  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_jacobianDeterminantFieldId);
  m_fieldIds.push_back(m_undamagedWeightedVolumeFieldId);
}

PeridigmNS::VonMisesStressDamageModel::~VonMisesStressDamageModel()
{
}

void
PeridigmNS::VonMisesStressDamageModel::initialize(const double dt,
                                                  const int numOwnedPoints,
                                                  const int* ownedIDs,
                                                  const int* neighborhoodList,
                                                  PeridigmNS::DataManager& dataManager) const
{
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
}

void
PeridigmNS::VonMisesStressDamageModel::computeDamage(const double dt,
                                                     const int numOwnedPoints,
                                                     const int* ownedIDs,
                                                     const int* neighborhoodList,
                                                     PeridigmNS::DataManager& dataManager) const
{

  double *vonMisesStress;
  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&vonMisesStress);

  double *brokenBondVolumeAveragedN, *brokenBondVolumeAveragedNP1;
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_NP1)->ExtractView(&brokenBondVolumeAveragedNP1);
  dataManager.getData(m_brokenBondVolumeAveragedFieldId, PeridigmField::STEP_N)->ExtractView(&brokenBondVolumeAveragedN);

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

  double tempDmg;

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

  for(int iID=0 ; iID<numOwnedPoints ; ++iID, ++flyingPointFlag, 
      ++delta, ++undamagedWeightedVolume, coord+=3, 
      ++brokenBondVolumeAveragedN, ++brokenBondVolumeAveragedNP1)
  {
    *brokenBondVolumeAveragedNP1 = 0.0;

    // Loop over the neighbors and compute bond damage
    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
      bondDamage++, vonMisesStress++){

      neighborIndex = *neighborListPtr;
      
      if(*vonMisesStress / m_criticalVonMisesStress >= m_criticalDamage)
        tempDmg  = 1.0;
      else if (*vonMisesStress / m_criticalVonMisesStress <= m_thresholdDamage)
        tempDmg = 0.0;
      else
        tempDmg = (*vonMisesStress / m_criticalVonMisesStress - m_thresholdDamage) / (m_criticalDamage - m_thresholdDamage);
      
      *bondDamage = max(*bondDamage, tempDmg); // non-decreasing damage 
      
      // Compute the weighted-volume broken bonds
      neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
      neighborCoord = coordinates + 3*neighborIndex;

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                deformedBondY*deformedBondY +
                                deformedBondZ*deformedBondZ);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

      *brokenBondVolumeAveragedNP1 += *bondDamage * omega / *undamagedWeightedVolume * neighborVolume;
    }

    if(*brokenBondVolumeAveragedNP1 > m_criticalDamageToNeglectMaterialPoint) // to prevent cases that only a couple of bonds are left, causing instability issues. 
      *flyingPointFlag = 1.0;
  }
}
