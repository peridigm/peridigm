/*! \file Peridigm_Compute_Energy.cpp */

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
// Michael L. Parks      parks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include <vector>

#include "Peridigm_Compute_Energy.hpp"
#include "Peridigm_Field.hpp"

//! Standard constructor.
PeridigmNS::Compute_Energy::Compute_Energy(Teuchos::RCP<const Teuchos::ParameterList> params,
                                           Teuchos::RCP<const Epetra_Comm> epetraComm_)
  : Compute(params, epetraComm_), m_volumeFieldId(-1), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_velocityFieldId(-1),
    m_weightedVolumeFieldId(-1), m_dilatationFieldId(-1), m_kineticEnergyFieldId(-1), m_strainEnergyDensityFieldId(-1), m_strainEnergyFieldId(-1),
    m_globalKineticEnergyFieldId(-1), m_globalStrainEnergyFieldId(-1)
{
  FieldManager& fieldManager = FieldManager::self();
  m_volumeFieldId = fieldManager.getFieldId("Volume");
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  m_velocityFieldId = fieldManager.getFieldId("Velocity");
  m_weightedVolumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
  m_dilatationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
  m_kineticEnergyFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Kinetic_Energy");
  m_strainEnergyDensityFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Strain_Energy_Density");
  m_strainEnergyFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Strain_Energy");
  m_globalKineticEnergyFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Global_Kinetic_Energy");
  m_globalStrainEnergyFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Global_Strain_Energy");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocityFieldId);
  m_fieldIds.push_back(m_weightedVolumeFieldId);
  m_fieldIds.push_back(m_dilatationFieldId);
  m_fieldIds.push_back(m_kineticEnergyFieldId);
  m_fieldIds.push_back(m_strainEnergyDensityFieldId);
  m_fieldIds.push_back(m_strainEnergyFieldId);
  m_fieldIds.push_back(m_globalKineticEnergyFieldId);
  m_fieldIds.push_back(m_globalStrainEnergyFieldId);
}

//! Destructor.
PeridigmNS::Compute_Energy::~Compute_Energy(){}

//! Fill the energy vectors
int PeridigmNS::Compute_Energy::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const
{
  int retval;
  double globalKE, globalSE, globalSEDensity;
  globalKE = globalSE = globalSEDensity = 0.0;
  Teuchos::RCP<Epetra_Vector> velocity, volume, force, ref, coord, w_volume, dilatation, numNeighbors, neighborID, kinetic_energy, strain_energy, strain_energy_density;
  std::vector<Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    Teuchos::RCP<NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    volume                = blockIt->getData(m_volumeFieldId, PeridigmField::STEP_NONE);
    ref                   = blockIt->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE);
    coord                 = blockIt->getData(m_coordinatesFieldId, PeridigmField::STEP_NP1);
    velocity              = blockIt->getData(m_velocityFieldId, PeridigmField::STEP_NP1);
    w_volume              = blockIt->getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE);
    dilatation            = blockIt->getData(m_dilatationFieldId, PeridigmField::STEP_NP1);
    kinetic_energy        = blockIt->getData(m_kineticEnergyFieldId, PeridigmField::STEP_NONE);
    strain_energy_density = blockIt->getData(m_strainEnergyDensityFieldId, PeridigmField::STEP_NONE);
    strain_energy         = blockIt->getData(m_strainEnergyFieldId, PeridigmField::STEP_NONE);
	
    // Sanity check
    if (velocity->Map().NumMyElements() != volume->Map().NumMyElements() || velocity->Map().NumMyElements() != ref->Map().NumMyElements())
      {
        retval = 1;
        return(retval);
      }
 	
    // Collect values
    double *volume_values = volume->Values();
    double *velocity_values = velocity->Values();
    double *ref_values = ref->Values();
    double *coord_values = coord->Values();
    double *w_volume_values = w_volume->Values();
    double *dilatation_values = dilatation->Values();
    double *kinetic_energy_values  = kinetic_energy->Values();
    double *strain_energy_values  = strain_energy->Values();
    double *strain_energy_density_values  = strain_energy_density->Values();
	
    // Get the material properties 
    double density  = blockIt->getMaterialModel()->Density();
    double SM = blockIt->getMaterialModel()->ShearModulus();
    double BM = blockIt->getMaterialModel()->BulkModulus();	

    // Initialize energy values
    double KE,  SE;
    KE = SE = 0.0;

    // Initialize strain energy density
    double W = 0.0; 

    // volume is a scalar and force a vector, so maps are different; must do multiplication on per-element basis
    int numElements = numOwnedPoints;

    double vol, vol2, w_vol;

    // Initialize local strain energy density
    double We;

    int neighborhoodListIndex = 0;
    for (int i=0;i<numElements;i++) 
      {
        int numNeighbors = neighborhoodList[neighborhoodListIndex++];
        int ID = ownedIDs[i];
        We = 0.0;
        vol = volume_values[ID];
        w_vol = w_volume_values[ID];
        double v1 = velocity_values[3*ID];
        double v2 = velocity_values[3*ID+1];
        double v3 = velocity_values[3*ID+2];
        for (int j=0; j<numNeighbors; j++)
          {
            int neighborID = neighborhoodList[neighborhoodListIndex++];
            TEUCHOS_TEST_FOR_EXCEPT_MSG(neighborID < 0, "Invalid neighbor list\n");
            int Ne = neighborID;
            vol2 = volume_values[Ne];
            double psi1 = (ref_values[3*Ne] - ref_values[3*i]);
            double psi2 = (ref_values[3*Ne+1] - ref_values[3*i+1]);
            double psi3 = (ref_values[3*Ne+2] - ref_values[3*i+2]);
            double eta1 = (coord_values[3*Ne] - coord_values[3*i]);
            double eta2 = (coord_values[3*Ne+1] - coord_values[3*i+1]);
            double eta3 = (coord_values[3*Ne+2] - coord_values[3*i+2]);
            
            // Compute the reference position
            double x = sqrt(psi1*psi1 + psi2*psi2 + psi3*psi3);
            
            // Compute the deformed position
            double y = sqrt(eta1*eta1 + eta2*eta2 + eta3*eta3);
            
            // Compute the extension
            double e = y-x;
            // \todo Generalize for different influence functions
            // Update the local strain energy density
            We = We + (1.0)*(e - dilatation_values[ID]*x/3)*vol2;
          }
        // Update the strain energy density                
        W = W + 0.5*BM*dilatation_values[ID]*dilatation_values[ID] + 0.5*(15.0*SM/w_vol)*We;		
        strain_energy_density_values[i] = 0.5*BM*dilatation_values[ID]*dilatation_values[ID] + 0.5*(15.0*SM/w_vol)*We;		
		
        // Update the strain energy
        strain_energy_values[i] = vol*strain_energy_density_values[i];
        SE = SE + strain_energy_values[i];
        // Update the global kinetic energy
        KE = KE + 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);
        kinetic_energy_values[i] = 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);	
      }
    
    // Update info across processors
    double localKE, localSE, globalBlockKE, globalBlockSE;
    localKE = KE;
    localSE = SE;
    
    epetraComm->SumAll(&localKE, &globalBlockKE, 1);
    epetraComm->SumAll(&localSE, &globalBlockSE, 1);

    globalKE += globalBlockKE;
    globalSE += globalBlockSE;
	}

  // Store global values
  Teuchos::RCP<Epetra_Vector> data;
  data = blocks->begin()->getData(m_globalKineticEnergyFieldId, PeridigmField::STEP_NONE);
  (*data)[0] = globalKE;
  data = blocks->begin()->getData(m_globalStrainEnergyFieldId, PeridigmField::STEP_NONE);
  (*data)[0] = globalSE;

  return(0);
}
