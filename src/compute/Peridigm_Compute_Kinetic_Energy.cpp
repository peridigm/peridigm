/*! \file Peridigm_Compute_Kinetic_Energy.cpp */

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

#include "Peridigm_Compute_Kinetic_Energy.hpp"
#include "../core/Peridigm.hpp"

//! Standard constructor.
PeridigmNS::Compute_Kinetic_Energy::Compute_Kinetic_Energy(PeridigmNS::Peridigm *peridigm_ ){peridigm = peridigm_;}

//! Destructor.
PeridigmNS::Compute_Kinetic_Energy::~Compute_Kinetic_Energy(){}


//! Returns the fieldspecs computed by this class
std::vector<Field_NS::FieldSpec> PeridigmNS::Compute_Kinetic_Energy::getFieldSpecs() const 
{
  	std::vector<Field_NS::FieldSpec> myFieldSpecs;

  	return myFieldSpecs;
}

//! Perform computation
int PeridigmNS::Compute_Kinetic_Energy::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const
{
	return 0;
}

int PeridigmNS::Compute_Kinetic_Energy::computeKineticEnergy( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, bool storeLocal ) const
{ 
	int retval;
	
	double globalKE = 0.0;
	Teuchos::RCP<Epetra_Vector> velocity, volume, force, numNeighbors, neighborID, kinetic_energy;
	std::vector<PeridigmNS::Block>::iterator blockIt;
	for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
	{
		Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
		Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
		const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
		const int* ownedIDs = neighborhoodData->OwnedIDs();

		velocity              = dataManager->getData(Field_NS::VELOC3D, Field_ENUM::STEP_NP1);
		volume                = dataManager->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE);
		kinetic_energy        = dataManager->getData(Field_NS::KINETIC_ENERGY, Field_ENUM::STEP_NP1);
		
		// Sanity check
		if (velocity->Map().NumMyElements() != volume->Map().NumMyElements())
		{
			retval = 1;
			return(retval);
		}
 	
		// Collect values
		double *volume_values = volume->Values();
		double *velocity_values = velocity->Values();
		double *kinetic_energy_values = kinetic_energy->Values();
	
		// Get the material properties 
		double density  = blockIt->getMaterialModel()->Density();
		//double SM = blockIt->getMaterialModel()->ShearModulus();
		//double BM = blockIt->getMaterialModel()->BulkModulus();	

		
		// Initialize global kinetic energy value
		double KE = 0.0;

		// volume is a scalar and force a vector, so maps are different; must do multiplication on per-element basis
		int numElements = numOwnedPoints;

		double vol;

		for (int i=0;i<numElements;i++) 
		{
			int ID = ownedIDs[i];
			vol = volume_values[ID];
			double v1 = velocity_values[3*ID];
			double v2 = velocity_values[3*ID+1];
			double v3 = velocity_values[3*ID+2];
			
			if (storeLocal)
				kinetic_energy_values[i] = 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);	//Store local kinetic energy	
			else
                                KE = KE + 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);      // Update the global kinetic energy
			
		}

		if (!storeLocal)
		{
			// Update info across processors
			double localKE, globalBlockKE;
			localKE = KE;

			peridigm->getEpetraComm()->SumAll(&localKE, &globalBlockKE, 1);

			globalKE += globalBlockKE;
		}
	}
	
	// Store global energy in block (block globals are static, so only need to assign data to first block)
	if (!storeLocal)
		blocks->begin()->getScalarData(Field_NS::GLOBAL_KINETIC_ENERGY) = globalKE;

	return(0);

}
