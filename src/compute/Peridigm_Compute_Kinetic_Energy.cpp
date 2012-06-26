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
  	myFieldSpecs.push_back(Field_NS::GLOBAL_KINETIC_ENERGY);
	myFieldSpecs.push_back(Field_NS::KINETIC_ENERGY);

  	return myFieldSpecs;
}



//! Fill the energy vectors
int PeridigmNS::Compute_Kinetic_Energy::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const
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
		double *kinetic_energy_values  = kinetic_energy->Values();
	
		// \todo Generalize this for multiple materials
		// Get the material properties 
		double density  = peridigm->getMaterialModels()->operator[](0)->Density();
		//double SM = peridigm->getMaterialModels()->operator[](0)->ShearModulus();
		//double BM = peridigm->getMaterialModels()->operator[](0)->BulkModulus();	

		// Initialize kinetic energy value
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
			
			// Update the global kinetic energy
			KE = KE + 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);
			//Store local kinetic energy
			kinetic_energy_values[i] = 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);	
		}

		// Update info across processors
		double localKE, globalBlockKE;
		localKE = KE;

		peridigm->getEpetraComm()->SumAll(&localKE, &globalBlockKE, 1);

/*
        	if (peridigm->getEpetraComm()->MyPID() == 1)
        	{
			std::cout << std::endl;	
			std::cout << "ref coords for node 0 = (" << ref_values[3*0] << ", " << ref_values[3*0+1] << ", " << ref_values[3*0+2] << ")" << "\n";
			std::cout << "ref coords for node 1 = (" << ref_values[3*1] << ", " << ref_values[3*1+1] << ", " << ref_values[3*1+2] << ")" << "\n";
			std::cout << "ref coords for node 2 = (" << ref_values[3*2] << ", " << ref_values[3*2+1] << ", " << ref_values[3*2+2] << ")" << "\n";
			std::cout << "ref coords for node 3 = (" << ref_values[3*3] << ", " << ref_values[3*3+1] << ", " << ref_values[3*3+2] << ")" << "\n" << "\n";
			std::cout << "cur coords for node 0 = (" << coord_values[3*0] << ", " << coord_values[3*0+1] << ", " << coord_values[3*0+2] << ")" << "\n";
                	std::cout << "cur coords for node 1 = (" << coord_values[3*1] << ", " << coord_values[3*1+1] << ", " << coord_values[3*1+2] << ")" << "\n";
                	std::cout << "cur coords for node 2 = (" << coord_values[3*2] << ", " << coord_values[3*2+1] << ", " << coord_values[3*2+2] << ")" << "\n";
                	std::cout << "cur coords for node 3 = (" << coord_values[3*3] << ", " << coord_values[3*3+1] << ", " << coord_values[3*3+2] << ")" << "\n" << "\n";
                	std::cout << "dilatation = (" << dilatation_values[0] << ", " << dilatation_values[1] << ", " << dilatation_values[2] << ", " << dilatation_values[3] << ")" << "\n" << "\n";
			std::cout << "weighted volume = (" << w_volume_values[0] << ", " << w_volume_values[1] << ", " << w_volume_values[2] << ", " << w_volume_values[3] << ")" << "\n" << "\n";

			std::cout << "Hello!" << std::endl;
			std::cout << std::endl;
			std::cout << "Total Kinetic Energy  =  " << globalKE << std::endl; 
        		std::cout << std::endl;
        		std::cout << "Total Strain Energy  =  " << globalSE << std::endl;
		}
*/
		globalKE += globalBlockKE;
	}
	
	// Store global energy in block (block globals are static, so only need to assign data to first block)
	blocks->begin()->getScalarData(Field_NS::GLOBAL_KINETIC_ENERGY) = globalKE;

	return(0);

}
