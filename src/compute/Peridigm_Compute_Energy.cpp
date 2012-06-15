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
#include "../core/Peridigm.hpp"

//! Standard constructor.
PeridigmNS::Compute_Energy::Compute_Energy(PeridigmNS::Peridigm *peridigm_ ){peridigm = peridigm_;}

//! Destructor.
PeridigmNS::Compute_Energy::~Compute_Energy(){}


//! Returns the fieldspecs computed by this class
std::vector<Field_NS::FieldSpec> PeridigmNS::Compute_Energy::getFieldSpecs() const 
{
  	std::vector<Field_NS::FieldSpec> myFieldSpecs;
  	myFieldSpecs.push_back(Field_NS::KINETIC_ENERGY);
	myFieldSpecs.push_back(Field_NS::STRAIN_ENERGY);
	myFieldSpecs.push_back(Field_NS::STRAIN_ENERGY_DENSITY);


  	return myFieldSpecs;
}



//! Fill the energy vectors
int PeridigmNS::Compute_Energy::compute(const int numOwnedPoints,
                                                  const int* ownedIDs,
                                                  const int* neighborhoodList,
                                                  PeridigmNS::DataManager& dataManager) const 
{
  
	int retval;

  	Teuchos::RCP<Epetra_Vector> velocity, volume, force, ref, coord, w_volume, dilatation, numNeighbors, neighborID, kinetic_energy, strain_energy, strain_energy_density;
  	velocity              = dataManager.getData(Field_NS::VELOC3D, Field_ENUM::STEP_NP1);
  	volume                = dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE);
        ref                   = dataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE);
        coord                 = dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1);
	w_volume              = dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE);
	dilatation            = dataManager.getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1);
	kinetic_energy        = dataManager.getData(Field_NS::KINETIC_ENERGY, Field_ENUM::STEP_NP1);
        strain_energy         = dataManager.getData(Field_NS::STRAIN_ENERGY, Field_ENUM::STEP_NP1);
        strain_energy_density = dataManager.getData(Field_NS::STRAIN_ENERGY_DENSITY, Field_ENUM::STEP_NP1);
	

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
	
	// \todo Generalize this for multiple materials and add accessor function for Shear Modulus
	// Get the material properties 
	double density  = peridigm->getMaterialModels()->operator[](0)->Density();
        double SM = peridigm->getMaterialModels()->operator[](0)->ShearModulus();
        double BM = peridigm->getMaterialModels()->operator[](0)->BulkModulus();	

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
		W = W + 0.5*BM*dilatation_values[ID]*dilatation_values[ID] + 0.5*(15*SM/w_vol)*We;		
		strain_energy_density_values[i] = 0.5*BM*dilatation_values[ID]*dilatation_values[ID] + 0.5*(15*SM/w_vol)*We;		
		
		//Update the strain energy
		strain_energy_values[i] = vol*strain_energy_density_values[i];
		SE = SE + strain_energy_values[i];
                // Update the global kinetic energy
    		KE = KE + 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);
  		kinetic_energy_values[i] = 0.5*vol*density*(v1*v1 + v2*v2 + v3*v3);	
	}

	// Update info across processors
	double localKE, localSE, globalKE, globalSE;
        localKE = KE;
        localSE = SE;

	peridigm->getEpetraComm()->SumAll(&localKE, &globalKE, 1);
	peridigm->getEpetraComm()->SumAll(&localSE, &globalSE, 1);


        if (peridigm->getEpetraComm()->MyPID() == 1)
        {
/*		std::cout << std::endl;	
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
*/
	}
	

	return(0);

}
