/*! \file utPeridigm_Strain_Energy.cpp */

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

#include <Peridigm_AbstractDiscretization.hpp>
#include "../Peridigm_Compute_Strain_Energy.hpp"
#include "../Peridigm_Compute_Local_Strain_Energy.hpp"
#include "../Peridigm_Compute_Global_Strain_Energy.hpp"
#include <Peridigm_DataManager.hpp>
#include <Peridigm_DiscretizationFactory.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#include <Epetra_Import.h>
#include <Teuchos_ParameterList.hpp>
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif
#include <vector>
#include "../../core/Peridigm.hpp"

using namespace boost::unit_test;
using namespace Teuchos;

Teuchos::RCP<PeridigmNS::Peridigm> createFourPointModel() {
  Teuchos::RCP<Epetra_Comm> comm;
#ifdef HAVE_MPI
  comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

  // set up parameter lists
  // these data would normally be read from an input xml file
  Teuchos::RCP<Teuchos::ParameterList> peridigmParams = rcp(new Teuchos::ParameterList());

  // material parameters
  Teuchos::ParameterList& materialParams = peridigmParams->sublist("Materials");
  Teuchos::ParameterList& linearElasticMaterialParams = materialParams.sublist("My Elastic Material");
  linearElasticMaterialParams.set("Material Model", "Elastic");
  linearElasticMaterialParams.set("Density", 7800.0);
  linearElasticMaterialParams.set("Bulk Modulus", 130.0e9);
  linearElasticMaterialParams.set("Shear Modulus", 78.0e9);

  // blocks
  Teuchos::ParameterList& blockParams = peridigmParams->sublist("Blocks");
  Teuchos::ParameterList& blockOneParams = blockParams.sublist("My Group of Blocks");
  blockOneParams.set("Block Names", "block_1");
  blockOneParams.set("Material", "My Elastic Material");

  // Set up discretization parameterlist
  Teuchos::ParameterList& discretizationParams = peridigmParams->sublist("Discretization");
  discretizationParams.set("Type", "PdQuickGrid");
  discretizationParams.set("Horizon", 5.0);
  // pdQuickGrid tensor product mesh generator parameters
  Teuchos::ParameterList& pdQuickGridParams = discretizationParams.sublist("TensorProduct3DMeshGenerator");
  pdQuickGridParams.set("Type", "PdQuickGrid");
  pdQuickGridParams.set("X Origin",  0.0);
  pdQuickGridParams.set("Y Origin",  0.0);
  pdQuickGridParams.set("Z Origin",  0.0);
  pdQuickGridParams.set("X Length",  6.0);
  pdQuickGridParams.set("Y Length",  1.0);
  pdQuickGridParams.set("Z Length",  1.0);
  pdQuickGridParams.set("Number Points X", 4);
  pdQuickGridParams.set("Number Points Y", 1);
  pdQuickGridParams.set("Number Points Z", 1);

  // output parameters (to force instantiation of data storage for compute classes in DataManager)
  Teuchos::ParameterList& outputParams = peridigmParams->sublist("Output");
  Teuchos::ParameterList& outputFields = outputParams.sublist("Output Variables");
  outputFields.set("Strain_Energy", true);
  outputFields.set("Global_Strain_Energy", true);

  // create the Peridigm object
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = Teuchos::rcp(new PeridigmNS::Peridigm(comm, peridigmParams));

  return peridigm;
}

void FourPointTest() 
{
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createFourPointModel();

  // Get the neighborhood data
  PeridigmNS::NeighborhoodData neighborhoodData = (*peridigm->getGlobalNeighborhoodData());
  // Get the data manager
  Teuchos::RCP<PeridigmNS::DataManager> dataManager = peridigm->getBlocks()->begin()->getDataManager();
  // Access the data we need
  Teuchos::RCP<Epetra_Vector> volume, ref, coords, dilatation, strain_energy;
  volume                = dataManager->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE);
  ref                   = dataManager->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE);
  coords                = dataManager->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1);
  dilatation            = dataManager->getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1);
  strain_energy         = dataManager->getData(Field_NS::STRAIN_ENERGY, Field_ENUM::STEP_NP1);
  // Get the neighborhood structure
  const int numOwnedPoints = (neighborhoodData.NumOwnedPoints());

  // Manufacture velocity data
  double *ref_values = ref->Values();
  double *coords_values = coords->Values();
  double *dilatation_values = dilatation->Values();
  int *myGIDs = ref->Map().MyGlobalElements();
  int numElements = numOwnedPoints;
  int numTotalElements = volume->Map().NumMyElements();
  for (int i=0;i<numTotalElements;i++) {
    int ID = myGIDs[i];
    coords_values[3*i] = ref_values[3*i] + 3.0*ID;
    coords_values[3*i+1] = ref_values[3*i+1] + 3.0*ID+1.0;
    coords_values[3*i+2] = ref_values[3*i+2] + 3.0*ID+2.0;
    dilatation_values[i] = 9.369;
  }

  // Get the blocks
  Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks = peridigm->getBlocks();

  // Fire the compute classes to fill the strain energy data
  peridigm->getComputeManager()->compute(blocks);
  
  // Now check that volumes and energy is correct
  double *volume_values = volume->Values();
  double *strain_energy_values  = strain_energy->Values();
  double globalSE = blocks->begin()->getScalarData(Field_NS::GLOBAL_STRAIN_ENERGY);
  BOOST_CHECK_CLOSE(globalSE, 4.0*8.559e12, 0.01);	// Check global scalar value
  for (int i=0;i<numElements;i++)
    BOOST_CHECK_CLOSE(volume_values[i], 1.5, 1.0e-15);
  for (int i=0;i<numElements;i++) {
    BOOST_CHECK_CLOSE(strain_energy_values[i],   8.559e12, 0.01);
  }
}

bool init_unit_test_suite() 
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPeridigm_Compute_Strain_Energy");
  proc->add(BOOST_TEST_CASE(&FourPointTest));
  framework::master_test_suite().add(proc);

  return success;
}

bool init_unit_test() 
{
  return init_unit_test_suite();
}

int main (int argc, char* argv[]) 
{
  int numProcs = 1;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

  int returnCode = -1;
  
  if(numProcs >= 1 && numProcs <= 4) {
    returnCode = unit_test_main(init_unit_test, argc, argv);
  }
  else {
    std::cerr << "Unit test runtime ERROR: utPeridigm_Compute_Strain_Energy only makes sense on 1 to 4 processors." << std::endl;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return returnCode;
}
