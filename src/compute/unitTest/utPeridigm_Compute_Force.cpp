/*! \file utPeridigm_Compute_Force.cpp */

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
#include "../Peridigm_Compute_Force.hpp"
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

using namespace boost::unit_test;

//! Create dataManager object with manufactured data, then check data filled by Compute_Force class 
Teuchos::RCP<PeridigmNS::DataManager> createDataManager(Teuchos::RCP<Teuchos::ParameterList> discretizationParams, Teuchos::RCP<PeridigmNS::Compute_Force> computeForce) {
  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = Teuchos::rcp(new Epetra_SerialComm);
  #endif

  // Initialize Discretization (initialize maps)
 
  PeridigmNS::DiscretizationFactory discFactory(discretizationParams);
  Teuchos::RCP<PeridigmNS::AbstractDiscretization> peridigmDisc = discFactory.create(comm);

  // oneDimensionalMap
   Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap = peridigmDisc->getMap(1);

  // oneDimensionalOverlapMap (includes ghosts)
  Teuchos::RCP<const Epetra_BlockMap> oneDimensionalOverlapMap = peridigmDisc->getOverlapMap(1);

  // threeDimensionalMap
  Teuchos::RCP<const Epetra_BlockMap> threeDimensionalMap = peridigmDisc->getMap(3);

  // threeDimensionalOverlapMap (includes ghosts)
  Teuchos::RCP<const Epetra_BlockMap> threeDimensionalOverlapMap = peridigmDisc->getOverlapMap(3);

  // bondConstitutiveDataMap (non-overlapping map)
  Teuchos::RCP<const Epetra_BlockMap> bondMap = peridigmDisc->getBondMap();

  // Set the initial positions
  Teuchos::RCP<Epetra_Vector> x = peridigmDisc->getInitialX();

  // Create the importers
  Teuchos::RCP<const Epetra_Import> oneDimensionalMapToOneDimensionalOverlapMapImporter = Teuchos::rcp(new Epetra_Import(*oneDimensionalOverlapMap, *oneDimensionalMap));
  Teuchos::RCP<const Epetra_Import> threeDimensionalMapToThreeDimensionalOverlapMapImporter = Teuchos::rcp(new Epetra_Import(*threeDimensionalOverlapMap, *threeDimensionalMap));

  // get the neighborlist from the discretization
   Teuchos::RCP<const PeridigmNS::NeighborhoodData> neighborhoodData = peridigmDisc->getNeighborhoodData();

  // Step #3: Initialize data manager
  Teuchos::RCP<PeridigmNS::DataManager> dataManager = Teuchos::rcp(new PeridigmNS::DataManager);
  dataManager->setMaps(oneDimensionalMap, threeDimensionalMap, oneDimensionalOverlapMap, threeDimensionalOverlapMap, bondMap);

  // Create a master list of variable specs
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > variableSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);

  // Fill list with specs utilized by Peridigm object
  variableSpecs->push_back(Field_NS::VOLUME);
  variableSpecs->push_back(Field_NS::COORD3D);
  variableSpecs->push_back(Field_NS::DISPL3D);
  variableSpecs->push_back(Field_NS::CURCOORD3D);
  variableSpecs->push_back(Field_NS::VELOC3D);
  variableSpecs->push_back(Field_NS::FORCE_DENSITY3D);
  variableSpecs->push_back(Field_NS::CONTACT_FORCE_DENSITY3D);

  // Don't add variable specs requested materials -- there are no material models used in this unit test

  // Now add the variable specs requested by the compute class
  std::vector<Field_NS::FieldSpec> computeSpecs = computeForce->getFieldSpecs();
  for (unsigned int i=0; i < computeSpecs.size(); i++) {
     variableSpecs->push_back(computeSpecs[i]);
  }

  // Remove duplicates
  std::unique(variableSpecs->begin(), variableSpecs->end());

  // Allocate data in the dataManager
  dataManager->allocateData(variableSpecs);

  // Fill the dataManager with data from the discretization
  dataManager->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->Import(*(peridigmDisc->getCellVolume()), *oneDimensionalMapToOneDimensionalOverlapMapImporter, Insert);
  dataManager->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->Import(*x, *threeDimensionalMapToThreeDimensionalOverlapMapImporter, Insert);
  dataManager->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_N)->Import(*x, *threeDimensionalMapToThreeDimensionalOverlapMapImporter, Insert);
  dataManager->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->Import(*x, *threeDimensionalMapToThreeDimensionalOverlapMapImporter, Insert);

  return dataManager;

}

void FourPointTest() {

  // Set up discretization parameterlist
  Teuchos::RCP<Teuchos::ParameterList> discretizationParams = rcp(new Teuchos::ParameterList("Discretization"));
  discretizationParams->set("Type", "PdQuickGrid");
  discretizationParams->set("Horizon", 5.0);
  // pdQuickGrid tensor product mesh generator parameters
  Teuchos::ParameterList& pdQuickGridParams = discretizationParams->sublist("TensorProduct3DMeshGenerator");
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

  // Create Compute_Force object
  // hand in NULL for parent pointer since not used by compute class
  Teuchos::RCP<PeridigmNS::Compute_Force> computeForce = Teuchos::rcp(new PeridigmNS::Compute_Force(NULL));

  // Create the data manager
  Teuchos::RCP<PeridigmNS::DataManager> dataManager = createDataManager(discretizationParams,computeForce);

  // Access the data we need
  Teuchos::RCP<Epetra_Vector> force, force_density, volume;
  force_density = dataManager->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1);
  force         = dataManager->getData(Field_NS::FORCE3D, Field_ENUM::STEP_NP1);
  volume        = dataManager->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE);

  // Manufacture force density data
  double *force_density_values  = force_density->Values();
  int numElements = volume->Map().NumMyElements();
  for (int i=0;i<numElements;i++) {
    force_density_values[3*i] = 3.0*i;
    force_density_values[3*i+1] = (3.0*i)+1.0;
    force_density_values[3*i+2] = (3.0*i)+2.0;
  }

  // Call the compute class
  int retval = computeForce->compute(dataManager);
  BOOST_CHECK_EQUAL( retval, 0 );

  // Now check that volumes and forces are correct
  double *volume_values = volume->Values();
  double *force_values  = force->Values();
  for (int i=0;i<numElements;i++)
    BOOST_CHECK_CLOSE(volume_values[i], 1.5, 1.0e-15);
  for (int i=0;i<numElements;i++) {
    BOOST_CHECK_CLOSE(force_values[3*i],   1.5*(3.0*i),   1.0e-15);
    BOOST_CHECK_CLOSE(force_values[3*i+1], 1.5*((3.0*i)+1.0), 1.0e-15);
    BOOST_CHECK_CLOSE(force_values[3*i+2], 1.5*((3.0*i)+2.0), 1.0e-15);
  }

}


bool init_unit_test_suite() {
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPeridigm_Compute_Force");
  proc->add(BOOST_TEST_CASE(&FourPointTest));
  framework::master_test_suite().add(proc);

  return success;
}

bool init_unit_test() {
  return init_unit_test_suite();
}

int main (int argc, char* argv[]) {
  int numProcs = 1;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

  int returnCode = -1;
  
  if(numProcs >= 1 && numProcs <= 4){
    returnCode = unit_test_main(init_unit_test, argc, argv);
  }
  else{
    std::cerr << "Unit test runtime ERROR: utPeridigm_Compute_Force only makes sense on 1 to 4 processors." << std::endl;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return returnCode;
}
