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
#include "../../core/Peridigm.hpp"

using namespace boost::unit_test;

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
  Teuchos::ParameterList& linearElasticMaterialParams = materialParams.sublist("My Linear Elastic Material");
  linearElasticMaterialParams.set("Material Model", "Linear Elastic");
  linearElasticMaterialParams.set("Density", 7800.0);
  linearElasticMaterialParams.set("Bulk Modulus", 130.0e9);
  linearElasticMaterialParams.set("Shear Modulus", 78.0e9);

  // blocks
  Teuchos::ParameterList& blockParams = peridigmParams->sublist("Blocks");
  Teuchos::ParameterList& blockOneParams = blockParams.sublist("My Group of Blocks");
  blockOneParams.set("Block Names", "block_1");
  blockOneParams.set("Material", "My Linear Elastic Material");

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
  outputFields.set("Force", true);

  // create the Peridigm object
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = Teuchos::rcp(new PeridigmNS::Peridigm(comm, peridigmParams));

  return peridigm;
}

void FourPointTest() {

  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createFourPointModel();

  // Get the data manager
  Teuchos::RCP<PeridigmNS::DataManager> dataManager = peridigm->getBlocks()->begin()->getDataManager();

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

  // Create Compute_Force object
  Teuchos::RCP<PeridigmNS::Compute_Force> computeForce = Teuchos::rcp(new PeridigmNS::Compute_Force( &(*peridigm) ) );

  // Get the blocks
  Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks = peridigm->getBlocks();

  // Call the compute class
  int retval = computeForce->compute( blocks );
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
