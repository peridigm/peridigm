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

#include <Peridigm_DataManager.hpp>
#include <Peridigm_Discretization.hpp>
#include "../Peridigm_Compute_Force.hpp"
#include <Peridigm_DiscretizationFactory.hpp>
#include "Peridigm_Field.hpp"
#include "Peridigm.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <vector>

#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif

using namespace PeridigmNS;

Teuchos::RCP<Peridigm> createFourPointModel(Teuchos::RCP<Epetra_Comm> comm) {
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
  blockOneParams.set("Horizon", 5.0);

  // Set up discretization parameterlist
  Teuchos::ParameterList& discretizationParams = peridigmParams->sublist("Discretization");
  discretizationParams.set("Type", "PdQuickGrid");

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
  Teuchos::RCP<Discretization> nullDiscretization;
  Teuchos::RCP<Peridigm> peridigm = Teuchos::rcp(new Peridigm(comm, peridigmParams, nullDiscretization));

  return peridigm;
}

// *************
// FourPointTest
// *************
TEUCHOS_UNIT_TEST(Compute_Force, FourPointTest) {

  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

  int numProcs = comm->NumProc();

  TEST_COMPARE(numProcs, <=, 4);
  if(numProcs > 4){
    std::cerr << "Unit test runtime ERROR: utPeridigm_Compute_Force only makes sense on 1 to 4 processors." << std::endl;
    return;
  }

  Teuchos::RCP<Peridigm> peridigm = createFourPointModel(comm);

  FieldManager& fieldManager = FieldManager::self();

  // Access the data we need
  Teuchos::RCP<Epetra_Vector> force, force_density, volume;
  force_density = peridigm->getBlocks()->begin()->getData(fieldManager.getFieldId("Force_Density"), PeridigmField::STEP_NP1);
  force         = peridigm->getBlocks()->begin()->getData(fieldManager.getFieldId("Force"), PeridigmField::STEP_NP1);
  volume        = peridigm->getBlocks()->begin()->getData(fieldManager.getFieldId("Volume"), PeridigmField::STEP_NONE);

  // Manufacture force density data
  double *force_density_values  = force_density->Values();
  int numElements = volume->Map().NumMyElements();
  for (int i=0;i<numElements;i++) {
    force_density_values[3*i] = 3.0*i;
    force_density_values[3*i+1] = (3.0*i)+1.0;
    force_density_values[3*i+2] = (3.0*i)+2.0;
  }

  // Create Compute_Force object
  Teuchos::RCP<Teuchos::ParameterList> params;
  Teuchos::RCP<Teuchos::ParameterList> computeClassGlobalData;
  Teuchos::RCP<Compute_Force> computeForce = Teuchos::rcp(new Compute_Force( params, peridigm->getEpetraComm(), computeClassGlobalData ));

  // Get the blocks
  Teuchos::RCP< std::vector<Block> > blocks = peridigm->getBlocks();

  // Call the compute class
  int retval = computeForce->compute( blocks );
  TEST_EQUALITY_CONST( retval, 0 );

  // Now check that volumes and forces are correct
  double *volume_values = volume->Values();
  double *force_values  = force->Values();
  for (int i=0;i<numElements;i++)
    TEST_FLOATING_EQUALITY(volume_values[i], 1.5, 1.0e-15);
  for (int i=0;i<numElements;i++) {
    TEST_FLOATING_EQUALITY(force_values[3*i], 1.5*(3.0*i), 1.0e-15); 
    TEST_FLOATING_EQUALITY(force_values[3*i+1], 1.5*((3.0*i)+1.0), 1.0e-15); 
    TEST_FLOATING_EQUALITY(force_values[3*i+2], 1.5*((3.0*i)+2.0), 1.0e-15); 
  }

}
