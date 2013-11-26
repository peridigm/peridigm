/*! \file utPeridigm_PdQuickGridDiscretization_MPI_np2.cpp */

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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_GlobalMPISession.hpp"
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#include <Epetra_MpiComm.h>
#include "Peridigm_PdQuickGridDiscretization.hpp"
#include "Peridigm_HorizonManager.hpp"


using namespace Teuchos;
using namespace PeridigmNS;


TEUCHOS_UNIT_TEST(PdQuickGridDiscretization_MPI_np2, SimpleTensorProductMeshTest) {


  Teuchos::RCP<Epetra_Comm> comm;
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));


  int numProcs = comm->NumProc();
  int rank     = comm->MyPID();

  TEST_COMPARE(numProcs, ==, 2);

  if(numProcs != 2){
     std::cerr << "Unit test runtime ERROR: utPeridigm_PdQuickGridDiscretization_MPI_np2 only makes sense on 2 processors" << std::endl;
     return;
  }

  RCP<ParameterList> discParams = rcp(new ParameterList);

  // create a 2x2x2 discretization
  // specify a spherical neighbor search with the horizon a tad longer than the mesh spacing
  discParams->set("Type", "PdQuickGrid");
  discParams->set("NeighborhoodType", "Spherical");
  ParameterList& quickGridParams = discParams->sublist("TensorProduct3DMeshGenerator");
  quickGridParams.set("Type", "PdQuickGrid");
  quickGridParams.set("X Origin", 0.0);
  quickGridParams.set("Y Origin", 0.0);
  quickGridParams.set("Z Origin", 0.0);
  quickGridParams.set("X Length", 1.0);
  quickGridParams.set("Y Length", 1.0);
  quickGridParams.set("Z Length", 1.0);
  quickGridParams.set("Number Points X", 2);
  quickGridParams.set("Number Points Y", 2);
  quickGridParams.set("Number Points Z", 2);

  // initialize the horizon manager and set the horizon to 0.501
  ParameterList blockParameterList;
  ParameterList& blockParams = blockParameterList.sublist("My Block");
  blockParams.set("Block Names", "block_1");
  blockParams.set("Horizon", 0.501);
  PeridigmNS::HorizonManager::self().loadHorizonInformationFromBlockParameters(blockParameterList);

  // create the discretization
  RCP<PdQuickGridDiscretization> discretization =
    rcp(new PdQuickGridDiscretization(comm, discParams));

  // sanity check, calling with a dimension other than 1 or 3 should throw an exception
  TEST_THROW(discretization->getGlobalOwnedMap(0), Teuchos::Exceptions::InvalidParameter);
  TEST_THROW(discretization->getGlobalOwnedMap(2), Teuchos::Exceptions::InvalidParameter);
  TEST_THROW(discretization->getGlobalOwnedMap(4), Teuchos::Exceptions::InvalidParameter);

  // basic checks on the 1d map
  Teuchos::RCP<const Epetra_BlockMap> map = discretization->getGlobalOwnedMap(1);
  TEST_ASSERT(map->NumGlobalElements() == 8);
  TEST_ASSERT(map->NumMyElements() == 4);
  TEST_ASSERT(map->ElementSize() == 1);
  TEST_ASSERT(map->IndexBase() == 0);
  TEST_ASSERT(map->UniqueGIDs() == true);
  int* myGlobalElements = map->MyGlobalElements();
  if(rank == 0){
    TEST_ASSERT(myGlobalElements[0] == 0);
    TEST_ASSERT(myGlobalElements[1] == 2);
    TEST_ASSERT(myGlobalElements[2] == 4);
    TEST_ASSERT(myGlobalElements[3] == 6);
  }
  if(rank == 1){
    TEST_ASSERT(myGlobalElements[0] == 5);
    TEST_ASSERT(myGlobalElements[1] == 7);
    TEST_ASSERT(myGlobalElements[2] == 1);
    TEST_ASSERT(myGlobalElements[3] == 3);
  }

  // check the 1d overlap map
  // for this simple discretization, everything should be ghosted on both processors
  Teuchos::RCP<const Epetra_BlockMap> overlapMap = discretization->getGlobalOverlapMap(1);
  TEST_ASSERT(overlapMap->NumGlobalElements() == 16);
  TEST_ASSERT(overlapMap->NumMyElements() == 8);
  TEST_ASSERT(overlapMap->ElementSize() == 1);
  TEST_ASSERT(overlapMap->IndexBase() == 0);
  TEST_ASSERT(overlapMap->UniqueGIDs() == false);
  myGlobalElements = overlapMap->MyGlobalElements();
  if(rank == 0){
    TEST_ASSERT(myGlobalElements[0] == 0);
    TEST_ASSERT(myGlobalElements[1] == 2);
    TEST_ASSERT(myGlobalElements[2] == 4);
    TEST_ASSERT(myGlobalElements[3] == 6);
    TEST_ASSERT(myGlobalElements[4] == 1);
    TEST_ASSERT(myGlobalElements[5] == 3);
    TEST_ASSERT(myGlobalElements[6] == 5);
    TEST_ASSERT(myGlobalElements[7] == 7);
  }
  if(rank == 1){
    TEST_ASSERT(myGlobalElements[0] == 5);
    TEST_ASSERT(myGlobalElements[1] == 7);
    TEST_ASSERT(myGlobalElements[2] == 1);
    TEST_ASSERT(myGlobalElements[3] == 3);
    TEST_ASSERT(myGlobalElements[4] == 0);
    TEST_ASSERT(myGlobalElements[5] == 2);
    TEST_ASSERT(myGlobalElements[6] == 4);
    TEST_ASSERT(myGlobalElements[7] == 6);
  }

  // same checks for 3d map
  map = discretization->getGlobalOwnedMap(3);
  TEST_ASSERT(map->NumGlobalElements() == 8);
  TEST_ASSERT(map->NumMyElements() == 4);
  TEST_ASSERT(map->ElementSize() == 3);
  TEST_ASSERT(map->IndexBase() == 0);
  TEST_ASSERT(map->UniqueGIDs() == true);
  myGlobalElements = map->MyGlobalElements();
  if(rank == 0){
    TEST_ASSERT(myGlobalElements[0] == 0);
    TEST_ASSERT(myGlobalElements[1] == 2);
    TEST_ASSERT(myGlobalElements[2] == 4);
    TEST_ASSERT(myGlobalElements[3] == 6);
  }
  if(rank == 1){
    TEST_ASSERT(myGlobalElements[0] == 5);
    TEST_ASSERT(myGlobalElements[1] == 7);
    TEST_ASSERT(myGlobalElements[2] == 1);
    TEST_ASSERT(myGlobalElements[3] == 3);
  }

  // check the 3d overlap map
  // for this simple discretization, everything should be ghosted on both processors
  overlapMap = discretization->getGlobalOverlapMap(3);
  TEST_ASSERT(overlapMap->NumGlobalElements() == 16);
  TEST_ASSERT(overlapMap->NumMyElements() == 8);
  TEST_ASSERT(overlapMap->ElementSize() == 3);
  TEST_ASSERT(overlapMap->IndexBase() == 0);
  TEST_ASSERT(overlapMap->UniqueGIDs() == false);
  myGlobalElements = overlapMap->MyGlobalElements();
  if(rank == 0){
    TEST_ASSERT(myGlobalElements[0] == 0);
    TEST_ASSERT(myGlobalElements[1] == 2);
    TEST_ASSERT(myGlobalElements[2] == 4);
    TEST_ASSERT(myGlobalElements[3] == 6);
    TEST_ASSERT(myGlobalElements[4] == 1);
    TEST_ASSERT(myGlobalElements[5] == 3);
    TEST_ASSERT(myGlobalElements[6] == 5);
    TEST_ASSERT(myGlobalElements[7] == 7);
  }
  if(rank == 1){
    TEST_ASSERT(myGlobalElements[0] == 5);
    TEST_ASSERT(myGlobalElements[1] == 7);
    TEST_ASSERT(myGlobalElements[2] == 1);
    TEST_ASSERT(myGlobalElements[3] == 3);
    TEST_ASSERT(myGlobalElements[4] == 0);
    TEST_ASSERT(myGlobalElements[5] == 2);
    TEST_ASSERT(myGlobalElements[6] == 4);
    TEST_ASSERT(myGlobalElements[7] == 6);
  }

  // check the bond map
  // the horizon was chosen such that each point should have three neighbors
  // note that if the NeighborhoodType parameter is not set to Spherical, this will fail
  Teuchos::RCP<const Epetra_BlockMap> bondMap = discretization->getGlobalBondMap();
  TEST_ASSERT(bondMap->NumGlobalElements() == 8);
  TEST_ASSERT(bondMap->NumMyElements() == 4);
  TEST_ASSERT(bondMap->IndexBase() == 0);
  TEST_ASSERT(bondMap->UniqueGIDs() == true);
  myGlobalElements = bondMap->MyGlobalElements();
  if(rank == 0){
    TEST_ASSERT(myGlobalElements[0] == 0);
    TEST_ASSERT(myGlobalElements[1] == 2);
    TEST_ASSERT(myGlobalElements[2] == 4);
    TEST_ASSERT(myGlobalElements[3] == 6);
  }
  if(rank == 1){
    TEST_ASSERT(myGlobalElements[0] == 5);
    TEST_ASSERT(myGlobalElements[1] == 7);
    TEST_ASSERT(myGlobalElements[2] == 1);
    TEST_ASSERT(myGlobalElements[3] == 3);
  }
  TEST_ASSERT(discretization->getNumBonds() == 4*3);

  // check the initial positions
  // all three coordinates are contained in a single vector
  Teuchos::RCP<Epetra_Vector> initialX = discretization->getInitialX();
  TEST_ASSERT(initialX->MyLength() == 4*3);
  TEST_ASSERT(initialX->GlobalLength() == 8*3);
  if(rank == 0){
    TEST_FLOATING_EQUALITY((*initialX)[0],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[1],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[2],  0.25, 1.0e-16);
 
    TEST_FLOATING_EQUALITY((*initialX)[3],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[4],  0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[5],  0.25, 1.0e-16);

    TEST_FLOATING_EQUALITY((*initialX)[6],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[7],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[8],  0.75, 1.0e-16);

    TEST_FLOATING_EQUALITY((*initialX)[9],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[10], 0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[11], 0.75, 1.0e-16);
  }
  if(rank == 1){
    TEST_FLOATING_EQUALITY((*initialX)[0],  0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[1],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[2],  0.75, 1.0e-16);

    TEST_FLOATING_EQUALITY((*initialX)[3],  0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[4],  0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[5],  0.75, 1.0e-16);

    TEST_FLOATING_EQUALITY((*initialX)[6],  0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[7],  0.25, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[8],  0.25, 1.0e-16);

    TEST_FLOATING_EQUALITY((*initialX)[9],  0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[10], 0.75, 1.0e-16);
    TEST_FLOATING_EQUALITY((*initialX)[11], 0.25, 1.0e-16);
  }

  // check cell volumes
  Teuchos::RCP<Epetra_Vector> volume = discretization->getCellVolume();
  TEST_ASSERT(volume->MyLength() == 4);
  TEST_ASSERT(volume->GlobalLength() == 8);
  for(int i=0 ; i<volume->MyLength() ; ++i)
    TEST_FLOATING_EQUALITY((*volume)[i], 0.125, 1.0e-16);

  // check the neighbor lists
  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = discretization->getNeighborhoodData();
  TEST_ASSERT(neighborhoodData->NumOwnedPoints() == 4);
  int* ownedIds = neighborhoodData->OwnedIDs();
  TEST_ASSERT(ownedIds[0] == 0);
  TEST_ASSERT(ownedIds[1] == 1);
  TEST_ASSERT(ownedIds[2] == 2);
  TEST_ASSERT(ownedIds[3] == 3);
  TEST_ASSERT(neighborhoodData->NeighborhoodListSize() == 16);
  int* neighborhood = neighborhoodData->NeighborhoodList();
  int* neighborhoodPtr = neighborhoodData->NeighborhoodPtr();
  // remember, these are local IDs on each processor, 
  // which includes both owned and ghost nodes (confusing!)
  if(rank == 0){
    TEST_ASSERT(neighborhoodPtr[0] == 0);
    TEST_ASSERT(neighborhood[0]    == 3);
    TEST_ASSERT(neighborhood[1]    == 4);
    TEST_ASSERT(neighborhood[2]    == 1);
    TEST_ASSERT(neighborhood[3]    == 2);

    TEST_ASSERT(neighborhoodPtr[1] == 4);
    TEST_ASSERT(neighborhood[4]    == 3);
    TEST_ASSERT(neighborhood[5]    == 0);
    TEST_ASSERT(neighborhood[6]    == 5);
    TEST_ASSERT(neighborhood[7]    == 3);

    TEST_ASSERT(neighborhoodPtr[2] == 8);
    TEST_ASSERT(neighborhood[8]    == 3);
    TEST_ASSERT(neighborhood[9]    == 0);
    TEST_ASSERT(neighborhood[10]   == 6);
    TEST_ASSERT(neighborhood[11]   == 3);

    TEST_ASSERT(neighborhoodPtr[3] == 12);
    TEST_ASSERT(neighborhood[12]   == 3);
    TEST_ASSERT(neighborhood[13]   == 1);
    TEST_ASSERT(neighborhood[14]   == 2);
    TEST_ASSERT(neighborhood[15]   == 7);
  }
  if(rank == 1){
    TEST_ASSERT(neighborhoodPtr[0] == 0);
    TEST_ASSERT(neighborhood[0]    == 3);
    TEST_ASSERT(neighborhood[1]    == 2);
    TEST_ASSERT(neighborhood[2]    == 6);
    TEST_ASSERT(neighborhood[3]    == 1);

    TEST_ASSERT(neighborhoodPtr[1] == 4);
    TEST_ASSERT(neighborhood[4]    == 3);
    TEST_ASSERT(neighborhood[5]    == 3);
    TEST_ASSERT(neighborhood[6]    == 0);
    TEST_ASSERT(neighborhood[7]    == 7);

    TEST_ASSERT(neighborhoodPtr[2] == 8);
    TEST_ASSERT(neighborhood[8]    == 3);
    TEST_ASSERT(neighborhood[9]    == 4);
    TEST_ASSERT(neighborhood[10]   == 3);
    TEST_ASSERT(neighborhood[11]   == 0);

    TEST_ASSERT(neighborhoodPtr[3] == 12);
    TEST_ASSERT(neighborhood[12]   == 3);
    TEST_ASSERT(neighborhood[13]   == 2);
    TEST_ASSERT(neighborhood[14]   == 5);
    TEST_ASSERT(neighborhood[15]   == 1);
  }
}



int main
(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
