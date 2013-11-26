/*! \file utPeridigm_PdQuickGridDiscretization.cpp */

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


#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif
#include "Peridigm_PdQuickGridDiscretization.hpp"
#include "Peridigm_HorizonManager.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif

using namespace Teuchos;
using namespace PeridigmNS;

TEUCHOS_UNIT_TEST(PdQuickGridDiscretization, SimpleTensorProductMeshTest) {


  Teuchos::RCP<const Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif
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
  TEST_ASSERT(map->NumMyElements() == 8);
  TEST_ASSERT(map->ElementSize() == 1);
  TEST_ASSERT(map->IndexBase() == 0);
  TEST_ASSERT(map->UniqueGIDs() == true);
  int* myGlobalElements = map->MyGlobalElements();
  for(int i=0 ; i<map->NumMyElements() ; ++i)
    TEST_ASSERT(myGlobalElements[i] == i);

  // for the serial case, the map and the overlap map should match
  Teuchos::RCP<const Epetra_BlockMap> overlapMap = discretization->getGlobalOverlapMap(1);
  TEST_ASSERT(map->SameAs(*overlapMap) == true);

  // same checks for 3d map
  map = discretization->getGlobalOwnedMap(3);
  TEST_ASSERT(map->NumGlobalElements() == 8);
  TEST_ASSERT(map->NumMyElements() == 8);
  TEST_ASSERT(map->ElementSize() == 3);
  TEST_ASSERT(map->IndexBase() == 0);
  TEST_ASSERT(map->UniqueGIDs() == true);
  myGlobalElements = map->MyGlobalElements();
  for(int i=0 ; i<map->NumMyElements() ; ++i)
    TEST_ASSERT(myGlobalElements[i] == i);

  // for the serial case, the map and the overlap map should match
  overlapMap = discretization->getGlobalOverlapMap(3);
  TEST_ASSERT(map->SameAs(*overlapMap) == true);

  // check the bond map
  // the horizon was chosen such that each point should have three neighbors
  // note that if the NeighborhoodType parameter is not set to Spherical, this will fail
  Teuchos::RCP<const Epetra_BlockMap> bondMap = discretization->getGlobalBondMap();
  TEST_ASSERT(bondMap->NumGlobalElements() == 8);
  TEST_ASSERT(bondMap->NumMyElements() == 8);
  TEST_ASSERT(bondMap->ConstantElementSize() == 0);
  TEST_ASSERT(bondMap->IndexBase() == 0);
  TEST_ASSERT(bondMap->UniqueGIDs() == true);
  myGlobalElements = bondMap->MyGlobalElements();
  for(int i=0 ; i<bondMap->NumMyElements() ; ++i)
    TEST_ASSERT(myGlobalElements[i] == i); 

  TEST_ASSERT(discretization->getNumBonds() == 8*3);

  // check the initial positions
  // all three coordinates are contained in a single vector
  Teuchos::RCP<Epetra_Vector> initialX = discretization->getInitialX();
  TEST_ASSERT(initialX->MyLength() == 8*3);
  TEST_ASSERT(initialX->GlobalLength() == 8*3);

  TEST_FLOATING_EQUALITY((*initialX)[0],  0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[1],  0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[2],  0.25, 1.0e-16);

  TEST_FLOATING_EQUALITY((*initialX)[3],  0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[4],  0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[5],  0.25, 1.0e-16);
  
  TEST_FLOATING_EQUALITY((*initialX)[6],  0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[7],  0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[8],  0.25, 1.0e-16);

  TEST_FLOATING_EQUALITY((*initialX)[9],  0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[10], 0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[11], 0.25, 1.0e-16);

  TEST_FLOATING_EQUALITY((*initialX)[12], 0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[13], 0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[14], 0.75, 1.0e-16);

  TEST_FLOATING_EQUALITY((*initialX)[15], 0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[16], 0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[17], 0.75, 1.0e-16);

  TEST_FLOATING_EQUALITY((*initialX)[18], 0.25, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[19], 0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[20], 0.75, 1.0e-16);

  TEST_FLOATING_EQUALITY((*initialX)[21], 0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[22], 0.75, 1.0e-16);
  TEST_FLOATING_EQUALITY((*initialX)[23], 0.75, 1.0e-16);

  // check cell volumes
  Teuchos::RCP<Epetra_Vector> volume = discretization->getCellVolume();
  TEST_ASSERT(volume->MyLength() == 8);
  TEST_ASSERT(volume->GlobalLength() == 8);
  for(int i=0 ; i<volume->MyLength() ; ++i)
    TEST_FLOATING_EQUALITY((*volume)[i], 0.125, 1.0e-16);

  // check the neighbor lists
  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = discretization->getNeighborhoodData();
  TEST_ASSERT(neighborhoodData->NumOwnedPoints() == 8);
  int* ownedIds = neighborhoodData->OwnedIDs();
  for(int i=0 ; i<neighborhoodData->NumOwnedPoints() ; ++i)
    TEST_ASSERT(ownedIds[i] == i);
  TEST_ASSERT(neighborhoodData->NeighborhoodListSize() == 32);
  int* neighborhood = neighborhoodData->NeighborhoodList();
  int* neighborhoodPtr = neighborhoodData->NeighborhoodPtr();

  TEST_ASSERT(neighborhoodPtr[0] == 0);
  TEST_ASSERT(neighborhood[0]    == 3);
  TEST_ASSERT(neighborhood[1]    == 1);
  TEST_ASSERT(neighborhood[2]    == 2);
  TEST_ASSERT(neighborhood[3]    == 4);

  TEST_ASSERT(neighborhoodPtr[1] == 4);
  TEST_ASSERT(neighborhood[4]    == 3);
  TEST_ASSERT(neighborhood[5]    == 0);
  TEST_ASSERT(neighborhood[6]    == 3);
  TEST_ASSERT(neighborhood[7]    == 5);

  TEST_ASSERT(neighborhoodPtr[2] == 8);
  TEST_ASSERT(neighborhood[8]    == 3);
  TEST_ASSERT(neighborhood[9]    == 0);
  TEST_ASSERT(neighborhood[10]   == 3);
  TEST_ASSERT(neighborhood[11]   == 6);

  TEST_ASSERT(neighborhoodPtr[3] == 12);
  TEST_ASSERT(neighborhood[12]   == 3);
  TEST_ASSERT(neighborhood[13]   == 1);
  TEST_ASSERT(neighborhood[14]   == 2);
  TEST_ASSERT(neighborhood[15]   == 7);

  TEST_ASSERT(neighborhoodPtr[4] == 16);
  TEST_ASSERT(neighborhood[16]   == 3);
  TEST_ASSERT(neighborhood[17]   == 0);
  TEST_ASSERT(neighborhood[18]   == 5);
  TEST_ASSERT(neighborhood[19]   == 6);

  TEST_ASSERT(neighborhoodPtr[5] == 20);
  TEST_ASSERT(neighborhood[20]   == 3);
  TEST_ASSERT(neighborhood[21]   == 1);
  TEST_ASSERT(neighborhood[22]   == 4);
  TEST_ASSERT(neighborhood[23]   == 7);

  TEST_ASSERT(neighborhoodPtr[6] == 24);
  TEST_ASSERT(neighborhood[24]   == 3);
  TEST_ASSERT(neighborhood[25]   == 2);
  TEST_ASSERT(neighborhood[26]   == 4);
  TEST_ASSERT(neighborhood[27]   == 7);

  TEST_ASSERT(neighborhoodPtr[7] == 28);
  TEST_ASSERT(neighborhood[28]   == 3);
  TEST_ASSERT(neighborhood[29]   == 3);
  TEST_ASSERT(neighborhood[30]   == 5);
  TEST_ASSERT(neighborhood[31]   == 6);
}


int main
(int argc, char* argv[])
{
 
  // Initialize UTF
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

 
}
