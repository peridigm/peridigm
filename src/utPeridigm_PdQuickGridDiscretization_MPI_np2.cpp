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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#include <Epetra_MpiComm.h>
#include "Peridigm_PdQuickGridDiscretization.hpp"

using namespace boost::unit_test;
using namespace Teuchos;
using namespace PeridigmNS;

void simpleTensorProductMesh()
{
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  Teuchos::RCP<Epetra_Comm> comm;
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  RCP<ParameterList> discParams = rcp(new ParameterList);

  // create a 2x2x2 discretization
  // specify a spherical neighbor search with the horizon a tad longer than the mesh spacing
  discParams->set("Type", "PdQuickGrid");
  discParams->set("NeighborhoodType", "Spherical");
  discParams->set("Horizon", 0.501);
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

  // create the discretization
  RCP<PdQuickGridDiscretization> discretization =
    rcp(new PdQuickGridDiscretization(comm, discParams));

  // sanity check, calling with a dimension other than 1 or 3 should throw an exception
  BOOST_CHECK_THROW(discretization->getGlobalOwnedMap(0), Teuchos::Exceptions::InvalidParameter);
  BOOST_CHECK_THROW(discretization->getGlobalOwnedMap(2), Teuchos::Exceptions::InvalidParameter);
  BOOST_CHECK_THROW(discretization->getGlobalOwnedMap(4), Teuchos::Exceptions::InvalidParameter);

  // basic checks on the 1d map
  Teuchos::RCP<const Epetra_BlockMap> map = discretization->getGlobalOwnedMap(1);
  BOOST_CHECK(map->NumGlobalElements() == 8);
  BOOST_CHECK(map->NumMyElements() == 4);
  BOOST_CHECK(map->ElementSize() == 1);
  BOOST_CHECK(map->IndexBase() == 0);
  BOOST_CHECK(map->UniqueGIDs() == true);
  int* myGlobalElements = map->MyGlobalElements();
  if(rank == 0){
    BOOST_CHECK(myGlobalElements[0] == 0);
    BOOST_CHECK(myGlobalElements[1] == 2);
    BOOST_CHECK(myGlobalElements[2] == 4);
    BOOST_CHECK(myGlobalElements[3] == 6);
  }
  if(rank == 1){
    BOOST_CHECK(myGlobalElements[0] == 5);
    BOOST_CHECK(myGlobalElements[1] == 7);
    BOOST_CHECK(myGlobalElements[2] == 1);
    BOOST_CHECK(myGlobalElements[3] == 3);
  }

  // check the 1d overlap map
  // for this simple discretization, everything should be ghosted on both processors
  Teuchos::RCP<const Epetra_BlockMap> overlapMap = discretization->getGlobalOverlapMap(1);
  BOOST_CHECK(overlapMap->NumGlobalElements() == 16);
  BOOST_CHECK(overlapMap->NumMyElements() == 8);
  BOOST_CHECK(overlapMap->ElementSize() == 1);
  BOOST_CHECK(overlapMap->IndexBase() == 0);
  BOOST_CHECK(overlapMap->UniqueGIDs() == false);
  myGlobalElements = overlapMap->MyGlobalElements();
  if(rank == 0){
    BOOST_CHECK(myGlobalElements[0] == 0);
    BOOST_CHECK(myGlobalElements[1] == 2);
    BOOST_CHECK(myGlobalElements[2] == 4);
    BOOST_CHECK(myGlobalElements[3] == 6);
    BOOST_CHECK(myGlobalElements[4] == 1);
    BOOST_CHECK(myGlobalElements[5] == 3);
    BOOST_CHECK(myGlobalElements[6] == 5);
    BOOST_CHECK(myGlobalElements[7] == 7);
  }
  if(rank == 1){
    BOOST_CHECK(myGlobalElements[0] == 5);
    BOOST_CHECK(myGlobalElements[1] == 7);
    BOOST_CHECK(myGlobalElements[2] == 1);
    BOOST_CHECK(myGlobalElements[3] == 3);
    BOOST_CHECK(myGlobalElements[4] == 0);
    BOOST_CHECK(myGlobalElements[5] == 2);
    BOOST_CHECK(myGlobalElements[6] == 4);
    BOOST_CHECK(myGlobalElements[7] == 6);
  }

  // same checks for 3d map
  map = discretization->getGlobalOwnedMap(3);
  BOOST_CHECK(map->NumGlobalElements() == 8);
  BOOST_CHECK(map->NumMyElements() == 4);
  BOOST_CHECK(map->ElementSize() == 3);
  BOOST_CHECK(map->IndexBase() == 0);
  BOOST_CHECK(map->UniqueGIDs() == true);
  myGlobalElements = map->MyGlobalElements();
  if(rank == 0){
    BOOST_CHECK(myGlobalElements[0] == 0);
    BOOST_CHECK(myGlobalElements[1] == 2);
    BOOST_CHECK(myGlobalElements[2] == 4);
    BOOST_CHECK(myGlobalElements[3] == 6);
  }
  if(rank == 1){
    BOOST_CHECK(myGlobalElements[0] == 5);
    BOOST_CHECK(myGlobalElements[1] == 7);
    BOOST_CHECK(myGlobalElements[2] == 1);
    BOOST_CHECK(myGlobalElements[3] == 3);
  }

  // check the 3d overlap map
  // for this simple discretization, everything should be ghosted on both processors
  overlapMap = discretization->getGlobalOverlapMap(3);
  BOOST_CHECK(overlapMap->NumGlobalElements() == 16);
  BOOST_CHECK(overlapMap->NumMyElements() == 8);
  BOOST_CHECK(overlapMap->ElementSize() == 3);
  BOOST_CHECK(overlapMap->IndexBase() == 0);
  BOOST_CHECK(overlapMap->UniqueGIDs() == false);
  myGlobalElements = overlapMap->MyGlobalElements();
  if(rank == 0){
    BOOST_CHECK(myGlobalElements[0] == 0);
    BOOST_CHECK(myGlobalElements[1] == 2);
    BOOST_CHECK(myGlobalElements[2] == 4);
    BOOST_CHECK(myGlobalElements[3] == 6);
    BOOST_CHECK(myGlobalElements[4] == 1);
    BOOST_CHECK(myGlobalElements[5] == 3);
    BOOST_CHECK(myGlobalElements[6] == 5);
    BOOST_CHECK(myGlobalElements[7] == 7);
  }
  if(rank == 1){
    BOOST_CHECK(myGlobalElements[0] == 5);
    BOOST_CHECK(myGlobalElements[1] == 7);
    BOOST_CHECK(myGlobalElements[2] == 1);
    BOOST_CHECK(myGlobalElements[3] == 3);
    BOOST_CHECK(myGlobalElements[4] == 0);
    BOOST_CHECK(myGlobalElements[5] == 2);
    BOOST_CHECK(myGlobalElements[6] == 4);
    BOOST_CHECK(myGlobalElements[7] == 6);
  }

  // check the bond map
  // the horizon was chosen such that each point should have three neighbors
  // note that if the NeighborhoodType parameter is not set to Spherical, this will fail
  Teuchos::RCP<const Epetra_BlockMap> bondMap = discretization->getGlobalBondMap();
  BOOST_CHECK(bondMap->NumGlobalElements() == 8);
  BOOST_CHECK(bondMap->NumMyElements() == 4);
  BOOST_CHECK(bondMap->IndexBase() == 0);
  BOOST_CHECK(bondMap->UniqueGIDs() == true);
  myGlobalElements = bondMap->MyGlobalElements();
  if(rank == 0){
    BOOST_CHECK(myGlobalElements[0] == 0);
    BOOST_CHECK(myGlobalElements[1] == 2);
    BOOST_CHECK(myGlobalElements[2] == 4);
    BOOST_CHECK(myGlobalElements[3] == 6);
  }
  if(rank == 1){
    BOOST_CHECK(myGlobalElements[0] == 5);
    BOOST_CHECK(myGlobalElements[1] == 7);
    BOOST_CHECK(myGlobalElements[2] == 1);
    BOOST_CHECK(myGlobalElements[3] == 3);
  }
  BOOST_CHECK(discretization->getNumBonds() == 4*3);

  // check the initial positions
  // all three coordinates are contained in a single vector
  Teuchos::RCP<Epetra_Vector> initialX = discretization->getInitialX();
  BOOST_CHECK(initialX->MyLength() == 4*3);
  BOOST_CHECK(initialX->GlobalLength() == 8*3);
  if(rank == 0){
    BOOST_CHECK_CLOSE((*initialX)[0],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[1],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[2],  0.25, 1.0e-16);
 
    BOOST_CHECK_CLOSE((*initialX)[3],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[4],  0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[5],  0.25, 1.0e-16);

    BOOST_CHECK_CLOSE((*initialX)[6],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[7],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[8],  0.75, 1.0e-16);

    BOOST_CHECK_CLOSE((*initialX)[9],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[10], 0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[11], 0.75, 1.0e-16);
  }
  if(rank == 1){
    BOOST_CHECK_CLOSE((*initialX)[0],  0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[1],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[2],  0.75, 1.0e-16);

    BOOST_CHECK_CLOSE((*initialX)[3],  0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[4],  0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[5],  0.75, 1.0e-16);

    BOOST_CHECK_CLOSE((*initialX)[6],  0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[7],  0.25, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[8],  0.25, 1.0e-16);

    BOOST_CHECK_CLOSE((*initialX)[9],  0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[10], 0.75, 1.0e-16);
    BOOST_CHECK_CLOSE((*initialX)[11], 0.25, 1.0e-16);
  }

  // check cell volumes
  Teuchos::RCP<Epetra_Vector> volume = discretization->getCellVolume();
  BOOST_CHECK(volume->MyLength() == 4);
  BOOST_CHECK(volume->GlobalLength() == 8);
  for(int i=0 ; i<volume->MyLength() ; ++i)
    BOOST_CHECK_CLOSE((*volume)[i], 0.125, 1.0e-16);

  // check the neighbor lists
  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = discretization->getNeighborhoodData();
  BOOST_CHECK(neighborhoodData->NumOwnedPoints() == 4);
  int* ownedIds = neighborhoodData->OwnedIDs();
  BOOST_CHECK(ownedIds[0] == 0);
  BOOST_CHECK(ownedIds[1] == 1);
  BOOST_CHECK(ownedIds[2] == 2);
  BOOST_CHECK(ownedIds[3] == 3);
  BOOST_CHECK(neighborhoodData->NeighborhoodListSize() == 16);
  int* neighborhood = neighborhoodData->NeighborhoodList();
  int* neighborhoodPtr = neighborhoodData->NeighborhoodPtr();
  // remember, these are local IDs on each processor, 
  // which includes both owned and ghost nodes (confusing!)
  if(rank == 0){
    BOOST_CHECK(neighborhoodPtr[0] == 0);
    BOOST_CHECK(neighborhood[0]    == 3);
    BOOST_CHECK(neighborhood[1]    == 4);
    BOOST_CHECK(neighborhood[2]    == 1);
    BOOST_CHECK(neighborhood[3]    == 2);

    BOOST_CHECK(neighborhoodPtr[1] == 4);
    BOOST_CHECK(neighborhood[4]    == 3);
    BOOST_CHECK(neighborhood[5]    == 0);
    BOOST_CHECK(neighborhood[6]    == 5);
    BOOST_CHECK(neighborhood[7]    == 3);

    BOOST_CHECK(neighborhoodPtr[2] == 8);
    BOOST_CHECK(neighborhood[8]    == 3);
    BOOST_CHECK(neighborhood[9]    == 0);
    BOOST_CHECK(neighborhood[10]   == 6);
    BOOST_CHECK(neighborhood[11]   == 3);

    BOOST_CHECK(neighborhoodPtr[3] == 12);
    BOOST_CHECK(neighborhood[12]   == 3);
    BOOST_CHECK(neighborhood[13]   == 1);
    BOOST_CHECK(neighborhood[14]   == 2);
    BOOST_CHECK(neighborhood[15]   == 7);
  }
  if(rank == 1){
    BOOST_CHECK(neighborhoodPtr[0] == 0);
    BOOST_CHECK(neighborhood[0]    == 3);
    BOOST_CHECK(neighborhood[1]    == 2);
    BOOST_CHECK(neighborhood[2]    == 6);
    BOOST_CHECK(neighborhood[3]    == 1);

    BOOST_CHECK(neighborhoodPtr[1] == 4);
    BOOST_CHECK(neighborhood[4]    == 3);
    BOOST_CHECK(neighborhood[5]    == 3);
    BOOST_CHECK(neighborhood[6]    == 0);
    BOOST_CHECK(neighborhood[7]    == 7);

    BOOST_CHECK(neighborhoodPtr[2] == 8);
    BOOST_CHECK(neighborhood[8]    == 3);
    BOOST_CHECK(neighborhood[9]    == 4);
    BOOST_CHECK(neighborhood[10]   == 3);
    BOOST_CHECK(neighborhood[11]   == 0);

    BOOST_CHECK(neighborhoodPtr[3] == 12);
    BOOST_CHECK(neighborhood[12]   == 3);
    BOOST_CHECK(neighborhood[13]   == 2);
    BOOST_CHECK(neighborhood[14]   == 5);
    BOOST_CHECK(neighborhood[15]   == 1);
  }
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success = true;

	test_suite* proc = BOOST_TEST_SUITE("utPeridigm_PdQuickGridDiscretization_MPI_np2");
	proc->add(BOOST_TEST_CASE(&simpleTensorProductMesh));
	framework::master_test_suite().add(proc);

	return success;
}

bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int return_code = -1;

  if(numProcs == 2){
	return_code = unit_test_main(init_unit_test, argc, argv);
  }
  else{
	std::cerr << "Unit test runtime ERROR: utPeridigm_PdQuickGridDiscretization_MPI_np2 only makes sense on 2 processors" << std::endl;
	std::cerr << "Re-run unit test $mpiexec -np 2 ./utPeridigm_PdQuickGridDiscretization_MPI_np2" << std::endl;
  }

  MPI_Finalize();

  return return_code;
}
