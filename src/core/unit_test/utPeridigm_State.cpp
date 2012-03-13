/*! \file utPeridigm_State.cpp */

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
#include <boost/test/parameterized_test.hpp>
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>
#include "Peridigm_State.hpp"
#include <vector>

using namespace boost::unit_test;
using namespace Teuchos;
using namespace PeridigmNS;

//! Create a two-point problem for testing.
PeridigmNS::State createTwoPointProblem()
{
  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

  // set up a hard-coded layout for two points
  int numCells = 2;

  // set up overlap maps, which include ghosted nodes
  int numGlobalElements(numCells), numMyElements(2), elementSize(1), indexBase(0);
  std::vector<int> myGlobalElements(numMyElements);
  for(int i=0; i<numMyElements ; ++i)
    myGlobalElements[i] = i;

  // overlapScalarPointMap
  // used for cell volumes and scalar constitutive data
  Teuchos::RCP<Epetra_BlockMap> overlapScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));
  // overlapVectorPointMap
  // used for positions, displacements, velocities and vector constitutive data
  elementSize = 3;
  Teuchos::RCP<Epetra_BlockMap> overlapVectorPointMap = 
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm)); 
  // ownedScalarBondMap
  // used for bond damage and bond constitutive data
  std::vector<int> bondElementSize(numMyElements, 1);
  Teuchos::RCP<Epetra_BlockMap> ownedScalarBondMap = 
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], &bondElementSize[0], indexBase, *comm));

  // create a state object
  State state;

  // create a list of scalar field specs and allocate the data
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > scalarPointFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  scalarPointFieldSpecs->push_back(Field_NS::VOLUME);
  scalarPointFieldSpecs->push_back(Field_NS::DAMAGE);
  scalarPointFieldSpecs->push_back(Field_NS::GID);
  scalarPointFieldSpecs->push_back(Field_NS::PROC_NUM);
  scalarPointFieldSpecs->push_back(Field_NS::WEIGHTED_VOLUME);
  scalarPointFieldSpecs->push_back(Field_NS::DILATATION);
  scalarPointFieldSpecs->push_back(Field_NS::NUM_NEIGHBORS);
  scalarPointFieldSpecs->push_back(Field_NS::LAMBDA);
  scalarPointFieldSpecs->push_back(Field_NS::NORM_TD);
  scalarPointFieldSpecs->push_back(Field_NS::DSF);
  scalarPointFieldSpecs->push_back(Field_NS::BC_MASK);
  state.allocateScalarPointData(scalarPointFieldSpecs, overlapScalarPointMap);
  BOOST_CHECK_EQUAL( state.getScalarPointMultiVector()->NumVectors(), (int)scalarPointFieldSpecs->size() );
  BOOST_CHECK_EQUAL( state.getScalarPointMultiVector()->MyLength(), overlapScalarPointMap->NumMyPoints() );
  BOOST_CHECK( state.getScalarPointMultiVector()->Map().SameAs( *overlapScalarPointMap ) );

  // create a list of vector field specs and allocate the data
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > vectorFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  vectorFieldSpecs->push_back(Field_NS::COORD3D);
  vectorFieldSpecs->push_back(Field_NS::DISPL3D);
  vectorFieldSpecs->push_back(Field_NS::CURCOORD3D);
  vectorFieldSpecs->push_back(Field_NS::VELOC3D);
  vectorFieldSpecs->push_back(Field_NS::ACCEL3D);
  vectorFieldSpecs->push_back(Field_NS::FORCE3D);
  vectorFieldSpecs->push_back(Field_NS::FORCE_DENSITY3D);
  vectorFieldSpecs->push_back(Field_NS::CONTACT_FORCE3D);
  vectorFieldSpecs->push_back(Field_NS::CONTACT_FORCE_DENSITY3D);
  vectorFieldSpecs->push_back(Field_NS::RESID3D);
  state.allocateVectorPointData(vectorFieldSpecs, overlapVectorPointMap);
  BOOST_CHECK_EQUAL( state.getVectorPointMultiVector()->NumVectors(), (int)vectorFieldSpecs->size() );
  BOOST_CHECK_EQUAL( state.getVectorPointMultiVector()->MyLength(), overlapVectorPointMap->NumMyPoints() );
  BOOST_CHECK( state.getVectorPointMultiVector()->Map().SameAs( *overlapVectorPointMap ) );

  // create a list of bond field specs and allocate the data
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > bondFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  bondFieldSpecs->push_back(Field_NS::BOND_DAMAGE);
  bondFieldSpecs->push_back(Field_NS::DEVIATORIC_PLASTIC_EXTENSION);
  bondFieldSpecs->push_back(Field_NS::DEVIATORIC_BACK_EXTENSION);
  state.allocateScalarBondData(bondFieldSpecs, ownedScalarBondMap);
  BOOST_CHECK_EQUAL( state.getScalarBondMultiVector()->NumVectors(), (int)bondFieldSpecs->size() );
  BOOST_CHECK_EQUAL( state.getScalarBondMultiVector()->MyLength(), ownedScalarBondMap->NumMyPoints() );
  BOOST_CHECK( state.getScalarBondMultiVector()->Map().SameAs( *ownedScalarBondMap ) );

  return state;
}

//! Create a three-point problem for testing.
PeridigmNS::State createThreePointProblem()
{
  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif
  int numProcs = comm->NumProc();
  int myPID = comm->MyPID();

  // set up a hard-coded layout for three points
  int numCells = 3;

  // set up overlap maps, which include ghosted nodes
  int numGlobalElements(numCells), numMyElements(0), elementSize(1), indexBase(0);
  std::vector<int> myGlobalElements;
  if(numProcs == 1){
    numMyElements = 3;
    myGlobalElements.resize(numMyElements);
    myGlobalElements[0] = 0;
    myGlobalElements[1] = 1;
    myGlobalElements[2] = 2;
  }
  else if(numProcs == 2){
    if(myPID == 0){
      numMyElements = 2;
      myGlobalElements.resize(numMyElements);
      myGlobalElements[0] = 0;
      myGlobalElements[1] = 1;
    }
    else if(myPID == 1){
      numMyElements = 2;
      myGlobalElements.resize(numMyElements);
      myGlobalElements[0] = 2;
      myGlobalElements[1] = 1;
    }
  }

  // overlapScalarPointMap
  // used for cell volumes and scalar constitutive data
  Teuchos::RCP<Epetra_BlockMap> overlapScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(-1, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));
  // overlapVectorPointMap
  // used for positions, displacements, velocities and vector constitutive data
  elementSize = 3;
  Teuchos::RCP<Epetra_BlockMap> overlapVectorPointMap = 
    Teuchos::rcp(new Epetra_BlockMap(-1, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm)); 
  // ownedScalarBondMap
  // used for bond damage and bond constitutive data
  std::vector<int> bondElementSize;
  if(numProcs == 1){
    numMyElements = 3;
    myGlobalElements.resize(numMyElements);
    myGlobalElements[0] = 0;
    myGlobalElements[1] = 1;
    myGlobalElements[2] = 2;
    bondElementSize.resize(numMyElements);
    bondElementSize[0] = 1;
    bondElementSize[1] = 2;
    bondElementSize[2] = 1;
  }
  else if(numProcs == 2){
    if(myPID == 0){
      numMyElements = 2;
      myGlobalElements.resize(numMyElements);
      myGlobalElements[0] = 0;
      myGlobalElements[1] = 1;
      bondElementSize.resize(numMyElements);
      bondElementSize[0] = 1;
      bondElementSize[1] = 2;
    }
    else if(myPID == 1){
      numMyElements = 1;
      myGlobalElements.resize(numMyElements);
      myGlobalElements[0] = 2;
      bondElementSize.resize(numMyElements);
      bondElementSize[0] = 1;
    }
  }
  Teuchos::RCP<Epetra_BlockMap> ownedScalarBondMap = 
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], &bondElementSize[0], indexBase, *comm));

  PeridigmNS::State state;

  // create a list of scalar field specs and allocate the data
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > scalarPointFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  scalarPointFieldSpecs->push_back(Field_NS::VOLUME);
  scalarPointFieldSpecs->push_back(Field_NS::DAMAGE);
  scalarPointFieldSpecs->push_back(Field_NS::GID);
  scalarPointFieldSpecs->push_back(Field_NS::PROC_NUM);
  scalarPointFieldSpecs->push_back(Field_NS::WEIGHTED_VOLUME);
  scalarPointFieldSpecs->push_back(Field_NS::DILATATION);
  scalarPointFieldSpecs->push_back(Field_NS::NUM_NEIGHBORS);
  scalarPointFieldSpecs->push_back(Field_NS::LAMBDA);
  scalarPointFieldSpecs->push_back(Field_NS::NORM_TD);
  scalarPointFieldSpecs->push_back(Field_NS::DSF);
  scalarPointFieldSpecs->push_back(Field_NS::BC_MASK);
  state.allocateScalarPointData(scalarPointFieldSpecs, overlapScalarPointMap);
  BOOST_CHECK_EQUAL( state.getScalarPointMultiVector()->NumVectors(), (int)scalarPointFieldSpecs->size() );
  BOOST_CHECK_EQUAL( state.getScalarPointMultiVector()->MyLength(), overlapScalarPointMap->NumMyPoints() );
  BOOST_CHECK( state.getScalarPointMultiVector()->Map().SameAs( *overlapScalarPointMap ) );

  // create a list of vector field specs and allocate the data
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > vectorFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  vectorFieldSpecs->push_back(Field_NS::COORD3D);
  vectorFieldSpecs->push_back(Field_NS::DISPL3D);
  vectorFieldSpecs->push_back(Field_NS::CURCOORD3D);
  vectorFieldSpecs->push_back(Field_NS::VELOC3D);
  vectorFieldSpecs->push_back(Field_NS::ACCEL3D);
  vectorFieldSpecs->push_back(Field_NS::FORCE3D);
  vectorFieldSpecs->push_back(Field_NS::FORCE_DENSITY3D);
  vectorFieldSpecs->push_back(Field_NS::CONTACT_FORCE3D);
  vectorFieldSpecs->push_back(Field_NS::CONTACT_FORCE_DENSITY3D);
  vectorFieldSpecs->push_back(Field_NS::RESID3D);
  state.allocateVectorPointData(vectorFieldSpecs, overlapVectorPointMap);
  BOOST_CHECK_EQUAL( state.getVectorPointMultiVector()->NumVectors(), (int)vectorFieldSpecs->size() );
  BOOST_CHECK_EQUAL( state.getVectorPointMultiVector()->MyLength(), overlapVectorPointMap->NumMyPoints() );
  BOOST_CHECK( state.getVectorPointMultiVector()->Map().SameAs( *overlapVectorPointMap ) );

  // create a list of bond field specs and allocate the data
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > bondFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  bondFieldSpecs->push_back(Field_NS::BOND_DAMAGE);
  bondFieldSpecs->push_back(Field_NS::DEVIATORIC_PLASTIC_EXTENSION);
  bondFieldSpecs->push_back(Field_NS::DEVIATORIC_BACK_EXTENSION);
  state.allocateScalarBondData(bondFieldSpecs, ownedScalarBondMap);
  BOOST_CHECK_EQUAL( state.getScalarBondMultiVector()->NumVectors(), (int)bondFieldSpecs->size() );
  BOOST_CHECK_EQUAL( state.getScalarBondMultiVector()->MyLength(), ownedScalarBondMap->NumMyPoints() );
  BOOST_CHECK( state.getScalarBondMultiVector()->Map().SameAs( *ownedScalarBondMap ) );

  return state;
}

//! Create a State object for a two-point problem, check data storage and basic functionality.
void twoPointProblem()
{
  PeridigmNS::State state = createTwoPointProblem();

  // check initialization of data to zero
  Teuchos::RCP<Epetra_Vector> coordinates = state.getData(Field_NS::COORD3D);
  for(int i=0 ; i<coordinates->MyLength() ; ++i)
    BOOST_CHECK_EQUAL( (*coordinates)[i], 0.0 );

  // set some data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(Field_NS::GID));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      ids[i] = i;
    // vector data
    Epetra_Vector& force = *(state.getData(Field_NS::FORCE3D));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      force[i*3] = i*3;
      force[i*3+1] = i*3+1;
      force[i*3+2] = i*3+2;
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(Field_NS::DAMAGE));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        bondDamage[firstPointInElement+j] = firstPointInElement+j;
    }
  }

  // check the data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(Field_NS::GID));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      BOOST_CHECK_CLOSE(ids[i], (double)(i), 1.0e-14);
    // vector data
    Epetra_Vector& force = *(state.getData(Field_NS::FORCE3D));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      BOOST_CHECK_CLOSE(force[i*3], (double)(i*3), 1.0e-14);
      BOOST_CHECK_CLOSE(force[i*3+1], (double)(i*3+1), 1.0e-14);
      BOOST_CHECK_CLOSE(force[i*3+2], (double)(i*3+2), 1.0e-14);
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(Field_NS::DAMAGE));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        BOOST_CHECK_CLOSE(bondDamage[firstPointInElement+j], (double)(firstPointInElement+j), 1.0e-14);
    }
  }
}

//! Create a State object for a three-point problem, check data storage and basic functionality.
void threePointProblem()
{
  PeridigmNS::State state = createThreePointProblem();

  // check initialization of data to zero
  Teuchos::RCP<Epetra_Vector> coordinates = state.getData(Field_NS::COORD3D);
  for(int i=0 ; i<coordinates->MyLength() ; ++i)
    BOOST_CHECK_EQUAL( (*coordinates)[i], 0.0 );

  // set some data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(Field_NS::GID));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      ids[i] = i;
    // vector data
    Epetra_Vector& force = *(state.getData(Field_NS::FORCE3D));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      force[i*3] = i*3;
      force[i*3+1] = i*3+1;
      force[i*3+2] = i*3+2;
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(Field_NS::DAMAGE));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        bondDamage[firstPointInElement+j] = firstPointInElement+j;
    }
  }

  // check the data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(Field_NS::GID));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      BOOST_CHECK_CLOSE(ids[i], (double)(i), 1.0e-14);
    // vector data
    Epetra_Vector& force = *(state.getData(Field_NS::FORCE3D));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      BOOST_CHECK_CLOSE(force[i*3], (double)(i*3), 1.0e-14);
      BOOST_CHECK_CLOSE(force[i*3+1], (double)(i*3+1), 1.0e-14);
      BOOST_CHECK_CLOSE(force[i*3+2], (double)(i*3+2), 1.0e-14);
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(Field_NS::DAMAGE));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        BOOST_CHECK_CLOSE(bondDamage[firstPointInElement+j], (double)(firstPointInElement+j), 1.0e-14);
    }
  }
}

//! Test ability to copy data from one State to another.
void copyFrom()
{
  PeridigmNS::State state = createThreePointProblem();

  // set some data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(Field_NS::GID));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      ids[i] = i;
    // vector data
    Epetra_Vector& force = *(state.getData(Field_NS::FORCE3D));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      force[i*3] = i*3;
      force[i*3+1] = i*3+1;
      force[i*3+2] = i*3+2;
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(Field_NS::DAMAGE));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j){
        bondDamage[firstPointInElement+j] = firstPointInElement+j;
      }
    }
  }

  // create a temporary State
  // this mimics what is done behind the scenes in the finite-difference Jacobian routine
  // there is one owned point, processor 0 has GID 0, and processor 1 (if it exists) has GID 2
  // in either case there is one ghosted point, its neighbor, GID = 1
  int procID = state.getScalarPointMultiVector()->Comm().MyPID();
  Epetra_SerialComm serialComm;
  int numOwnedIDs = 1;
  std::vector<int> tempMyGlobalIDs(2); // includes ghost
  // main owned point
  if(procID == 0)
    tempMyGlobalIDs[0] = 0;
  else if(procID == 1)
    tempMyGlobalIDs[0] = 2;
  // neighbor
  tempMyGlobalIDs[1] = 1;
  std::vector<int> bondElementSize(1);
  bondElementSize[0] = 1;
  Teuchos::RCP<Epetra_BlockMap> tempOverlapScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(tempMyGlobalIDs.size(), tempMyGlobalIDs.size(), &tempMyGlobalIDs[0], 1, 0, serialComm));
  Teuchos::RCP<Epetra_BlockMap> tempOverlapVectorPointMap =
    Teuchos::rcp(new Epetra_BlockMap(tempMyGlobalIDs.size(), tempMyGlobalIDs.size(), &tempMyGlobalIDs[0], 3, 0, serialComm));
  Teuchos::RCP<Epetra_BlockMap> tempOwnedScalarBondMap = 
    Teuchos::rcp(new Epetra_BlockMap(numOwnedIDs, numOwnedIDs, &tempMyGlobalIDs[0], &bondElementSize[0], 0, serialComm));

  PeridigmNS::State tempState;
  
  Teuchos::RCP<Field_ENUM::Relation> relation = Teuchos::rcp(new Field_ENUM::Relation);
  Teuchos::RCP<Field_ENUM::Length> length = Teuchos::rcp(new Field_ENUM::Length);

  // allocate data for the field specs, which are taken from the initial state object

  // scalar point data
  *relation = Field_ENUM::ELEMENT;
  *length = Field_ENUM::SCALAR;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > scalarPointFieldSpecs = state.getFieldSpecs(relation, length);
  tempState.allocateScalarPointData(scalarPointFieldSpecs, tempOverlapScalarPointMap);
  BOOST_CHECK_EQUAL( tempState.getScalarPointMultiVector()->NumVectors(), state.getScalarPointMultiVector()->NumVectors());
  BOOST_CHECK_EQUAL( tempState.getScalarPointMultiVector()->MyLength(), tempOverlapScalarPointMap->NumMyPoints() );
  BOOST_CHECK( tempState.getScalarPointMultiVector()->Map().SameAs( *tempOverlapScalarPointMap ) );

  // vector point data
  *relation = Field_ENUM::NODE;
  *length = Field_ENUM::VECTOR3D;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > vectorFieldSpecs = state.getFieldSpecs(relation, length);
  tempState.allocateVectorPointData(vectorFieldSpecs, tempOverlapVectorPointMap);
  BOOST_CHECK_EQUAL( tempState.getVectorPointMultiVector()->NumVectors(), state.getVectorPointMultiVector()->NumVectors());
  BOOST_CHECK_EQUAL( tempState.getVectorPointMultiVector()->MyLength(), tempOverlapVectorPointMap->NumMyPoints() );
  BOOST_CHECK( tempState.getVectorPointMultiVector()->Map().SameAs( *tempOverlapVectorPointMap ) );

  // scalar bond data
  *relation = Field_ENUM::BOND;
  *length = Field_ENUM::SCALAR;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > bondFieldSpecs = state.getFieldSpecs(relation, length);
  tempState.allocateScalarBondData(bondFieldSpecs, tempOwnedScalarBondMap);
  BOOST_CHECK_EQUAL( tempState.getScalarBondMultiVector()->NumVectors(), state.getScalarBondMultiVector()->NumVectors());
  BOOST_CHECK_EQUAL( tempState.getScalarBondMultiVector()->MyLength(), tempOwnedScalarBondMap->NumMyPoints() );
  BOOST_CHECK( tempState.getScalarBondMultiVector()->Map().SameAs( *tempOwnedScalarBondMap ) );

  // copy the data from the initial State to the smaller, temporary state
  tempState.copyLocallyOwnedDataFromState(Teuchos::RCP<PeridigmNS::State>(&state, false));

  // check the data
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > tempFieldSpecs = tempState.getFieldSpecs();
  for(unsigned int iSpec=0 ; iSpec<tempFieldSpecs->size() ; ++iSpec){
    Teuchos::RCP<Epetra_Vector> data = state.getData( (*tempFieldSpecs)[iSpec] );
    Teuchos::RCP<Epetra_Vector> tempData = tempState.getData( (*tempFieldSpecs)[iSpec] );
    for(int iLID=0 ; iLID<tempData->Map().NumMyElements() ; ++iLID){
      int globalID = tempData->Map().GID(iLID);
      BOOST_CHECK(globalID != -1);
      int stateLID = data->Map().LID(globalID);
      BOOST_CHECK(stateLID != -1);
      int tempElementSize = tempData->Map().ElementSize(iLID);
      int stateElementSize = data->Map().ElementSize(stateLID);
      BOOST_CHECK_EQUAL( tempElementSize, stateElementSize );
      int tempFirstPointInElement = tempData->Map().FirstPointInElement(iLID);
      int stateFirstPointInElement = data->Map().FirstPointInElement(stateLID);
      for(int i=0 ; i<tempElementSize ; ++i){
        BOOST_CHECK_EQUAL( (*tempData)[tempFirstPointInElement + i], (*data)[stateFirstPointInElement + i] );
      }
    }
  }
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success = true;

	test_suite* proc = BOOST_TEST_SUITE("utPeridigm_State");
	proc->add(BOOST_TEST_CASE(&twoPointProblem));
	proc->add(BOOST_TEST_CASE(&threePointProblem));
	proc->add(BOOST_TEST_CASE(&copyFrom));
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
  int numProcs = 1;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

  int returnCode = -1;
  
  if(numProcs == 1 || numProcs == 2){
    returnCode = unit_test_main(init_unit_test, argc, argv);
  }
  else{
    std::cerr << "Unit test runtime ERROR: utPeridigm_State only makes sense on 1 or 2 processors." << std::endl;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return returnCode;
}
