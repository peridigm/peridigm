/*! \file utPeridigm_State.cpp  with Teuchos Unit test Library*/

//@HEADER
// ************************************************************************
//
// ************************************************************************
//@HEADER 

#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#include <Epetra_SerialComm.h>
#include "Peridigm_State.hpp"
#include <vector>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif

using namespace Teuchos;
using namespace PeridigmNS;
using namespace std;

//! Create a two-point problem for testing.
PeridigmNS::State createTwoPointProblem(Teuchos::RCP<Epetra_Comm> comm, Teuchos::RCP<Epetra_BlockMap> &overlapScalarPointMap, Teuchos::RCP<Epetra_BlockMap> &overlapVectorPointMap, Teuchos::RCP<Epetra_BlockMap> &ownedScalarBondMap, vector<int> &scalarPointFieldIds, vector<int> &vectorPointFieldIds, vector<int> &bondFieldIds)
{
  
  // set up a hard-coded layout for two points
  int numCells = 2;

  // set up overlap maps, which include ghosted nodes
  int numGlobalElements(numCells), numMyElements(2), elementSize(1), indexBase(0);
  std::vector<int> myGlobalElements(numMyElements);
  for(int i=0; i<numMyElements ; ++i)
    myGlobalElements[i] = i;

  // overlapScalarPointMap
  // used for cell volumes and scalar constitutive data
  overlapScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));
  // overlapVectorPointMap
  // used for positions, displacements, velocities and vector constitutive data
  elementSize = 3;
  overlapVectorPointMap = 
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm)); 
  // ownedScalarBondMap
  // used for bond damage and bond constitutive data
  std::vector<int> bondElementSize(numMyElements, 1);
  ownedScalarBondMap = 
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], &bondElementSize[0], indexBase, *comm));

  // create a state object
  State state;

  FieldManager& fm = FieldManager::self();

  // create a list of scalar field specs and allocate the data
  
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Element_Id") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Proc_Num") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Number_Of_Neighbors") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Lambda") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Surface_Correction_Factor") );

  state.allocatePointData(PeridigmField::SCALAR, scalarPointFieldIds, overlapScalarPointMap);

  // create a list of vector field specs and allocate the data
  
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Contact_Force_Density") );

  state.allocatePointData(PeridigmField::VECTOR,vectorPointFieldIds, overlapVectorPointMap);

  // create a list of bond field specs and allocate the data
  
  bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage") );
  bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deviatoric_Plastic_Extension") );
  state.allocateBondData(bondFieldIds, ownedScalarBondMap);

  return state;
}

//! Create a State object for a two-point problem, check data storage and basic functionality.

TEUCHOS_UNIT_TEST(State, TwoPointTest) {

  Teuchos::RCP<Epetra_Comm> comm;
  
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

 
  PeridigmNS::State state;
  Teuchos::RCP<Epetra_BlockMap> overlapScalarPointMap;
  Teuchos::RCP<Epetra_BlockMap> overlapVectorPointMap;
  Teuchos::RCP<Epetra_BlockMap> ownedScalarBondMap;
  vector<int> scalarPointFieldIds;
  vector<int> vectorPointFieldIds;
  vector<int> bondFieldIds;
  
  state = createTwoPointProblem(comm, overlapScalarPointMap, overlapVectorPointMap, ownedScalarBondMap, scalarPointFieldIds, vectorPointFieldIds, bondFieldIds);


  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::SCALAR)->NumVectors(), (int)scalarPointFieldIds.size() );
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::SCALAR)->MyLength(), overlapScalarPointMap->NumMyPoints() );
  TEST_ASSERT( state.getPointMultiVector(PeridigmField::SCALAR)->Map().SameAs( *overlapScalarPointMap ) );

  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::VECTOR)->NumVectors(), (int)vectorPointFieldIds.size() );
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::VECTOR)->MyLength(), overlapVectorPointMap->NumMyPoints() );
  TEST_ASSERT( state.getPointMultiVector(PeridigmField::VECTOR)->Map().SameAs( *overlapVectorPointMap ) );

  TEST_EQUALITY( state.getBondMultiVector()->NumVectors(), (int)bondFieldIds.size() );
  TEST_EQUALITY( state.getBondMultiVector()->MyLength(), ownedScalarBondMap->NumMyPoints() );
  TEST_ASSERT( state.getBondMultiVector()->Map().SameAs( *ownedScalarBondMap ) );

  FieldManager& fm = FieldManager::self();
  int coordinatesFieldId = fm.getFieldId("Coordinates");
  int elementIdFieldId = fm.getFieldId("Element_Id");
  int forceDensityFieldId = fm.getFieldId("Force_Density");
  int damageFieldId = fm.getFieldId("Damage");

  // check initialization of data to zero
  Teuchos::RCP<Epetra_Vector> coordinates = state.getData(coordinatesFieldId);
  for(int i=0 ; i<coordinates->MyLength() ; ++i)
      TEST_EQUALITY_CONST( (*coordinates)[i], 0.0 );

  // set some data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      ids[i] = i;
    // vector data
    Epetra_Vector& force = *(state.getData(forceDensityFieldId));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      force[i*3] = i*3;
      force[i*3+1] = i*3+1;
      force[i*3+2] = i*3+2;
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        bondDamage[firstPointInElement+j] = firstPointInElement+j;
    }
  }

  // check the data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    for(int i=0 ; i<ids.MyLength() ; ++i)
        TEST_FLOATING_EQUALITY(ids[i], (double)(i), 1.0e-14);
    // vector data
    Epetra_Vector& force = *(state.getData(forceDensityFieldId));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
        TEST_FLOATING_EQUALITY(force[i*3], (double)(i*3), 1.0e-14);
        TEST_FLOATING_EQUALITY(force[i*3+1], (double)(i*3+1), 1.0e-14);
        TEST_FLOATING_EQUALITY(force[i*3+2], (double)(i*3+2), 1.0e-14);
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        TEST_FLOATING_EQUALITY(bondDamage[firstPointInElement+j], (double)(firstPointInElement+j), 1.0e-14);
    }
  }
}

//! Create a three-point problem for testing.

PeridigmNS::State createThreePointProblem(Teuchos::RCP<Epetra_Comm> comm, Teuchos::RCP<Epetra_BlockMap> &overlapScalarPointMap, Teuchos::RCP<Epetra_BlockMap> &overlapVectorPointMap, Teuchos::RCP<Epetra_BlockMap> &ownedScalarBondMap, vector<int> &scalarPointFieldIds, vector<int> &vectorPointFieldIds, vector<int> &bondFieldIds)
{
  
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
  overlapScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(-1, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));
  // overlapVectorPointMap
  // used for positions, displacements, velocities and vector constitutive data
  elementSize = 3;
  overlapVectorPointMap = 
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
  ownedScalarBondMap = 
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], &bondElementSize[0], indexBase, *comm));

  PeridigmNS::State state;

  FieldManager& fm = FieldManager::self();

  // create a list of scalar field specs and allocate the data
  
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Element_Id") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Proc_Num") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Number_Of_Neighbors") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Lambda") );
  scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Surface_Correction_Factor") );
  state.allocatePointData(PeridigmField::SCALAR, scalarPointFieldIds, overlapScalarPointMap);
  
  // create a list of vector field specs and allocate the data
  
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density") );
  vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Contact_Force_Density") );
  state.allocatePointData(PeridigmField::VECTOR,vectorPointFieldIds, overlapVectorPointMap);
  
  // create a list of bond field specs and allocate the data
  
  bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage") );
  bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deviatoric_Plastic_Extension") );
  state.allocateBondData(bondFieldIds, ownedScalarBondMap);

  return state;
}

//! Create a State object for a three-point problem, check data storage and basic functionality.

TEUCHOS_UNIT_TEST(State, ThreePointTest) {


  Teuchos::RCP<Epetra_Comm> comm;

  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

  PeridigmNS::State state;
  Teuchos::RCP<Epetra_BlockMap> overlapScalarPointMap;
  Teuchos::RCP<Epetra_BlockMap> overlapVectorPointMap;
  Teuchos::RCP<Epetra_BlockMap> ownedScalarBondMap;
  vector<int> scalarPointFieldIds;
  vector<int> vectorPointFieldIds;
  vector<int> bondFieldIds;
  
  state = createThreePointProblem(comm, overlapScalarPointMap, overlapVectorPointMap, ownedScalarBondMap, scalarPointFieldIds, vectorPointFieldIds, bondFieldIds);
  
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::SCALAR)->NumVectors(), (int)scalarPointFieldIds.size() );
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::SCALAR)->MyLength(), overlapScalarPointMap->NumMyPoints() );
  TEST_ASSERT( state.getPointMultiVector(PeridigmField::SCALAR)->Map().SameAs( *overlapScalarPointMap ) );

  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::VECTOR)->NumVectors(), (int)vectorPointFieldIds.size() );
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::VECTOR)->MyLength(), overlapVectorPointMap->NumMyPoints() );
  TEST_ASSERT( state.getPointMultiVector(PeridigmField::VECTOR)->Map().SameAs( *overlapVectorPointMap ) );

  TEST_EQUALITY( state.getBondMultiVector()->NumVectors(), (int)bondFieldIds.size() );
  TEST_EQUALITY( state.getBondMultiVector()->MyLength(), ownedScalarBondMap->NumMyPoints() );
  TEST_ASSERT( state.getBondMultiVector()->Map().SameAs( *ownedScalarBondMap ) );

  FieldManager& fm = FieldManager::self();
  int coordinatesFieldId = fm.getFieldId("Coordinates");
  int elementIdFieldId = fm.getFieldId("Element_Id");
  int forceDensityFieldId = fm.getFieldId("Force_Density");
  int damageFieldId = fm.getFieldId("Damage");

  // check initialization of data to zero
  Teuchos::RCP<Epetra_Vector> coordinates = state.getData(coordinatesFieldId);
  for(int i=0 ; i<coordinates->MyLength() ; ++i)
    TEST_EQUALITY_CONST( (*coordinates)[i], 0.0 );

  // set some data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      ids[i] = i;
    // vector data
    Epetra_Vector& force = *(state.getData(forceDensityFieldId));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      force[i*3] = i*3;
      force[i*3+1] = i*3+1;
      force[i*3+2] = i*3+2;
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        bondDamage[firstPointInElement+j] = firstPointInElement+j;
    }
  }

  // check the data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      TEST_FLOATING_EQUALITY(ids[i], (double)(i), 1.0e-14);
    // vector data
    Epetra_Vector& force = *(state.getData(forceDensityFieldId));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      TEST_FLOATING_EQUALITY(force[i*3], (double)(i*3), 1.0e-14);
      TEST_FLOATING_EQUALITY(force[i*3+1], (double)(i*3+1), 1.0e-14);
      TEST_FLOATING_EQUALITY(force[i*3+2], (double)(i*3+2), 1.0e-14);
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
        TEST_FLOATING_EQUALITY(bondDamage[firstPointInElement+j], (double)(firstPointInElement+j), 1.0e-14);
    }
  }
}

//! Test ability to copy data from one State to another.

TEUCHOS_UNIT_TEST(State, CopyFrom) {

  Teuchos::RCP<Epetra_Comm> comm;

  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

  PeridigmNS::State state;
  Teuchos::RCP<Epetra_BlockMap> overlapScalarPointMap;
  Teuchos::RCP<Epetra_BlockMap> overlapVectorPointMap;
  Teuchos::RCP<Epetra_BlockMap> ownedScalarBondMap;
  vector<int> scalarPointFieldIds;
  vector<int> vectorPointFieldIds;
  vector<int> bondFieldIds;
  
  state = createThreePointProblem(comm, overlapScalarPointMap, overlapVectorPointMap, ownedScalarBondMap, scalarPointFieldIds, vectorPointFieldIds, bondFieldIds);
  
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::SCALAR)->NumVectors(), (int)scalarPointFieldIds.size() );
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::SCALAR)->MyLength(), overlapScalarPointMap->NumMyPoints() );
  TEST_ASSERT( state.getPointMultiVector(PeridigmField::SCALAR)->Map().SameAs( *overlapScalarPointMap ) );

  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::VECTOR)->NumVectors(), (int)vectorPointFieldIds.size() );
  TEST_EQUALITY( state.getPointMultiVector(PeridigmField::VECTOR)->MyLength(), overlapVectorPointMap->NumMyPoints() );
  TEST_ASSERT( state.getPointMultiVector(PeridigmField::VECTOR)->Map().SameAs( *overlapVectorPointMap ) );

  TEST_EQUALITY( state.getBondMultiVector()->NumVectors(), (int)bondFieldIds.size() );
  TEST_EQUALITY( state.getBondMultiVector()->MyLength(), ownedScalarBondMap->NumMyPoints() );
  TEST_ASSERT( state.getBondMultiVector()->Map().SameAs( *ownedScalarBondMap ) );

  FieldManager& fm = FieldManager::self();
  int elementIdFieldId = fm.getFieldId("Element_Id");
  int forceDensityFieldId = fm.getFieldId("Force_Density");
  int damageFieldId = fm.getFieldId("Damage");

  // set some data
  {
    // scalar data
    Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    for(int i=0 ; i<ids.MyLength() ; ++i)
      ids[i] = i;
    // vector data
    Epetra_Vector& force = *(state.getData(forceDensityFieldId));
    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      force[i*3] = i*3;
      force[i*3+1] = i*3+1;
      force[i*3+2] = i*3+2;
    }
    // bond data
    Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
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
  int procID = state.getPointMultiVector(PeridigmField::SCALAR)->Comm().MyPID();
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
  
  //PeridigmField::Relation relation;
  //PeridigmField::Length length;

  //allocate data for the field specs, which are taken from the initial state object

  //scalar point data
  //relation = PeridigmField::ELEMENT;
  //length = PeridigmField::SCALAR;
  //vector<int> scalarPointFieldIds = state.getFieldIds(relation, length);

  tempState.allocatePointData(PeridigmField::SCALAR, scalarPointFieldIds, tempOverlapScalarPointMap);
  TEST_EQUALITY( tempState.getPointMultiVector(PeridigmField::SCALAR)->NumVectors(), state.getPointMultiVector(PeridigmField::SCALAR)->NumVectors());
  TEST_EQUALITY( tempState.getPointMultiVector(PeridigmField::SCALAR)->MyLength(), tempOverlapScalarPointMap->NumMyPoints() );
  TEST_ASSERT( tempState.getPointMultiVector(PeridigmField::SCALAR)->Map().SameAs( *tempOverlapScalarPointMap ) );

  //vector point data
  //relation = PeridigmField::NODE;
  //length = PeridigmField::VECTOR;
  //vector<int> vectorFieldIds = state.getFieldIds(relation, length);

  tempState.allocatePointData(PeridigmField::VECTOR,vectorPointFieldIds, tempOverlapVectorPointMap);
  TEST_EQUALITY( tempState.getPointMultiVector(PeridigmField::VECTOR)->NumVectors(), state.getPointMultiVector(PeridigmField::VECTOR)->NumVectors());
  TEST_EQUALITY( tempState.getPointMultiVector(PeridigmField::VECTOR)->MyLength(), tempOverlapVectorPointMap->NumMyPoints() );
  TEST_ASSERT( tempState.getPointMultiVector(PeridigmField::VECTOR)->Map().SameAs( *tempOverlapVectorPointMap ) );

  //scalar bond data
  //relation = PeridigmField::BOND;
  //length = PeridigmField::SCALAR;
  //vector<int> bondFieldIds = state.getFieldIds(relation, length);

  tempState.allocateBondData(bondFieldIds, tempOwnedScalarBondMap);
  TEST_EQUALITY( tempState.getBondMultiVector()->NumVectors(), state.getBondMultiVector()->NumVectors());
  TEST_EQUALITY( tempState.getBondMultiVector()->MyLength(), tempOwnedScalarBondMap->NumMyPoints() );
  TEST_ASSERT( tempState.getBondMultiVector()->Map().SameAs( *tempOwnedScalarBondMap ) );

  // copy the data from the initial State to the smaller, temporary state
  tempState.copyLocallyOwnedDataFromState(Teuchos::RCP<PeridigmNS::State>(&state, false));

  // check the data
  vector<int> tempFieldIds = tempState.getFieldIds();
  for(unsigned int iId=0 ; iId<tempFieldIds.size() ; ++iId){
    Teuchos::RCP<Epetra_Vector> data = state.getData( tempFieldIds[iId] );
    Teuchos::RCP<Epetra_Vector> tempData = tempState.getData( tempFieldIds[iId] );
    for(int iLID=0 ; iLID<tempData->Map().NumMyElements() ; ++iLID){
      int globalID = tempData->Map().GID(iLID);
      TEST_ASSERT(globalID != -1);
      int stateLID = data->Map().LID(globalID);
      TEST_ASSERT(stateLID != -1);
      int tempElementSize = tempData->Map().ElementSize(iLID);
      int stateElementSize = data->Map().ElementSize(stateLID);
      TEST_EQUALITY( tempElementSize, stateElementSize );
      int tempFirstPointInElement = tempData->Map().FirstPointInElement(iLID);
      int stateFirstPointInElement = data->Map().FirstPointInElement(stateLID);
      for(int i=0 ; i<tempElementSize ; ++i){
        TEST_EQUALITY( (*tempData)[tempFirstPointInElement + i], (*data)[stateFirstPointInElement + i] );
      }
    }
  }
}





int main( int argc, char* argv[] ) {

    int numProcs = 1;

    int returnCode = -1;
   
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
   
    if(numProcs == 1 || numProcs == 2){
       returnCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
    }
    else{
       std::cerr << "Unit test runtime ERROR: utPeridigm_State only makes sense on 1 or 2 processors." << std::endl;
    }

    return returnCode;

}

