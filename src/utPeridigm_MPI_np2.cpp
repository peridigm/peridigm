/*! \file utPeridigm_MPI_np2.cpp */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? 
// David J. Littlewood   djlittl@sandia.gov 
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#include <Epetra_MpiComm.h>
#include "Peridigm.hpp"

using namespace boost::unit_test;
using namespace Teuchos;
using namespace PeridigmNS;

Teuchos::RCP<PeridigmNS::Peridigm> createTwoPointModel()
{
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  Teuchos::RCP<Epetra_Comm> comm;
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

  // set up parameter lists
  // these data would normally be read from an input xml file
  Teuchos::RCP<Teuchos::ParameterList> peridigmParams = rcp(new Teuchos::ParameterList());

  // problem parameters
  ParameterList& problemParams = peridigmParams->sublist("Problem");
  problemParams.set("Verbose", false);

  // material parameters
  ParameterList& materialParams = problemParams.sublist("Material");
  ParameterList& linearElasticMaterialParams = materialParams.sublist("Linear Elastic");
  linearElasticMaterialParams.set("Density", 7800.0);
  linearElasticMaterialParams.set("Bulk Modulus", 130.0e9);
  linearElasticMaterialParams.set("Shear Modulus", 78.0e9);

  // discretization parameters
  ParameterList& discretizationParams = problemParams.sublist("Discretization");
  discretizationParams.set("Type", "PdQuickGrid");
  discretizationParams.set("Horizon", 2.0);

  // pdQuickGrid tensor product mesh generator parameters
  ParameterList& pdQuickGridParams = discretizationParams.sublist("TensorProduct3DMeshGenerator");
  pdQuickGridParams.set("Type", "PdQuickGrid");
  pdQuickGridParams.set("X Origin", -2.0);
  pdQuickGridParams.set("Y Origin", -0.5);
  pdQuickGridParams.set("Z Origin", -0.5);
  pdQuickGridParams.set("X Length",  4.0);
  pdQuickGridParams.set("Y Length",  1.0);
  pdQuickGridParams.set("Z Length",  1.0);
  pdQuickGridParams.set("Number Points X", 2);
  pdQuickGridParams.set("Number Points Y", 1);
  pdQuickGridParams.set("Number Points Z", 1);

  // boundary conditions
  ParameterList& bcParams = problemParams.sublist("Boundary Conditions");
  // node sets
  // these sets associate a name with a list of node ids, stored as a string
  // in this case there's only one node per set
  bcParams.set("Min X Node Set", "0");
  bcParams.set("Max X Node Set", "1");
  // initial velocity boundary conditions
  // each boundary condition is associated with a node set, defined above
  ParameterList& initialVelocityMinXFace = bcParams.sublist("Initial Velocity Min X Face");
  initialVelocityMinXFace.set("Type", "Initial Velocity");
  initialVelocityMinXFace.set("Node Set", "Min X Node Set");
  initialVelocityMinXFace.set("Coordinate", "x");
  initialVelocityMinXFace.set("Value", -1.0);
  ParameterList& initialVelocityMaxXFace = bcParams.sublist("Initial Velocity Max X Face");
  initialVelocityMaxXFace.set("Type", "Initial Velocity");
  initialVelocityMaxXFace.set("Node Set", "Max X Node Set");
  initialVelocityMaxXFace.set("Coordinate", "x");
  initialVelocityMaxXFace.set("Value", 1.0);

  // solver parameters
  ParameterList& solverParams = peridigmParams->sublist("Solver");
  solverParams.set("Verbose", "false");
  ParameterList& verletParams = solverParams.sublist("Verlet");
  verletParams.set("Initial Time", 0.0);
  verletParams.set("Final Time", 1.0);
  verletParams.set("Fixed dt", 1.0);

  // create the Peridigm object
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = rcp(new PeridigmNS::Peridigm(comm, peridigmParams));

  return peridigm;
}

Teuchos::RCP<PeridigmNS::Peridigm> createEightPointModel()
{
  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

  // set up parameter lists
  // these data would normally be read from an input xml file
  Teuchos::RCP<Teuchos::ParameterList> peridigmParams = rcp(new Teuchos::ParameterList());

  // problem parameters
  ParameterList& problemParams = peridigmParams->sublist("Problem");
  problemParams.set("Verbose", false);

  // material parameters
  ParameterList& materialParams = problemParams.sublist("Material");
  ParameterList& linearElasticMaterialParams = materialParams.sublist("Linear Elastic");
  linearElasticMaterialParams.set("Density", 7800.0);
  linearElasticMaterialParams.set("Bulk Modulus", 130.0e9);
  linearElasticMaterialParams.set("Shear Modulus", 78.0e9);

  // discretization parameters
  ParameterList& discretizationParams = problemParams.sublist("Discretization");
  discretizationParams.set("Type", "PdQuickGrid");
  discretizationParams.set("Horizon", 2.0);

  // pdQuickGrid tensor product mesh generator parameters
  ParameterList& pdQuickGridParams = discretizationParams.sublist("TensorProduct3DMeshGenerator");
  pdQuickGridParams.set("Type", "PdQuickGrid");
  pdQuickGridParams.set("X Origin", -2.0);
  pdQuickGridParams.set("Y Origin", -2.0);
  pdQuickGridParams.set("Z Origin", -2.0);
  pdQuickGridParams.set("X Length",  4.0);
  pdQuickGridParams.set("Y Length",  4.0);
  pdQuickGridParams.set("Z Length",  4.0);
  pdQuickGridParams.set("Number Points X", 2);
  pdQuickGridParams.set("Number Points Y", 2);
  pdQuickGridParams.set("Number Points Z", 2);

  // solver parameters
  ParameterList& solverParams = peridigmParams->sublist("Solver");
  solverParams.set("Verbose", "false");
  ParameterList& verletParams = solverParams.sublist("Verlet");
  verletParams.set("Initial Time", 0.0);
  verletParams.set("Final Time", 1.0);
  verletParams.set("Fixed dt", 1.0);

  // create the Peridigm object
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = rcp(new PeridigmNS::Peridigm(comm, peridigmParams));

  return peridigm;
}

void initialize()
{
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createTwoPointModel();

  BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->NumMyElements(), 1);
  BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->ElementSize(), 1);
  BOOST_CHECK_EQUAL(peridigm->getOneDimensionalOverlapMap()->NumMyElements(), 2);
  BOOST_CHECK_EQUAL(peridigm->getOneDimensionalOverlapMap()->ElementSize(), 1);
  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalMap()->NumMyElements(), 1);
  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalMap()->ElementSize(), 3);
  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalOverlapMap()->NumMyElements(), 2);
  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalOverlapMap()->ElementSize(), 3);
  BOOST_CHECK_EQUAL(peridigm->getBondMap()->NumMyElements(), 1);

  // \todo Write additional asserts
}

void rebalanceTwoPointModel()
{
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createTwoPointModel();

  // Make copies of everything so that we can identify any changes
  // that might occur during rebalance (there should be none)
  Epetra_BlockMap oneDimensionalMap(*peridigm->getOneDimensionalMap());
  Epetra_BlockMap oneDimensionalOverlapMap(*peridigm->getOneDimensionalOverlapMap());
  Epetra_BlockMap threeDimensionalMap(*peridigm->getThreeDimensionalMap());
  Epetra_BlockMap threeDimensionalOverlapMap(*peridigm->getThreeDimensionalOverlapMap());
  Epetra_BlockMap bondMap(*peridigm->getBondMap());
  Epetra_Vector initialX(*peridigm->getX());
  Epetra_Vector initialU(*peridigm->getU());
  Epetra_Vector initialY(*peridigm->getY());
  Epetra_Vector initialV(*peridigm->getV());
  Epetra_Vector initialA(*peridigm->getA());
  Epetra_Vector initialForce(*peridigm->getForce());
  Epetra_Vector volume(*peridigm->getDataManager()->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE));
  Epetra_Vector coord3d(*peridigm->getDataManager()->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE));
  Epetra_Vector weightedVolume(*peridigm->getDataManager()->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE));
  Epetra_Vector displ3dN(*peridigm->getDataManager()->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector displ3dNP1(*peridigm->getDataManager()->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector curcoord3dN(*peridigm->getDataManager()->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector curcoord3dNP1(*peridigm->getDataManager()->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector veloc3dN(*peridigm->getDataManager()->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector veloc3dNP1(*peridigm->getDataManager()->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector force3dN(*peridigm->getDataManager()->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector force3dNP1(*peridigm->getDataManager()->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector dilatationN(*peridigm->getDataManager()->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector dilatationNP1(*peridigm->getDataManager()->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector damageN(*peridigm->getDataManager()->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector damageNP1(*peridigm->getDataManager()->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector bondDamageN(*peridigm->getDataManager()->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector bondDamageNP1(*peridigm->getDataManager()->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1));
  PeridigmNS::NeighborhoodData neighborhoodData(*peridigm->getNeighborhoodData());
  //PeridigmNS::NeighborhoodData contactNeighborhoodData(*peridigm->getContactNeighborhoodData());

  // call the rebalance function, which should produce to changes in serial
  peridigm->rebalance();

  // check everything to make sure nothing changed
  // check maps
  BOOST_CHECK(peridigm->getOneDimensionalMap()->SameAs(oneDimensionalMap));
  BOOST_CHECK(peridigm->getOneDimensionalOverlapMap()->SameAs(oneDimensionalOverlapMap));
  BOOST_CHECK(peridigm->getThreeDimensionalMap()->SameAs(threeDimensionalMap));
  BOOST_CHECK(peridigm->getThreeDimensionalOverlapMap()->SameAs(threeDimensionalOverlapMap));
  BOOST_CHECK(peridigm->getBondMap()->SameAs(bondMap));
  // check mothership vectors
  for(int i=0 ; i<initialX.MyLength(); ++i){
    BOOST_CHECK_CLOSE(initialX[i], (*peridigm->getX())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialU[i], (*peridigm->getU())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialY[i], (*peridigm->getY())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialV[i], (*peridigm->getV())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialA[i], (*peridigm->getA())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialForce[i], (*peridigm->getForce())[i], 1.0e-15);
  }
  // check data in DataManager
  Teuchos::RCP<PeridigmNS::DataManager> dataManager = peridigm->getDataManager();
  for(int i=0 ; i<dataManager->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(volume[i], (*dataManager->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(coord3d[i], (*dataManager->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(weightedVolume[i], (*dataManager->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(displ3dN[i], (*dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(displ3dNP1[i], (*dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(curcoord3dN[i], (*dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(curcoord3dNP1[i], (*dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(veloc3dN[i], (*dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(veloc3dNP1[i], (*dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(force3dN[i], (*dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(force3dNP1[i], (*dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(dilatationN[i], (*dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(dilatationNP1[i], (*dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(damageN[i], (*dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(damageNP1[i], (*dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(bondDamageN[i], (*dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(bondDamageNP1[i], (*dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  // check neighborhood data
  BOOST_CHECK_EQUAL(neighborhoodData.NumOwnedPoints(), peridigm->getNeighborhoodData()->NumOwnedPoints());
  BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodListSize(), peridigm->getNeighborhoodData()->NeighborhoodListSize());
  for(int i=0 ; i<peridigm->getNeighborhoodData()->NumOwnedPoints() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.OwnedIDs()[i], peridigm->getNeighborhoodData()->OwnedIDs()[i]);
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodPtr()[i], peridigm->getNeighborhoodData()->NeighborhoodPtr()[i]);
  }
  for(int i=0 ; i<peridigm->getNeighborhoodData()->NeighborhoodListSize() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodList()[i], peridigm->getNeighborhoodData()->NeighborhoodList()[i]);
  }
}

void rebalanceEightPointModel()
{
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createEightPointModel();

  // Make copies of everything so that we can identify any changes
  // that might occur during rebalance (there should be none)
  Epetra_BlockMap oneDimensionalMap(*peridigm->getOneDimensionalMap());
  Epetra_BlockMap oneDimensionalOverlapMap(*peridigm->getOneDimensionalOverlapMap());
  Epetra_BlockMap threeDimensionalMap(*peridigm->getThreeDimensionalMap());
  Epetra_BlockMap threeDimensionalOverlapMap(*peridigm->getThreeDimensionalOverlapMap());
  Epetra_BlockMap bondMap(*peridigm->getBondMap());
  Epetra_Vector initialX(*peridigm->getX());
  Epetra_Vector initialU(*peridigm->getU());
  Epetra_Vector initialY(*peridigm->getY());
  Epetra_Vector initialV(*peridigm->getV());
  Epetra_Vector initialA(*peridigm->getA());
  Epetra_Vector initialForce(*peridigm->getForce());
  Epetra_Vector volume(*peridigm->getDataManager()->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE));
  Epetra_Vector coord3d(*peridigm->getDataManager()->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE));
  Epetra_Vector weightedVolume(*peridigm->getDataManager()->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE));
  Epetra_Vector displ3dN(*peridigm->getDataManager()->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector displ3dNP1(*peridigm->getDataManager()->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector curcoord3dN(*peridigm->getDataManager()->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector curcoord3dNP1(*peridigm->getDataManager()->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector veloc3dN(*peridigm->getDataManager()->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector veloc3dNP1(*peridigm->getDataManager()->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector force3dN(*peridigm->getDataManager()->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector force3dNP1(*peridigm->getDataManager()->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector dilatationN(*peridigm->getDataManager()->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector dilatationNP1(*peridigm->getDataManager()->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector damageN(*peridigm->getDataManager()->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector damageNP1(*peridigm->getDataManager()->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1));
  Epetra_Vector bondDamageN(*peridigm->getDataManager()->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_N));
  Epetra_Vector bondDamageNP1(*peridigm->getDataManager()->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1));
  PeridigmNS::NeighborhoodData neighborhoodData(*peridigm->getNeighborhoodData());
  //PeridigmNS::NeighborhoodData contactNeighborhoodData(*peridigm->getContactNeighborhoodData());

  // call the rebalance function, which should produce to changes in serial
  peridigm->rebalance();

  // check everything to make sure nothing changed
  // check maps
  BOOST_CHECK(peridigm->getOneDimensionalMap()->SameAs(oneDimensionalMap));
  BOOST_CHECK(peridigm->getOneDimensionalOverlapMap()->SameAs(oneDimensionalOverlapMap));
  BOOST_CHECK(peridigm->getThreeDimensionalMap()->SameAs(threeDimensionalMap));
  BOOST_CHECK(peridigm->getThreeDimensionalOverlapMap()->SameAs(threeDimensionalOverlapMap));
  BOOST_CHECK(peridigm->getBondMap()->SameAs(bondMap));
  // check mothership vectors
  for(int i=0 ; i<initialX.MyLength(); ++i){
    BOOST_CHECK_CLOSE(initialX[i], (*peridigm->getX())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialU[i], (*peridigm->getU())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialY[i], (*peridigm->getY())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialV[i], (*peridigm->getV())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialA[i], (*peridigm->getA())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialForce[i], (*peridigm->getForce())[i], 1.0e-15);
  }
  // check data in DataManager
  Teuchos::RCP<PeridigmNS::DataManager> dataManager = peridigm->getDataManager();
  for(int i=0 ; i<dataManager->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(volume[i], (*dataManager->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(coord3d[i], (*dataManager->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(weightedVolume[i], (*dataManager->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(displ3dN[i], (*dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(displ3dNP1[i], (*dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(curcoord3dN[i], (*dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(curcoord3dNP1[i], (*dataManager->getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(veloc3dN[i], (*dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(veloc3dNP1[i], (*dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(force3dN[i], (*dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(force3dNP1[i], (*dataManager->getData(Field_NS::FORCE3D, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(dilatationN[i], (*dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(dilatationNP1[i], (*dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(damageN[i], (*dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(damageNP1[i], (*dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_N)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(bondDamageN[i], (*dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_N))[i], 1.0e-15);
  for(int i=0 ; i<dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1)->MyLength() ; ++i)
    BOOST_CHECK_CLOSE(bondDamageNP1[i], (*dataManager->getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1))[i], 1.0e-15);
  // check neighborhood data
  BOOST_CHECK_EQUAL(neighborhoodData.NumOwnedPoints(), peridigm->getNeighborhoodData()->NumOwnedPoints());
  BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodListSize(), peridigm->getNeighborhoodData()->NeighborhoodListSize());
  for(int i=0 ; i<peridigm->getNeighborhoodData()->NumOwnedPoints() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.OwnedIDs()[i], peridigm->getNeighborhoodData()->OwnedIDs()[i]);
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodPtr()[i], peridigm->getNeighborhoodData()->NeighborhoodPtr()[i]);
  }
  for(int i=0 ; i<peridigm->getNeighborhoodData()->NeighborhoodListSize() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodList()[i], peridigm->getNeighborhoodData()->NeighborhoodList()[i]);
  }
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success = true;

	test_suite* proc = BOOST_TEST_SUITE("utPeridigm__MPI_np2");
	proc->add(BOOST_TEST_CASE(&initialize));
    proc->add(BOOST_TEST_CASE(&rebalanceTwoPointModel));
    proc->add(BOOST_TEST_CASE(&rebalanceEightPointModel));
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
	std::cerr << "Unit test runtime ERROR: utPeridigm_MPI_np2 only makes sense on 2 processors" << std::endl;
	std::cerr << "Re-run unit test $mpiexec -np 2 ./utPeridigm_MPI_np2" << std::endl;
  }

  MPI_Finalize();

  return return_code;
}
