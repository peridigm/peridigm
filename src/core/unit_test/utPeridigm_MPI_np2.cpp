/*! \file utPeridigm_MPI_np2.cpp */

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

  // discretization parameters
  ParameterList& discretizationParams = peridigmParams->sublist("Discretization");
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

  // material parameters
  ParameterList& materialParams = peridigmParams->sublist("Materials");
  ParameterList& linearElasticMaterialParams = materialParams.sublist("My Linear Elastic Material");
  linearElasticMaterialParams.set("Material Model", "Linear Elastic");
  linearElasticMaterialParams.set("Density", 7800.0);
  linearElasticMaterialParams.set("Bulk Modulus", 130.0e9);
  linearElasticMaterialParams.set("Shear Modulus", 78.0e9);

  // block parameters
  ParameterList& blockParams = peridigmParams->sublist("Blocks");
  ParameterList& blockGroupParams = blockParams.sublist("My Group of Blocks");
  blockGroupParams.set("Block Names", "block_1");
  blockGroupParams.set("Material", "My Linear Elastic Material");

  // boundary conditions
  ParameterList& bcParams = peridigmParams->sublist("Boundary Conditions");
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
  initialVelocityMinXFace.set("Value", "-1.0");
  ParameterList& initialVelocityMaxXFace = bcParams.sublist("Initial Velocity Max X Face");
  initialVelocityMaxXFace.set("Type", "Initial Velocity");
  initialVelocityMaxXFace.set("Node Set", "Max X Node Set");
  initialVelocityMaxXFace.set("Coordinate", "x");
  initialVelocityMaxXFace.set("Value", "1.0");

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

  // discretization parameters
  ParameterList& discretizationParams = peridigmParams->sublist("Discretization");
  discretizationParams.set("Type", "PdQuickGrid");
  discretizationParams.set("Horizon", 4.1);

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

  // material parameters
  ParameterList& materialParams = peridigmParams->sublist("Materials");
  ParameterList& linearElasticMaterialParams = materialParams.sublist("My Linear Elastic Material");
  linearElasticMaterialParams.set("Material Model", "Linear Elastic");
  linearElasticMaterialParams.set("Density", 7800.0);
  linearElasticMaterialParams.set("Bulk Modulus", 130.0e9);
  linearElasticMaterialParams.set("Shear Modulus", 78.0e9);

  // block parameters
  ParameterList& blockParams = peridigmParams->sublist("Blocks");
  ParameterList& blockGroupParams = blockParams.sublist("My Group of Blocks");
  blockGroupParams.set("Block Names", "block_1");
  blockGroupParams.set("Material", "My Linear Elastic Material");

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
  Epetra_BlockMap bondMap(*peridigm->getBondMap());
  Epetra_Vector initialX(*peridigm->getX());
  Epetra_Vector initialU(*peridigm->getU());
  Epetra_Vector initialY(*peridigm->getY());
  Epetra_Vector initialV(*peridigm->getV());
  Epetra_Vector initialA(*peridigm->getA());
  Epetra_Vector initialForce(*peridigm->getForce());
  Epetra_Vector volume( *peridigm->getBlock(0)->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE) );
  Epetra_Vector coord3d( *peridigm->getBlock(0)->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE) );
  Epetra_Vector weightedVolume( *peridigm->getBlock(0)->getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE) );
  Epetra_Vector displ3dN( *peridigm->getBlock(0)->getData(Field_NS::DISPL3D, Field_ENUM::STEP_N) );
  Epetra_Vector displ3dNP1( *peridigm->getBlock(0)->getData(Field_NS::DISPL3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector curcoord3dN( *peridigm->getBlock(0)->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_N) );
  Epetra_Vector curcoord3dNP1( *peridigm->getBlock(0)->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector veloc3dN( *peridigm->getBlock(0)->getData(Field_NS::VELOC3D, Field_ENUM::STEP_N) );
  Epetra_Vector veloc3dNP1( *peridigm->getBlock(0)->getData(Field_NS::VELOC3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector force3dN( *peridigm->getBlock(0)->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_N) );
  Epetra_Vector force3dNP1( *peridigm->getBlock(0)->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector dilatationN( *peridigm->getBlock(0)->getData(Field_NS::DILATATION, Field_ENUM::STEP_N) );
  Epetra_Vector dilatationNP1( *peridigm->getBlock(0)->getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1) );
  Epetra_Vector damageN( *peridigm->getBlock(0)->getData(Field_NS::DAMAGE, Field_ENUM::STEP_N) );
  Epetra_Vector damageNP1( *peridigm->getBlock(0)->getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1) );
  Epetra_Vector bondDamageN( *peridigm->getBlock(0)->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_N) );
  Epetra_Vector bondDamageNP1( *peridigm->getBlock(0)->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1) );
  PeridigmNS::NeighborhoodData neighborhoodData(*peridigm->getGlobalNeighborhoodData() );
  //PeridigmNS::NeighborhoodData contactNeighborhoodData(*peridigm->getContactNeighborhoodData());

  // call the rebalance function, which should produce no changes
  peridigm->rebalance();

  // check everything to make sure nothing changed
  // check maps
  BOOST_CHECK(peridigm->getOneDimensionalMap()->SameAs(oneDimensionalMap));
  BOOST_CHECK(peridigm->getOneDimensionalOverlapMap()->SameAs(oneDimensionalOverlapMap));
  BOOST_CHECK(peridigm->getThreeDimensionalMap()->SameAs(threeDimensionalMap));
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
  // check field data
  Teuchos::RCP<PeridigmNS::Block> block = peridigm->getBlock(0);
  for(int i=0 ; i<block->numPoints() ; ++i){
    BOOST_CHECK_CLOSE(volume[i], (*block->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(coord3d[i], (*block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(weightedVolume[i], (*block->getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(displ3dN[i], (*block->getData(Field_NS::DISPL3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(displ3dNP1[i], (*block->getData(Field_NS::DISPL3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(curcoord3dN[i], (*block->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(curcoord3dNP1[i], (*block->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(veloc3dN[i], (*block->getData(Field_NS::VELOC3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(veloc3dNP1[i], (*block->getData(Field_NS::VELOC3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(force3dN[i], (*block->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(force3dNP1[i], (*block->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(dilatationN[i], (*block->getData(Field_NS::DILATATION, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(dilatationNP1[i], (*block->getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(damageN[i], (*block->getData(Field_NS::DAMAGE, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(damageNP1[i], (*block->getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(bondDamageN[i], (*block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(bondDamageNP1[i], (*block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1))[i], 1.0e-15);
  }
  // check neighborhood data
  BOOST_CHECK_EQUAL(neighborhoodData.NumOwnedPoints(), peridigm->getGlobalNeighborhoodData()->NumOwnedPoints());
  BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodListSize(), peridigm->getGlobalNeighborhoodData()->NeighborhoodListSize());
  for(int i=0 ; i<peridigm->getGlobalNeighborhoodData()->NumOwnedPoints() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.OwnedIDs()[i], peridigm->getGlobalNeighborhoodData()->OwnedIDs()[i]);
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodPtr()[i], peridigm->getGlobalNeighborhoodData()->NeighborhoodPtr()[i]);
  }
  for(int i=0 ; i<peridigm->getGlobalNeighborhoodData()->NeighborhoodListSize() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodList()[i], peridigm->getGlobalNeighborhoodData()->NeighborhoodList()[i]);
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
  Epetra_BlockMap bondMap(*peridigm->getBondMap());
  Epetra_Vector initialX(*peridigm->getX());
  Epetra_Vector initialU(*peridigm->getU());
  Epetra_Vector initialY(*peridigm->getY());
  Epetra_Vector initialV(*peridigm->getV());
  Epetra_Vector initialA(*peridigm->getA());
  Epetra_Vector initialForce(*peridigm->getForce());
  Epetra_Vector volume( *peridigm->getBlock(0)->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE) );
  Epetra_Vector coord3d( *peridigm->getBlock(0)->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE) );
  Epetra_Vector weightedVolume( *peridigm->getBlock(0)->getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE) );
  Epetra_Vector displ3dN( *peridigm->getBlock(0)->getData(Field_NS::DISPL3D, Field_ENUM::STEP_N) );
  Epetra_Vector displ3dNP1( *peridigm->getBlock(0)->getData(Field_NS::DISPL3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector curcoord3dN( *peridigm->getBlock(0)->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_N) );
  Epetra_Vector curcoord3dNP1( *peridigm->getBlock(0)->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector veloc3dN( *peridigm->getBlock(0)->getData(Field_NS::VELOC3D, Field_ENUM::STEP_N) );
  Epetra_Vector veloc3dNP1( *peridigm->getBlock(0)->getData(Field_NS::VELOC3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector force3dN( *peridigm->getBlock(0)->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_N) );
  Epetra_Vector force3dNP1( *peridigm->getBlock(0)->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1) );
  Epetra_Vector dilatationN( *peridigm->getBlock(0)->getData(Field_NS::DILATATION, Field_ENUM::STEP_N) );
  Epetra_Vector dilatationNP1( *peridigm->getBlock(0)->getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1) );
  Epetra_Vector damageN( *peridigm->getBlock(0)->getData(Field_NS::DAMAGE, Field_ENUM::STEP_N) );
  Epetra_Vector damageNP1( *peridigm->getBlock(0)->getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1) );
  Epetra_Vector bondDamageN( *peridigm->getBlock(0)->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_N) );
  Epetra_Vector bondDamageNP1( *peridigm->getBlock(0)->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1) );
  PeridigmNS::NeighborhoodData neighborhoodData(*peridigm->getGlobalNeighborhoodData() );
  //PeridigmNS::NeighborhoodData contactNeighborhoodData(*peridigm->getContactNeighborhoodData());

  // call the rebalance function, which should produce no changes
  peridigm->rebalance();

  // check everything to make sure nothing changed
  // check maps
  BOOST_CHECK(peridigm->getOneDimensionalMap()->SameAs(oneDimensionalMap));
  BOOST_CHECK(peridigm->getOneDimensionalOverlapMap()->SameAs(oneDimensionalOverlapMap));
  BOOST_CHECK(peridigm->getThreeDimensionalMap()->SameAs(threeDimensionalMap));
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
  // check field data
  Teuchos::RCP<PeridigmNS::Block> block = peridigm->getBlock(0);
  for(int i=0 ; i<block->numPoints() ; ++i){
    BOOST_CHECK_CLOSE(volume[i], (*block->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(coord3d[i], (*block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(weightedVolume[i], (*block->getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(displ3dN[i], (*block->getData(Field_NS::DISPL3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(displ3dNP1[i], (*block->getData(Field_NS::DISPL3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(curcoord3dN[i], (*block->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(curcoord3dNP1[i], (*block->getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(veloc3dN[i], (*block->getData(Field_NS::VELOC3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(veloc3dNP1[i], (*block->getData(Field_NS::VELOC3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(force3dN[i], (*block->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(force3dNP1[i], (*block->getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(dilatationN[i], (*block->getData(Field_NS::DILATATION, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(dilatationNP1[i], (*block->getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(damageN[i], (*block->getData(Field_NS::DAMAGE, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(damageNP1[i], (*block->getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(bondDamageN[i], (*block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_N))[i], 1.0e-15);
    BOOST_CHECK_CLOSE(bondDamageNP1[i], (*block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1))[i], 1.0e-15);
  }
  // check neighborhood data
  BOOST_CHECK_EQUAL(neighborhoodData.NumOwnedPoints(), peridigm->getGlobalNeighborhoodData()->NumOwnedPoints());
  BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodListSize(), peridigm->getGlobalNeighborhoodData()->NeighborhoodListSize());
  for(int i=0 ; i<peridigm->getGlobalNeighborhoodData()->NumOwnedPoints() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.OwnedIDs()[i], peridigm->getGlobalNeighborhoodData()->OwnedIDs()[i]);
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodPtr()[i], peridigm->getGlobalNeighborhoodData()->NeighborhoodPtr()[i]);
  }
  for(int i=0 ; i<peridigm->getGlobalNeighborhoodData()->NeighborhoodListSize() ; ++i){
    BOOST_CHECK_EQUAL(neighborhoodData.NeighborhoodList()[i], peridigm->getGlobalNeighborhoodData()->NeighborhoodList()[i]);
  }
}

void rebalanceEightPointModelSwitchCorners()
{
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createEightPointModel();

  // make sure the points have ended up where we expect them to be
  // there is more than one 'correct' answer here, but for the purpose of this test we want
  // to check the decomposition to see that it's the same as when the test was set up
  if(rank == 0){
    // global ID 0
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(0), 0);
    BOOST_CHECK_CLOSE((*peridigm->getX())[0], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[2], -1.0, 1.0e-15);
    // global ID 2
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(2), 1);
    BOOST_CHECK_CLOSE((*peridigm->getX())[3], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[4], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[5], -1.0, 1.0e-15);
    // global ID 4
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(4), 2);
    BOOST_CHECK_CLOSE((*peridigm->getX())[6], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[7], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[8], 1.0, 1.0e-15);
    // global ID 6
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(6), 3);
    BOOST_CHECK_CLOSE((*peridigm->getX())[9], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[11], 1.0, 1.0e-15);
  }
  else if(rank == 1){
    // global ID 5
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(5), 0);
    BOOST_CHECK_CLOSE((*peridigm->getX())[0], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[2], 1.0, 1.0e-15);
    // global ID 7
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(7), 1);
    BOOST_CHECK_CLOSE((*peridigm->getX())[3], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[4], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[5], 1.0, 1.0e-15);
    // global ID 1
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(1), 2);
    BOOST_CHECK_CLOSE((*peridigm->getX())[6], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[7], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[8], -1.0, 1.0e-15);
    // global ID 3
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(3), 3);
    BOOST_CHECK_CLOSE((*peridigm->getX())[9], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[11], -1.0, 1.0e-15);
  }

  // switch the positions of the points with global IDs 2 and 7
  // initial position of global ID 2 is (-1, 1, -1)
  // initial position of global ID 7 is ( 1, 1,  1)
  Teuchos::RCP<const Epetra_Vector> x = peridigm->getX();
  Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp_const_cast<Epetra_Vector>(peridigm->getU());
  Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp_const_cast<Epetra_Vector>(peridigm->getY());
  if(rank == 0){
    // displacement of global ID 2
    (*u)[3] = 2.0;
    (*u)[4] = 0.0;
    (*u)[5] = 2.0;
  }
  else if(rank == 1){
    // displacement of global ID 7
    (*u)[3] = -2.0;
    (*u)[4] = 0.0;
    (*u)[5] = -2.0;
  }
  // update y with new displacement
  for(int i=0 ; i<y->MyLength() ; ++i)
    (*y)[i] = (*x)[i] + (*u)[i];

  // set BOND_DAMAGE to indicate broken bonds
  // this will allow bond data rebalance to be checked
  // \todo The following is a great illustration of why the neighborhood data structure needs work...
  Teuchos::RCP<PeridigmNS::Block> block = peridigm->getBlock(0);
  Teuchos::RCP<Epetra_Vector> bondDamageN = block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_N);
  Teuchos::RCP<Epetra_Vector> bondDamageNP1 = block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1);
  if(rank == 0){
    // break the second bond for the point with global ID 2
    Teuchos::RCP<const PeridigmNS::NeighborhoodData> neighborhoodData = peridigm->getGlobalNeighborhoodData();
    int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    int* ownedIDs = neighborhoodData->OwnedIDs();
    int* neighborhoodList = neighborhoodData->NeighborhoodList();
    int neighborhoodIndex = 0;
    int dataIndex = 0;
    for(int i=0 ; i<numOwnedPoints ; ++i){
      int localID = ownedIDs[i];
      int globalID = peridigm->getOneDimensionalOverlapMap()->GID(localID);
      int numNeighbors = neighborhoodList[neighborhoodIndex++];
      BOOST_CHECK_EQUAL(numNeighbors, 7);
      for(int j=0 ; j<numNeighbors ; ++j){
        if(globalID == 2 && j == 1){
          // break the bond for the NP1 state
          (*bondDamageNP1)[dataIndex] = 1.0;
        }
        neighborhoodIndex++;
        dataIndex++;
      }
    }
  }
  else if(rank == 1){
    // break the seventh bond for the point with global ID 7
    Teuchos::RCP<const PeridigmNS::NeighborhoodData> neighborhoodData = peridigm->getGlobalNeighborhoodData();
    int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    int* ownedIDs = neighborhoodData->OwnedIDs();
    int* neighborhoodList = neighborhoodData->NeighborhoodList();
    int neighborhoodIndex = 0;
    int dataIndex = 0;
    for(int i=0 ; i<numOwnedPoints ; ++i){
      int localID = ownedIDs[i];
      int globalID = peridigm->getOneDimensionalOverlapMap()->GID(localID);
      int numNeighbors = neighborhoodList[neighborhoodIndex++];
      BOOST_CHECK_EQUAL(numNeighbors, 7);
      for(int j=0 ; j<numNeighbors ; ++j){
        if(globalID == 7 && j == 6){
          // break the bond for both states
          (*bondDamageN)[dataIndex] = 1.0;
          (*bondDamageNP1)[dataIndex] = 1.0;
        }
        neighborhoodIndex++;
        dataIndex++;
      }
    }
  }

  // before rebalance the global IDs are distributed as follows:
  // processor 0: 0 2 4 6
  // processor 1: 5 7 1 3

  // call the rebalance function
  peridigm->rebalance();

  // after rebalance the global IDs are distributed as follows:
  // processor 0: 0 4 6 7
  // processor 1: 5 1 3 2

  // check global IDs
  // the points with global IDs 2 and 7 should be swapped relative to where they started
  if(rank == 0){
    // global ID 0
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(0), 0);
    BOOST_CHECK_CLOSE((*peridigm->getX())[0], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[2], -1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[0], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[1], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[2], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[0], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[2], -1.0, 1.0e-15);
    // global ID 4
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(4), 1);
    BOOST_CHECK_CLOSE((*peridigm->getX())[3], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[4], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[5], 1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[3], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[4], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[5], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[3], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[4], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[5], 1.0, 1.0e-15);
    // global ID 6
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(6), 2);
    BOOST_CHECK_CLOSE((*peridigm->getX())[6], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[7], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[8], 1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[6], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[7], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[8], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[6], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[7], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[8], 1.0, 1.0e-15);
    // global ID 7 (should be where global ID 2 was originally)
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(7), 3);
    BOOST_CHECK_CLOSE((*peridigm->getX())[9], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[11], 1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[9], -2.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[10], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[11], -2.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[9], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[11], -1.0, 1.0e-15);
  }
  else if(rank == 1){
    // global ID 5
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(5), 0);
    BOOST_CHECK_CLOSE((*peridigm->getX())[0], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[2], 1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[0], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[1], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[2], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[0], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[2], 1.0, 1.0e-15);
    // global ID 1
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(1), 1);
    BOOST_CHECK_CLOSE((*peridigm->getX())[3], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[4], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[5], -1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[3], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[4], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[5], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[3], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[4], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[5], -1.0, 1.0e-15);
    // global ID 3
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(3), 2);
    BOOST_CHECK_CLOSE((*peridigm->getX())[6], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[7], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[8], -1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[6], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[7], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[8], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[6], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[7], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[8], -1.0, 1.0e-15);
    // global ID 2 (should be where global ID 7 was originally)
    BOOST_CHECK_EQUAL(peridigm->getOneDimensionalMap()->LID(2), 3);
    BOOST_CHECK_CLOSE((*peridigm->getX())[9], -1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getX())[11], -1.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getU())[9], 2.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[10], 0.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getU())[11], 2.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getY())[9], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE((*peridigm->getY())[11], 1.0, 1.0e-15);
  }

  // check acceleration and force mothership vectors
  for(int i=0 ; i<peridigm->getX()->MyLength(); ++i){
    BOOST_CHECK_CLOSE((*peridigm->getV())[i], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getA())[i], 0.0, 1.0e-15);
    BOOST_CHECK_CLOSE((*peridigm->getForce())[i], 0.0, 1.0e-15);
  }
  // check field data
  if(rank == 0){
    // length of the overlap vectors should be 8*3 = 24, all the off-processor points are ghosted
    BOOST_CHECK_EQUAL(block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->MyLength(), 24);
    double* coord3d;
    block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&coord3d);
    // global ID 0
    BOOST_CHECK_CLOSE(coord3d[0], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[2], -1.0, 1.0e-15);
    // global ID 4
    BOOST_CHECK_CLOSE(coord3d[3], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[4], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[5], 1.0, 1.0e-15);
    // global ID 6
    BOOST_CHECK_CLOSE(coord3d[6], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[7], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[8], 1.0, 1.0e-15);
    // global ID 7
    BOOST_CHECK_CLOSE(coord3d[9], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[11], 1.0, 1.0e-15);
    // global ID 1, ghosted
    BOOST_CHECK_CLOSE(coord3d[12], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[13], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[14], -1.0, 1.0e-15);    
    // global ID 2, ghosted
    BOOST_CHECK_CLOSE(coord3d[15], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[16], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[17], -1.0, 1.0e-15);    
    // global ID 3, ghosted
    BOOST_CHECK_CLOSE(coord3d[18], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[19], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[20], -1.0, 1.0e-15);    
    // global ID 5, ghosted
    BOOST_CHECK_CLOSE(coord3d[21], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[22], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[23], 1.0, 1.0e-15);    
  }
  else if(rank == 1){
    // length of the overlap vectors should be 8*3 = 24, all the off-processor points are ghosted
    BOOST_CHECK_EQUAL(block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->MyLength(), 24);
    double* coord3d;
    block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&coord3d);
    // global ID 5
    BOOST_CHECK_CLOSE(coord3d[0], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[1], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[2], 1.0, 1.0e-15);
    // global ID 1
    BOOST_CHECK_CLOSE(coord3d[3], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[4], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[5], -1.0, 1.0e-15);
    // global ID 3
    BOOST_CHECK_CLOSE(coord3d[6], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[7], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[8], -1.0, 1.0e-15);
    // global ID 2
    BOOST_CHECK_CLOSE(coord3d[9], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[10], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[11], -1.0, 1.0e-15);
    // global ID 0, ghosted
    BOOST_CHECK_CLOSE(coord3d[12], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[13], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[14], -1.0, 1.0e-15);    
    // global ID 4, ghosted
    BOOST_CHECK_CLOSE(coord3d[15], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[16], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[17], 1.0, 1.0e-15);    
    // global ID 6, ghosted
    BOOST_CHECK_CLOSE(coord3d[18], -1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[19], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[20], 1.0, 1.0e-15);    
    // global ID 7, ghosted
    BOOST_CHECK_CLOSE(coord3d[21], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[22], 1.0, 1.0e-15); BOOST_CHECK_CLOSE(coord3d[23], 1.0, 1.0e-15);    
  }

  // check BOND_DAMAGE
  bondDamageN = block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_N);
  bondDamageNP1 = block->getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1);
  if(rank == 0){
    // all bonds should be intact except for the seventh bond for the point with global ID 7.
    Teuchos::RCP<const PeridigmNS::NeighborhoodData> neighborhoodData = peridigm->getGlobalNeighborhoodData();
    int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    int* ownedIDs = neighborhoodData->OwnedIDs();
    int* neighborhoodList = neighborhoodData->NeighborhoodList();
    int neighborhoodIndex = 0;
    int dataIndex = 0;
    for(int i=0 ; i<numOwnedPoints ; ++i){
      int localID = ownedIDs[i];
      int globalID = peridigm->getOneDimensionalOverlapMap()->GID(localID);
      int numNeighbors = neighborhoodList[neighborhoodIndex++];
      BOOST_CHECK_EQUAL(numNeighbors, 7);
      for(int j=0 ; j<numNeighbors ; ++j){
        if(globalID == 7 && j == 6){
          // bond should be broken for both states
          BOOST_CHECK_CLOSE((*bondDamageN)[dataIndex], 1.0, 1.0e-14);
          BOOST_CHECK_CLOSE((*bondDamageNP1)[dataIndex], 1.0, 1.0e-14);
        }
        else{
          BOOST_CHECK_CLOSE((*bondDamageN)[dataIndex], 0.0, 1.0e-14);
          BOOST_CHECK_CLOSE((*bondDamageNP1)[dataIndex], 0.0, 1.0e-14);
        }
        neighborhoodIndex++;
        dataIndex++;
      }
    }
  }
  else if(rank == 1){
    // all bonds should be intact except for the second bond for the point with global ID 2.
    Teuchos::RCP<const PeridigmNS::NeighborhoodData> neighborhoodData = peridigm->getGlobalNeighborhoodData();
    int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    int* ownedIDs = neighborhoodData->OwnedIDs();
    int* neighborhoodList = neighborhoodData->NeighborhoodList();
    int neighborhoodIndex = 0;
    int dataIndex = 0;
    for(int i=0 ; i<numOwnedPoints ; ++i){
      int localID = ownedIDs[i];
      int globalID = peridigm->getOneDimensionalOverlapMap()->GID(localID);
      int numNeighbors = neighborhoodList[neighborhoodIndex++];
      BOOST_CHECK_EQUAL(numNeighbors, 7);
      for(int j=0 ; j<numNeighbors ; ++j){
        if(globalID == 2 && j == 1){
          // bond should be broken for NP1 state
          BOOST_CHECK_CLOSE((*bondDamageN)[dataIndex], 0.0, 1.0e-14);
          BOOST_CHECK_CLOSE((*bondDamageNP1)[dataIndex], 1.0, 1.0e-14);
        }
        else{
          BOOST_CHECK_CLOSE((*bondDamageN)[dataIndex], 0.0, 1.0e-14);
          BOOST_CHECK_CLOSE((*bondDamageNP1)[dataIndex], 0.0, 1.0e-14);
        }
        neighborhoodIndex++;
        dataIndex++;
      }
    }
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
    proc->add(BOOST_TEST_CASE(&rebalanceEightPointModelSwitchCorners));
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
