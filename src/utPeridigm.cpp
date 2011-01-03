/*! \file utPeridigm.cpp */

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
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif
#include "Peridigm.hpp"

using namespace boost::unit_test;
using namespace Teuchos;
using namespace PeridigmNS;

Teuchos::RCP<PeridigmNS::Peridigm> createTwoPointModel()
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

void initialize()
{
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createTwoPointModel();

  // \todo Write meaningful asserts.
}

//! This is a one-dimensional rebalance test; the rebalance should have no effect.
void rebalance()
{
  Teuchos::RCP<PeridigmNS::Peridigm> peridigm = createTwoPointModel();

  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalMap()->NumMyElements(), 2);
  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalMap()->ElementSize(), 3);
  Epetra_Vector initialX(*peridigm->getX());
  Epetra_Vector initialU(*peridigm->getU());
  Epetra_Vector initialY(*peridigm->getY());
  Epetra_Vector initialV(*peridigm->getV());
  Epetra_Vector initialA(*peridigm->getA());
  Epetra_Vector initialForce(*peridigm->getForce());

  peridigm->rebalance();

  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalMap()->NumMyElements(), 2);
  BOOST_CHECK_EQUAL(peridigm->getThreeDimensionalMap()->ElementSize(), 3);
  for(int i=0 ; i<initialX.MyLength(); ++i){
    BOOST_CHECK_CLOSE(initialX[i], (*peridigm->getX())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialU[i], (*peridigm->getU())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialY[i], (*peridigm->getY())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialV[i], (*peridigm->getV())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialA[i], (*peridigm->getA())[i], 1.0e-15);
    BOOST_CHECK_CLOSE(initialForce[i], (*peridigm->getForce())[i], 1.0e-15);
  }

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success = true;

	test_suite* proc = BOOST_TEST_SUITE("utPeridigm");
	proc->add(BOOST_TEST_CASE(&initialize));
	proc->add(BOOST_TEST_CASE(&rebalance));
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
  #ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
  #endif

  // Initialize UTF
  return unit_test_main(init_unit_test, argc, argv);

  #ifdef HAVE_MPI
    MPI_Finalize();
  #endif
}
