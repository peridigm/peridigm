/*! \file utPHAL_UpdateForceState_MPI_np2.cpp */

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
#include <Phalanx_DataLayout_MDALayout.hpp>
#include <Epetra_MpiComm.h>
#include "PHAL_FactoryTraits.hpp"
#include "PHAL_Dimension.hpp" // needed for Dummy Data Layout
#include "PHAL_UpdateForceState.hpp"
#include "../materials/Peridigm_LinearElasticIsotropicMaterial.hpp"

#include <mpi.h>

using namespace boost::unit_test;

void testInitialize()
{
  // set up a parameter list that will be passed to the evaluator constructor
  Teuchos::RCP<Teuchos::ParameterList> p = rcp(new Teuchos::ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_update_force_state;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", false);

  //! \todo Figure out the purpose of the dummy layout code.
  Teuchos::RCP<PHX::DataLayout> dummy = Teuchos::rcp(new PHX::MDALayout<Dummy>(0));
  p->set< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout", dummy);

  // instantiate the evaluator
  UpdateForceState<PHAL::PeridigmTraits::Residual, PHAL::PeridigmTraits> evaluator(*p);

  // assert the evaluator name
  BOOST_CHECK(evaluator.getName() == "UpdateForceState");

  // assert the vector of fields that the evaluator evaluates
  const std::vector<Teuchos::RCP< PHX::FieldTag > > & evaluatedFields = evaluator.evaluatedFields();
  BOOST_REQUIRE(evaluatedFields.size() == 1);
  BOOST_CHECK(evaluatedFields[0]->name() == "UpdateForceState");
  BOOST_CHECK(evaluatedFields[0]->dataLayout() == *dummy);

  // assert the vector of fields that the evaluator is dependent upon
  const std::vector<Teuchos::RCP< PHX::FieldTag > > & dependentFields = evaluator.dependentFields();
  BOOST_CHECK(dependentFields.size() == 0);
}

void testTwoPts_MPI_np2()
{
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // set up a hard-coded layout for two points
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  // node map
  int numGlobalElements = 2;
  int numMyElements;
  int* myGlobalElements = 0;
  int indexBase = 0;
  if(rank == 0){
	numMyElements = 1;
	myGlobalElements = new int[numMyElements];
	myGlobalElements[0] = 0;
  }
  else if(rank == 1){
	numMyElements = 1;
	myGlobalElements = new int[numMyElements];
	myGlobalElements[0] = 1;
  }
  Epetra_Map node_map(numGlobalElements, numMyElements, myGlobalElements, indexBase, comm); 
  if(myGlobalElements != 0){
	delete[] myGlobalElements;
	myGlobalElements = 0;
  }

  // overlap node map
  // for this problem it works out the the overlap map is the same for both processors
  numGlobalElements = -1;
  numMyElements = 2;
  myGlobalElements = new int[numMyElements];
  myGlobalElements[0] = 0;
  myGlobalElements[1] = 1;
  Epetra_Map node_overlap_map(numGlobalElements, numMyElements, myGlobalElements, indexBase, comm); 
  if(myGlobalElements != 0){
	delete[] myGlobalElements;
	myGlobalElements = 0;
  }

  // overlap unknown map
  numGlobalElements = 12;
  if(rank == 0){
	numMyElements = 6;
	myGlobalElements = new int[numMyElements];
	myGlobalElements[0] = 0;
	myGlobalElements[1] = 1;
	myGlobalElements[2] = 2;
	myGlobalElements[3] = 6;
	myGlobalElements[4] = 7;
	myGlobalElements[5] = 8;
  }
  else if(rank == 1){
	numMyElements = 6;
	myGlobalElements = new int[numMyElements];
	myGlobalElements[0] = 3;
	myGlobalElements[1] = 4;
	myGlobalElements[2] = 5;
	myGlobalElements[3] = 9;
	myGlobalElements[4] = 10;
	myGlobalElements[5] = 11;
  }
  Epetra_Map unknown_map(numGlobalElements, numMyElements, myGlobalElements, indexBase, comm); 
  if(myGlobalElements != 0){
	delete[] myGlobalElements;
	myGlobalElements = 0;
  }

  // overlap unknown map
  // for this problem it works out the the overlap map is the same for both processors
  numGlobalElements = -1;
  numMyElements = 12;
  myGlobalElements = new int[numMyElements];
  for(int i=0 ; i<12 ; ++i){
	myGlobalElements[i] = i;
  }
  Epetra_Map unknown_overlap_map(numGlobalElements, numMyElements, myGlobalElements, indexBase, comm); 
  if(myGlobalElements != 0){
	delete[] myGlobalElements;
	myGlobalElements = 0;
  }

  // create a linear elastic isotropic peridynamic solid  material model
  Teuchos::ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  Peridigm::LinearElasticIsotropicMaterial mat(params);

  // set up the vectors
  Epetra_Vector initial_x(unknown_map);
  Epetra_Vector x(unknown_map);
  Epetra_Vector cell_volume(node_map);
  int numStateVariables = mat.NumConstitutiveVariables();
  BOOST_CHECK(mat.NumConstitutiveVariables() == 2);
  Epetra_MultiVector force_state_data(node_map, numStateVariables);  

  // two-point discretization
  if(rank == 0){
	// this is the first point
	initial_x[0] = 0.0; initial_x[1] = 0.0; initial_x[2] = 0.0;
	x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
	cell_volume[0] = 1.0;
  }
  else if(rank == 1){
	// this is the second point
	initial_x[0] = 1.0; initial_x[1] = 0.0; initial_x[2] = 0.0;
	x[0] = 2.0; x[1] = 0.0; x[2] = 0.0;
	cell_volume[0] = 1.0;
  }
  
  force_state_data.PutScalar(0.0);

  // both points are neighbors of each other
  map< int, Teuchos::RCP<vector<int> > > neighbor_list;
  Teuchos::RCP<vector<int> > n_list = Teuchos::rcp(new vector<int>);
  if(rank == 0){
	n_list->push_back(1);
	neighbor_list[0] = n_list;
  }
  else if(rank == 1){
	n_list->push_back(0);
	neighbor_list[1] = n_list;
  }

  // create the importer and exporter objects
  Epetra_Import node_importer(node_overlap_map, node_map);
  Epetra_Import unknown_importer(unknown_overlap_map, unknown_map);
  Epetra_Export node_exporter(node_overlap_map, node_map);
  Epetra_Export unknown_exporter(unknown_overlap_map, unknown_map);

  // create and populate the overlap maps
  Epetra_Vector initial_x_overlap(unknown_overlap_map);
  initial_x_overlap.Import(initial_x, unknown_importer, Insert);
  Epetra_Vector x_overlap(unknown_overlap_map);
  x_overlap.Import(x, unknown_importer, Insert);
  Epetra_MultiVector force_state_data_overlap(node_overlap_map, numStateVariables);
  force_state_data_overlap.PutScalar(0.0);
  Epetra_Vector cell_volume_overlap(node_overlap_map);
  cell_volume_overlap.Import(cell_volume, node_importer, Insert);
  Epetra_Vector x_dot_overlap(unknown_overlap_map);
  x_dot_overlap.PutScalar(0.0);

  // create a workset with rcps to the relevant data
  PHAL::Workset workset;
  workset.initial_x_overlap = Teuchos::rcp(&initial_x_overlap, false);
  workset.x_overlap = Teuchos::rcp(&x_overlap, false);
  workset.force_state_data = Teuchos::rcp(&force_state_data, false);
  workset.force_state_data_overlap = Teuchos::rcp(&force_state_data_overlap, false);
  workset.cell_volume_overlap = Teuchos::rcp(&cell_volume_overlap, false);
  workset.x_dot_overlap = Teuchos::rcp(&x_dot_overlap, false);
  workset.node_importer = Teuchos::rcp(&node_importer, false);
  workset.node_exporter = Teuchos::rcp(&node_exporter, false);
  workset.neighbor_list = Teuchos::rcp(&neighbor_list, false);
  workset.materials.push_back(Teuchos::rcp(&mat, false));
  workset.myPID = comm.MyPID();

  // set up a parameter list that will be passed to the evaluator constructor
  Teuchos::RCP<Teuchos::ParameterList> p = rcp(new Teuchos::ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_update_force_state;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", false);
  Teuchos::RCP<PHX::DataLayout> dummy = Teuchos::rcp(new PHX::MDALayout<Dummy>(0));
  p->set< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout", dummy);

  // instantiate the evaluator
  UpdateForceState<PHAL::PeridigmTraits::Residual, PHAL::PeridigmTraits> evaluator(*p);

  // make a call to the evaluateFields() function in the evaluator
  // this is the workhorse function that calls the material model
  // and updates the values in the force_state_data vector
  evaluator.evaluateFields(workset);

  // the evaluator should have scattered the force_state_data
  // such that each processor has access to all the data
  
  // assert the weighted volumes
  double weighted_volume = force_state_data_overlap[0][0];
  BOOST_CHECK_CLOSE(weighted_volume, 1.0, 1.0e-15);
  weighted_volume = force_state_data_overlap[0][1];
  BOOST_CHECK_CLOSE(weighted_volume, 1.0, 1.0e-15);
  weighted_volume = force_state_data[0][0];
  BOOST_CHECK_CLOSE(weighted_volume, 1.0, 1.0e-15);

  // assert the dilatations
  double dilatation = force_state_data_overlap[1][0];
  BOOST_CHECK_CLOSE(dilatation, 3.0, 1.0e-15);
  dilatation = force_state_data_overlap[1][1];
  BOOST_CHECK_CLOSE(dilatation, 3.0, 1.0e-15);
  dilatation = force_state_data[1][0];
  BOOST_CHECK_CLOSE(dilatation, 3.0, 1.0e-15);
}

bool init_unit_test_suite()
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPHAL_UpdateForceState");
  proc->add(BOOST_TEST_CASE(&testInitialize));
  proc->add(BOOST_TEST_CASE(&testTwoPts_MPI_np2));
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
	std::cerr << "Unit test runtime ERROR: utPHAL_UpdateForceState_MPI_np2 only makes sense on 2 processors" << std::endl;
	std::cerr << "Re-run unit test $mpiexec -np 2 ./utPHAL_UpdateForceState_MPI_np2" << std::endl;
  }

  MPI_Finalize();

  return return_code;
}
