/*! \file utPHAL_EvaluateForce_MPI_np2.cpp */

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
#include <Phalanx_DataLayout_MDALayout.hpp>
#include <Epetra_MpiComm.h>
#include "PHAL_FactoryTraits.hpp"
#include "PHAL_Dimension.hpp" // needed for Dummy Data Layout
#include "PHAL_EvaluateForce.hpp"
#include "../materials/Peridigm_LinearElasticIsotropicMaterial.hpp"

#include <mpi.h>


TEUCHOS_UNIT_TEST(PHPD_EvaluateForce_MPI_np2, TestInitializeTest) 

{
  // set up a parameter list that will be passed to the evaluator constructor
  Teuchos::RCP<Teuchos::ParameterList> p = rcp(new Teuchos::ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_evaluate_force;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", false);

  //! \todo Figure out the purpose of the dummy layout code.
  Teuchos::RCP<PHX::DataLayout> dummy = Teuchos::rcp(new PHX::MDALayout<Dummy>(0));
  p->set< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout", dummy);

  // instantiate the evaluator
  EvaluateForce<PHAL::PeridigmTraits::Residual, PHAL::PeridigmTraits> evaluator(*p);

  // assert the evaluator name
  TEST_ASSERT(evaluator.getName() == "EvaluateForce");

  // assert the vector of fields that the evaluator evaluates
  const std::vector<Teuchos::RCP< PHX::FieldTag > > & evaluatedFields = evaluator.evaluatedFields();
  TEST_ASSERT(evaluatedFields.size() == 1);
  TEST_ASSERT(evaluatedFields[0]->name() == "EvaluateForce");
  TEST_ASSERT(evaluatedFields[0]->dataLayout() == *dummy);

  // assert the vector of fields that the evaluator is dependent upon
  const std::vector<Teuchos::RCP< PHX::FieldTag > > & dependentFields = evaluator.dependentFields();
  TEST_ASSERT(dependentFields.size() == 1);
  TEST_ASSERT(dependentFields[0]->name() == "UpdateForceState");
  TEST_ASSERT(dependentFields[0]->dataLayout() == *dummy);
}


TEUCHOS_UNIT_TEST(PHPD_EvaluateForce_MPI_np2, TestTwoPts_MPI_np2Test) 

{
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // set up a hard-coded layout for two points
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  TEST_COMPARE(numProcs, ==, 2);
  if(numProcs != 2){
    std::cerr << "Unit test runtime ERROR: utPHAL_EvaluateForce_MPI_np2 only makes sense on 2 processors" << std::endl;
    return;
  }


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
  TEST_ASSERT(mat.NumConstitutiveVariables() == 2);
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

  // create place for the evaluator to put the solution
  Epetra_Vector x_dot(unknown_map);
  x_dot.PutScalar(-1.0);

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
  workset.x_dot = Teuchos::rcp(&x_dot, false);
  workset.x_dot_overlap = Teuchos::rcp(&x_dot_overlap, false);
  workset.node_importer = Teuchos::rcp(&node_importer, false);
  workset.node_exporter = Teuchos::rcp(&node_exporter, false);
  workset.unknown_exporter = Teuchos::rcp(&unknown_exporter, false);
  workset.neighbor_list = Teuchos::rcp(&neighbor_list, false);
  workset.materials.push_back(Teuchos::rcp(&mat, false));
  workset.myPID = comm.MyPID();

  // fill in the force_state_data directly, as opposed to calling
  // the UpdateForceState evaluator

  // it turns out that things are the same on both processors,
  // but this is not true in general so handle each seperately

  if(rank == 0){
	// weighted volume
	force_state_data[0][0] = 1.0;
	force_state_data_overlap[0][0] = 1.0;
	force_state_data_overlap[0][1] = 1.0;
	// dilatation
	force_state_data[1][0] = 3.0;
	force_state_data_overlap[1][0] = 3.0;
	force_state_data_overlap[1][1] = 3.0;
  }
  else if(rank == 1){
	// weighted volume
	force_state_data[0][0] = 1.0;
	force_state_data_overlap[0][0] = 1.0;
	force_state_data_overlap[0][1] = 1.0;
	// dilatation
	force_state_data[1][0] = 3.0;
	force_state_data_overlap[1][0] = 3.0;
	force_state_data_overlap[1][1] = 3.0;
  }

  // set up a parameter list that will be passed to the evaluator constructor
  Teuchos::RCP<Teuchos::ParameterList> p = rcp(new Teuchos::ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_evaluate_force;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", false);
  Teuchos::RCP<PHX::DataLayout> dummy = Teuchos::rcp(new PHX::MDALayout<Dummy>(0));
  p->set< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout", dummy);

  // instantiate the evaluator
  EvaluateForce<PHAL::PeridigmTraits::Residual, PHAL::PeridigmTraits> evaluator(*p);

  // make a call to the evaluateFields() function in the evaluator
  // this is the workhorse function that calls the material model,
  // evaluates the pairwise forces, and updates x_dot_overlap
  evaluator.evaluateFields(workset);

  // assert the data in x_dot_overlap
  double node0_velocity_x_overlap = x_dot_overlap[0];
  TEST_COMPARE(node0_velocity_x_overlap, <=, 1.0e-14);
  double node0_velocity_y_overlap = x_dot_overlap[1];
  TEST_COMPARE(node0_velocity_y_overlap, <=,  1.0e-14);
  double node0_velocity_z_overlap = x_dot_overlap[2];
  TEST_COMPARE(node0_velocity_z_overlap, <=, 1.0e-14);
  double node1_velocity_x_overlap = x_dot_overlap[3];
  TEST_COMPARE(node1_velocity_x_overlap,  <=, 1.0e-14);
  double node1_velocity_y_overlap = x_dot_overlap[4];
  TEST_COMPARE(node1_velocity_y_overlap, <=, 1.0e-14);
  double node1_velocity_z_overlap = x_dot_overlap[5];
  TEST_COMPARE(node1_velocity_z_overlap, <=, 1.0e-14);
  double node0_acceleration_x_overlap = x_dot_overlap[6];
  TEST_FLOATING_EQUALITY(node0_acceleration_x_overlap, 1.5e8, 1.0e-14);
  double node0_acceleration_y_overlap = x_dot_overlap[7];
  TEST_COMPARE(node0_acceleration_y_overlap, 1.0e-14);
  double node0_acceleration_z_overlap = x_dot_overlap[8];
  TEST_COMPARE(node0_acceleration_z_overlap, <=, 1.0e-14);
  double node1_acceleration_x_overlap = x_dot_overlap[9];
  TEST_FLOATING_EQUALITY(node1_acceleration_x_overlap, -1.5e8, 1.0e-14);
  double node1_acceleration_y_overlap = x_dot_overlap[10];
  TEST_COMPARE(node1_acceleration_y_overlap, 1.0e-14);
  double node1_acceleration_z_overlap = x_dot_overlap[11];
  TEST_COMPARE(node1_acceleration_z_overlap, 1.0e-14);

  if(rank == 0){
	double node0_velocity_x = x_dot[0];
	TEST_COMPARE(node0_velocity_x, <=, 1.0e-14);
	double node0_velocity_y = x_dot[1];
	TEST_COMPARE(node0_velocity_y, <= , 1.0e-14);
	double node0_velocity_z = x_dot[2];
	TEST_COMPARE(node0_velocity_z, <=, 1.0e-14);
	double node0_acceleration_x = x_dot[3];
	TEST_FLOATING_EQUALITY(node0_acceleration_x, 3.0e8, 1.0e-14);
	double node0_acceleration_y = x_dot[4];
	TEST_COMPARE(node0_acceleration_y, <=, 1.0e-14);
	double node0_acceleration_z = x_dot[5];
	TEST_COMPARE(node0_acceleration_z, <=, 1.0e-14);
  }
  else if(rank == 1){
	double node1_velocity_x = x_dot[0];
	TEST_COMPARE(node1_velocity_x, <=, 1.0e-14);
	double node1_velocity_y = x_dot[1];
	TEST_COMPARE(node1_velocity_y, <=, 1.0e-14);
	double node1_velocity_z = x_dot[2];
	TEST_COMPARE(node1_velocity_z, 1.0e-14);
	double node1_acceleration_x = x_dot[3];
	TEST_FLOATING_EQUALITY(node1_acceleration_x, -3.0e8, 1.0e-14);
	double node1_acceleration_y = x_dot[4];
        TEST_COMPARE(node1_acceleration_y, <=, 1.0e-14);
	double node1_acceleration_z = x_dot[5];
	TEST_COMPARE(node1_acceleration_z, <=,1.0e-14);
  }
}


int main
(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
