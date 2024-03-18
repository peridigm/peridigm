/*! \file utPeridigm_MultiphysicsElasticMaterial.cpp */

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
#include "Teuchos_UnitTestRepository.hpp"
#include "Peridigm_MultiphysicsElasticMaterial.hpp"
#include "Peridigm_ElasticMaterial.hpp"
#include "Peridigm_SerialMatrix.hpp"
#include "Peridigm_Field.hpp"
#include <Epetra_SerialComm.h>

using namespace PeridigmNS;
using namespace Teuchos;

//! Tests state variable count and name accessor functions.

TEUCHOS_UNIT_TEST(MultiphysicsElasticMaterial, testStateVariableAccessors) {

  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("Horizon", 10.0);
	params.set("Apply Automatic Differentiation Jacobian", true);
  params.set("Permeability", 1.0);
	params.set("Fluid density", 1000.0);
	params.set("Fluid dynamic viscosity", 1.0);
	params.set("Fluid compressibility", 1.0);
	params.set("Fluid Reynolds viscosity temperature effect", 0.0);
	params.set("Fluid linear thermal expansion", 0.0);
	params.set("Permeability curve inflection damage", .50);
	params.set("Permeability alpha", .25);
	params.set("Max permeability", 1.e-4);

  MultiphysicsElasticMaterial mat(params);
  // \todo check field specs
}

//! Tests two-point system under compression against hand calculations.

TEUCHOS_UNIT_TEST(MultiphysicsElasticMaterial, testTwoPts) {

  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("Horizon", 10.0);
	params.set("Apply Automatic Differentiation Jacobian", true);
  params.set("Permeability", 1.0);
	params.set("Fluid density", 1000.0);
	params.set("Fluid dynamic viscosity", 1.0);
	params.set("Fluid compressibility", 1.0);
	params.set("Fluid Reynolds viscosity temperature effect", 0.0);
	params.set("Fluid linear thermal expansion", 0.0);
	params.set("Permeability curve inflection damage", .50);
	params.set("Permeability alpha", .25);
	params.set("Max permeability", 1.e-4);

  MultiphysicsElasticMaterial mat(params);

	
  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(2, 0, comm);
  Epetra_Map unknownMap(6, 0, comm);
  Epetra_Map bondMap(2, 0, comm);
  double dt = 1.0;
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;
  // \todo check field specs

  // set up discretization
  numOwnedPoints = 2;
  ownedIDs = new int[numOwnedPoints];
  ownedIDs[0] = 0;
  ownedIDs[1] = 1;
  neighborhoodList = new int[4];
  neighborhoodList[0] = 1;
  neighborhoodList[1] = 1;
  neighborhoodList[2] = 1;
  neighborhoodList[3] = 0;

  // create the data manager
  // in serial, the overlap and non-overlap maps are the same
  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&bondMap, false));
  dataManager.allocateData(mat.FieldIds());

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  int fluidPressureYFieldId = fieldManager.getFieldId("Fluid_Pressure_Y");
  int volumeFieldId = fieldManager.getFieldId("Volume");
  int weightedVolumeFieldId = fieldManager.getFieldId("Weighted_Volume");
  int dilatationFieldId = fieldManager.getFieldId("Dilatation");
  int bondDamageFieldId = fieldManager.getFieldId("Bond_Damage");
  int forceDensityFieldId = fieldManager.getFieldId("Force_Density");
//  int flowDensityFieldId = fieldManager.getFieldId("Flux_Density");

  Epetra_Vector& x = *dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE);
  Epetra_Vector& y = *dataManager.getData(coordinatesFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& fluidPressureY = *dataManager.getData(fluidPressureYFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& cellVolume = *dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE);

  x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
  x[3] = 1.0; x[4] = 0.0; x[5] = 0.0;
  y[0] = 0.0; y[1] = 0.0; y[2] = 0.0;
  y[3] = 2.0; y[4] = 0.0; y[5] = 0.0;

  // current fluid pressure is set to zero so that the solid mechanics model is not interfered with
  for(int i=0 ; i<fluidPressureY.MyLength() ; ++i){
    fluidPressureY[i] = 0.0;
  }
	  
  for(int i=0; i<cellVolume.MyLength(); ++i){
	cellVolume[i] = 1.0;
  }

	  
  mat.initialize(dt, 
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 dataManager);
  
  mat.computeForce(dt, 
				   numOwnedPoints,
				   ownedIDs,
				   neighborhoodList,
                   dataManager);

  double currentPositionX1 = y[0];
  TEST_COMPARE(currentPositionX1, <=, 1.0e-14);
  double currentPositionY1 = y[1];
  TEST_COMPARE(currentPositionY1, <=, 1.0e-14);
  double currentPositionZ1 = y[2];
  TEST_COMPARE(currentPositionZ1, <=, 1.0e-14);
  double currentPositionX2 = y[3];
  TEST_FLOATING_EQUALITY(currentPositionX2, 2.0, 1.0e-12);
  double currentPositionY2 = y[4];
  TEST_COMPARE(currentPositionY2, <=, 1.0e-14);
  double currentPositionZ2 = y[5];
  TEST_COMPARE(currentPositionZ2, <=,  1.0e-14);
  Epetra_Vector& weightedVolume = *dataManager.getData(weightedVolumeFieldId, PeridigmField::STEP_NONE);
  TEST_FLOATING_EQUALITY(weightedVolume[0], 1.0, 1.0e-15);
  TEST_FLOATING_EQUALITY(weightedVolume[1], 1.0, 1.0e-15);
  Epetra_Vector& dilatation = *dataManager.getData(dilatationFieldId, PeridigmField::STEP_NP1);
  TEST_FLOATING_EQUALITY(dilatation[0], 3.0, 1.0e-15);
  TEST_FLOATING_EQUALITY(dilatation[1], 3.0, 1.0e-15);
  Epetra_Vector& bondDamage = *dataManager.getData(bondDamageFieldId, PeridigmField::STEP_NP1);
  TEST_COMPARE(bondDamage[0], <=, 1.0e-15);
  TEST_COMPARE(bondDamage[1], <=, 1.0e-15);

  Epetra_Vector& force = *dataManager.getData(forceDensityFieldId, PeridigmField::STEP_NP1);
  TEST_FLOATING_EQUALITY(force[0], 2.34e+12, 1.0e-2);
  TEST_COMPARE(force[1], <=, 1.0e-14);
  TEST_COMPARE(force[2], <=, 1.0e-14);
  TEST_FLOATING_EQUALITY(force[3], -2.34e+12, 1.0e-2);
  TEST_COMPARE(force[4], <=, 1.0e-14);
  TEST_COMPARE(force[5], <=, 1.0e-14);

  delete[] ownedIDs;
  delete[] neighborhoodList;
}

//! Tests eight-cell block under compression against hand calculations.

TEUCHOS_UNIT_TEST(MultiphysicsElasticMaterial, testEightPts) {
  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("Horizon", 10.0);
	params.set("Apply Automatic Differentiation Jacobian", true);
  params.set("Permeability", 1.0);
	params.set("Fluid density", 1000.0);
	params.set("Fluid dynamic viscosity", 1.0);
	params.set("Fluid compressibility", 1.0);
	params.set("Fluid Reynolds viscosity temperature effect", 0.0);
	params.set("Fluid linear thermal expansion", 0.0);
	params.set("Permeability curve inflection damage", .50);
	params.set("Permeability alpha", .25);
	params.set("Max permeability", 1.e-4);

  MultiphysicsElasticMaterial mat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(8, 0, comm);
  Epetra_Map unknownMap(24, 0, comm);
  Epetra_Map bondMap(56, 0, comm); // total number of bonds = 8(7) = 56
  double dt = 1.0;
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;

  // set up discretization
  // all cells are neighbors of each other
  numOwnedPoints = 8;
  ownedIDs = new int[numOwnedPoints];
  for(int i=0 ; i<numOwnedPoints; ++i){
	ownedIDs[i] = i;
  }
  // the neighborhood list has the format
  // numNeighborsNode1 nID1 nID2 ... nIDn numNeighborsNode2 nID1 nID2 ... nIDn ...
  // there are 8 cells, each has 7 neighbors
  // total length of neighborhoodList = 8(1+7) = 64
  neighborhoodList = new int[64];
  int neighborhoodListIndex = 0;
  for(int i=0 ; i<numOwnedPoints; ++i){
	neighborhoodList[neighborhoodListIndex++] = 7;
	for(int j=0 ; j<8 ; ++j){
	  if(i != j){
		neighborhoodList[neighborhoodListIndex++] = j;
	  }
	}
  }

  // create the data manager
  // in serial, the overlap and non-overlap maps are the same
  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&bondMap, false));
  dataManager.allocateData(mat.FieldIds());
 
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
	int fluidPressureYFieldId = fieldManager.getFieldId("Fluid_Pressure_Y");
  int volumeFieldId = fieldManager.getFieldId("Volume");
  int weightedVolumeFieldId = fieldManager.getFieldId("Weighted_Volume");
  int dilatationFieldId = fieldManager.getFieldId("Dilatation");
  int bondDamageFieldId = fieldManager.getFieldId("Bond_Damage");
  int forceDensityFieldId = fieldManager.getFieldId("Force_Density");
//  int flowDensityFieldId = fieldManager.getFieldId("Flux_Density");

  Epetra_Vector& x = *dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE);
  Epetra_Vector& y = *dataManager.getData(coordinatesFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& fluidPressureY = *dataManager.getData(fluidPressureYFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& cellVolume = *dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE);

  // initial positions
  x[0]  = 0.0; x[1]  = 0.0; x[2]  = 0.0;
  x[3]  = 0.0; x[4]  = 0.0; x[5]  = 1.0;
  x[6]  = 0.0; x[7]  = 1.0; x[8]  = 0.0;
  x[9]  = 0.0; x[10] = 1.0; x[11] = 1.0;
  x[12] = 1.0; x[13] = 0.0; x[14] = 0.0;
  x[15] = 1.0; x[16] = 0.0; x[17] = 1.0;
  x[18] = 1.0; x[19] = 1.0; x[20] = 0.0;
  x[21] = 1.0; x[22] = 1.0; x[23] = 1.0;

  // current positions
  y[0]  = 0.0; y[1]  = 0.0; y[2]  = 0.0;
  y[3]  = 0.0; y[4]  = 0.0; y[5]  = 0.98;
  y[6]  = 0.0; y[7]  = 1.0; y[8]  = 0.0;
  y[9]  = 0.0; y[10] = 1.0; y[11] = 0.98;
  y[12] = 1.0; y[13] = 0.0; y[14] = 0.0;
  y[15] = 1.0; y[16] = 0.0; y[17] = 0.98;
  y[18] = 1.0; y[19] = 1.0; y[20] = 0.0;
  y[21] = 1.0; y[22] = 1.0; y[23] = 0.98;

	// current fluid pressures
	// they are set to zero so that the solid mechanics model is not interfered with
	for(int i=0; i<8; ++i){
		fluidPressureY[i] = 0.0;
	}

  // cell volumes
  for(int i=0; i<cellVolume.MyLength(); ++i){
	cellVolume[i] = 1.0;
  }

  mat.initialize(dt, 
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 dataManager);

  mat.computeForce(dt, 
				   numOwnedPoints,
				   ownedIDs,
				   neighborhoodList,
                   dataManager);

  double currentPosition;
  currentPosition = y[0];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[1];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[2];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[3];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[4];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[5];
  TEST_FLOATING_EQUALITY(currentPosition, 0.98, 1.0e-15);
  currentPosition = y[6];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[7];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[8];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[9];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[10];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[11];
  TEST_FLOATING_EQUALITY(currentPosition, 0.98, 1.0e-15);
  currentPosition = y[12];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[13];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[14];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[15];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[16];
  TEST_COMPARE(currentPosition, <=,  1.0e-14);
  currentPosition = y[17];
  TEST_FLOATING_EQUALITY(currentPosition, 0.98, 1.0e-15);
  currentPosition = y[18];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[19];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[20];
  TEST_COMPARE(currentPosition, <=, 1.0e-14);
  currentPosition = y[21];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[22];
  TEST_FLOATING_EQUALITY(currentPosition, 1.0, 1.0e-15);
  currentPosition = y[23];
  TEST_FLOATING_EQUALITY(currentPosition, 0.98, 1.0e-15);

  // the weighted volumes and dilatations are the
  // same for all points in this test problem
  Epetra_Vector& weightedVolume = *dataManager.getData(weightedVolumeFieldId, PeridigmField::STEP_NONE);
  Epetra_Vector& dilatation = *dataManager.getData(dilatationFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& bondDamage = *dataManager.getData(bondDamageFieldId, PeridigmField::STEP_NP1);
  for(int i=0; i<8; ++i){
	TEST_FLOATING_EQUALITY(weightedVolume[i], 12.0, 1.0e-12);
	TEST_FLOATING_EQUALITY(dilatation[i], -0.01991593994643333, 1.0e-12);
  }

  // the bond damage should be all zeros (no damage)
  for(int i=0 ; i<bondDamage.MyLength() ; ++i){
      TEST_COMPARE(bondDamage[i], <=, 1.0e-15);
  }

  // all the cells experience the same force magnitude, but
  // it is applied in different directions

  // the hand solution for the pairwise forces between cells is
  // t = (y-x)*9.75e10

  // in this model, the forces between cells are always
  // equal and opposite (hence the factor of 2 in the
  // calculations below)

  // calculate by hand the force on cells 0
  // this can then be adjusted for the other cells by
  // rotating the force vector
  
  double t, unit_vec[3], vec_mag;
  double f[3] = {0.0, 0.0, 0.0};

  // force on cell 0 that results from interaction with cell 1
  t = (0.98 - 1.0)*9.75e10;
  vec_mag = 0.98;
  unit_vec[0] = 0.0; unit_vec[1] = 0.0; unit_vec[2] = 0.98/vec_mag;
  f[0] += 2.0*unit_vec[0]*t; f[1] += 2.0*unit_vec[1]*t; f[2] += 2.0*unit_vec[2]*t;
  
  // force on cell 0 that results from interaction with cell 2
  t = 0.0;

  // force on cell 0 that results from interaction with cell 3
  t = (sqrt(1.0+0.98*0.98) - sqrt(2.0))*9.75e10;
  vec_mag = sqrt(1.0 + 0.98*0.98);
  unit_vec[0] = 0.0; unit_vec[1] = 1.0/vec_mag; unit_vec[2] = 0.98/vec_mag;
  f[0] += 2.0*unit_vec[0]*t; f[1] += 2.0*unit_vec[1]*t; f[2] += 2.0*unit_vec[2]*t;

  // force on cell 0 that results from interaction with cell 4
  t = 0.0;

  // force on cell 0 that results from interaction with cell 5
  t = (sqrt(1.0+0.98*0.98) - sqrt(2.0))*9.75e10;
  vec_mag = sqrt(1.0 + 0.98*0.98);
  unit_vec[0] = 1.0/vec_mag; unit_vec[1] = 0.0; unit_vec[2] = 0.98/vec_mag;
  f[0] += 2.0*unit_vec[0]*t; f[1] += 2.0*unit_vec[1]*t; f[2] += 2.0*unit_vec[2]*t;

  // force on cell 0 that results from interaction with cell 6
  t = 0.0;

  // force on cell 0 that results from interaction with cell 7
  t = (sqrt(2.0 + 0.98*0.98) - sqrt(3.0))*9.75e10;
  vec_mag = sqrt(2.0 + 0.98*0.98);
  unit_vec[0] = 1.0/vec_mag; unit_vec[1] = 1.0/vec_mag; unit_vec[2] = 0.98/vec_mag;
  f[0] += 2.0*unit_vec[0]*t; f[1] += 2.0*unit_vec[1]*t; f[2] += 2.0*unit_vec[2]*t;

  Epetra_Vector& force = *dataManager.getData(forceDensityFieldId, PeridigmField::STEP_NP1);

  // assert the net force on cell 0
  TEST_FLOATING_EQUALITY(force[0], f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[1], f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[2], f[2], 1.0e-11);

  // assert the net forces on the other cells
  // these forces differ only by direction

  // force on cell 1
  TEST_FLOATING_EQUALITY(force[3],       f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[4],       f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[5],  f[2], 1.0e-11);

  // force on cell 2
  TEST_FLOATING_EQUALITY(force[6],       f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[7],  f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[8],       f[2], 1.0e-11);

  // force on cell 3
  TEST_FLOATING_EQUALITY(force[9],       f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[10], f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[11], f[2], 1.0e-11);

  // force on cell 4
  TEST_FLOATING_EQUALITY(-1.0*force[12], f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[13],      f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[14],      f[2], 1.0e-11);

  // force on cell 5
  TEST_FLOATING_EQUALITY(-1.0*force[15], f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[16],      f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[17], f[2], 1.0e-11);

  // force on cell 6
  TEST_FLOATING_EQUALITY(-1.0*force[18], f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[19], f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(force[20],      f[2], 1.0e-11);

  // force on cell 7
  TEST_FLOATING_EQUALITY(-1.0*force[21], f[0], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[22], f[1], 1.0e-11);
  TEST_FLOATING_EQUALITY(-1.0*force[23], f[2], 1.0e-11);

  delete[] ownedIDs;
  delete[] neighborhoodList;
}

//! Tests arbitrary three-cell system under arbitrary deformation against hand calculations.
TEUCHOS_UNIT_TEST(MultiphysicsElasticMaterial, testThreePts) {

  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("Horizon", 10.0);
	params.set("Apply Automatic Differentiation Jacobian", true);
  params.set("Permeability", 1.0);
	params.set("Fluid density", 1000.0);
	params.set("Fluid dynamic viscosity", 1.0);
	params.set("Fluid compressibility", 1.0);
	params.set("Fluid Reynolds viscosity temperature effect", 0.0);
	params.set("Fluid linear thermal expansion", 0.0);
	params.set("Permeability curve inflection damage", .50);
	params.set("Permeability alpha", .25);
	params.set("Max permeability", 1.e-4);

  MultiphysicsElasticMaterial mat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(3, 0, comm);
  Epetra_Map unknownMap(9, 0, comm);
  Epetra_Map bondMap(6, 0, comm);    // total number of bonds = 3(2) = 6
  double dt = 1.0;
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;

  // set up discretization
  // all cells are neighbors of each other
  numOwnedPoints = 3;
  ownedIDs = new int[numOwnedPoints];
  for(int i=0 ; i<numOwnedPoints; ++i){
	ownedIDs[i] = i;
  }
  // the neighborhood list has the format
  // numNeighborsNode1 nID1 nID2 ... nIDn numNeighborsNode2 nID1 nID2 ... nIDn ...
  // there are 3 cells, each has 2 neighbors
  // total length of neighborhoodList = 3(1+2) = 9
  neighborhoodList = new int[9];
  int neighborhoodListIndex = 0;
  for(int i=0 ; i<numOwnedPoints; ++i){
	neighborhoodList[neighborhoodListIndex++] = 2;
	for(int j=0 ; j<3 ; ++j){
	  if(i != j){
		neighborhoodList[neighborhoodListIndex++] = j;
	  }
	}
  }

  // create the data manager
  // in serial, the overlap and non-overlap maps are the same
  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&bondMap, false));
  dataManager.allocateData(mat.FieldIds());

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
	int modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
	int fluidPressureYFieldId = fieldManager.getFieldId("Fluid_Pressure_Y");
  int volumeFieldId = fieldManager.getFieldId("Volume");
  int weightedVolumeFieldId = fieldManager.getFieldId("Weighted_Volume");
  int dilatationFieldId = fieldManager.getFieldId("Dilatation");
  int bondDamageFieldId = fieldManager.getFieldId("Bond_Damage");
  int forceDensityFieldId = fieldManager.getFieldId("Force_Density");
//  int flowDensityFieldId = fieldManager.getFieldId("Flux_Density");

  Epetra_Vector& x = *dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE);
  Epetra_Vector& y = *dataManager.getData(coordinatesFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& fluidPressureY = *dataManager.getData(fluidPressureYFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& cellVolume = *dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE);

  // initial positions
  x[0] =  1.1; x[1] = 2.6;  x[2] = -0.1;
  x[3] = -2.0; x[4] = 0.9;  x[5] = -0.3;
  x[6] =  0.0; x[7] = 0.01; x[8] =  1.8;

  // current positions
  y[0] = 1.2;  y[1] = 2.4;  y[2] = -0.1;
  y[3] = -1.9; y[4] = 0.7;  y[5] = -0.8;
  y[6] = 0.1;  y[7] = 0.21; y[8] =  1.6;

	// current fluid pressures
	for(int i=0; i<3; i++){
		fluidPressureY[i] = 0.0;
	}

  // cell volumes
  cellVolume[0] = 0.9;
  cellVolume[1] = 1.1;
  cellVolume[2] = 0.8;

  mat.initialize(dt, 
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 dataManager);

  mat.computeForce(dt, 
				   numOwnedPoints,
				   ownedIDs,
				   neighborhoodList,
           dataManager);

  double currentPosition;
  currentPosition = y[0];
  TEST_FLOATING_EQUALITY(currentPosition, 1.2, 1.0e-13);
  currentPosition = y[1];
  TEST_FLOATING_EQUALITY(currentPosition, 2.4, 1.0e-13);
  currentPosition = y[2];
  TEST_FLOATING_EQUALITY(currentPosition, -0.1, 1.0e-13);
  currentPosition = y[3];
  TEST_FLOATING_EQUALITY(currentPosition, -1.9, 1.0e-13);
  currentPosition = y[4];
  TEST_FLOATING_EQUALITY(currentPosition, 0.7, 1.0e-13);
  currentPosition = y[5];
  TEST_FLOATING_EQUALITY(currentPosition, -0.8, 1.0e-13);
  currentPosition = y[6];
  TEST_FLOATING_EQUALITY(currentPosition, 0.1, 1.0e-13);
  currentPosition = y[7];
  TEST_FLOATING_EQUALITY(currentPosition, 0.21, 1.0e-13);
  currentPosition = y[8];
  TEST_FLOATING_EQUALITY(currentPosition, 1.6, 1.0e-13);

  // check the weighted volume and dilatation
  // against hand calculations
  Epetra_Vector& weightedVolume = *dataManager.getData(weightedVolumeFieldId, PeridigmField::STEP_NONE);
  Epetra_Vector& dilatation = *dataManager.getData(dilatationFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& bondDamage = *dataManager.getData(bondDamageFieldId, PeridigmField::STEP_NP1);
  TEST_FLOATING_EQUALITY(weightedVolume[0], 23.016479999999931, 1.0e-12);
  TEST_FLOATING_EQUALITY(dilatation[0], -0.114127034572639, 1.0e-11);
  TEST_FLOATING_EQUALITY(weightedVolume[1],18.64767999999997, 1.0e-12);
  TEST_FLOATING_EQUALITY(dilatation[1], 0.08257537372985, 1.0e-11);
  TEST_FLOATING_EQUALITY(weightedVolume[2], 20.497599999999963, 1.0e-12);
  TEST_FLOATING_EQUALITY(dilatation[2], -0.12166177890553, 1.0e-11);

  // the bond state should be all zeros (no damage)
  for(int i=0 ; i<bondDamage.MyLength() ; ++i){
      TEST_COMPARE(bondDamage[i], <=,  1.0e-15);
  }

  Epetra_Vector& force = *dataManager.getData(forceDensityFieldId, PeridigmField::STEP_NP1);

  // check the net forces against hand calculations
  double ref_soln_x, ref_soln_y, ref_soln_z;

  // cell 0
  // force on cell 0 due to interaction with cell 1
  // t_0 < x_1 - x_0 > dV_1
  ref_soln_x = -2753550531.144094*cellVolume[1];
  ref_soln_y = -1510011581.595148*cellVolume[1];
  ref_soln_z = -621769474.7744727*cellVolume[1];
  // - t_1 < x_0 - x_1 > dV_1
  ref_soln_x -= 3398655528.6806289*cellVolume[1];
  ref_soln_y -= 1863778838.308732*cellVolume[1];
  ref_soln_z -= 767438345.18594845*cellVolume[1];
  // add force on cell 0 due to interaction with cell 2
  // t_0 < x_2 - x_0 > dV_2
  ref_soln_x += 7736514797.0571524*cellVolume[2];
  ref_soln_y += 15402697641.41379*cellVolume[2];
  ref_soln_z += -11956431959.088315*cellVolume[2];
  // - t_2 < x_0 - x_2 > dV_2
  ref_soln_x -= -8687228655.850972*cellVolume[2];
  ref_soln_y -= -17295482505.73966*cellVolume[2];
  ref_soln_z -= 13425717013.587864*cellVolume[2];
  // assert the values of net force on cell 0
  TEST_FLOATING_EQUALITY(force[0], ref_soln_x, 1.0e-11);
  TEST_FLOATING_EQUALITY(force[1], ref_soln_y, 1.0e-11);
  TEST_FLOATING_EQUALITY(force[2], ref_soln_z, 1.0e-11);

  // cell 1
  // force on cell 1 due to interaction with cell 0
  // t_1 < x_0 - x_1 > dV_0
  ref_soln_x = 3398655528.6806289*cellVolume[0];
  ref_soln_y = 1863778838.308732*cellVolume[0];
  ref_soln_z = 767438345.18594845*cellVolume[0];
  // t_0 < x_1 - x_0 > dV_0
  ref_soln_x -= -2753550531.144094*cellVolume[0];
  ref_soln_y -= -1510011581.595148*cellVolume[0];
  ref_soln_z -= -621769474.7744727*cellVolume[0];
  // add force on cell 1 due to interaction with cell 2
  // t_1 < x_2 - x_1 >
  ref_soln_x += 5110873051.7924152*cellVolume[2];
  ref_soln_y += -1252163897.689138*cellVolume[2];
  ref_soln_z += 6133047662.15088*cellVolume[2];
  // - t_2 < x_1 - x_2 >
  ref_soln_x -= -4649613866.523322*cellVolume[2];
  ref_soln_y -= 1139155397.2982139*cellVolume[2];
  ref_soln_z -= -5579536639.827986*cellVolume[2];
  // assert the values of net force on cell 1
  TEST_FLOATING_EQUALITY(force[3], ref_soln_x, 1.0e-11);
  TEST_FLOATING_EQUALITY(force[4], ref_soln_y, 1.0e-10);
  TEST_FLOATING_EQUALITY(force[5], ref_soln_z, 1.0e-11);

  // cell 2
  // add force on cell 2 due to interaction with cell 0
  // t_2 < x_0 - x_2 > dV_2
  ref_soln_x = -8687228655.850972*cellVolume[0];
  ref_soln_y = -17295482505.73966*cellVolume[0];
  ref_soln_z = 13425717013.587864*cellVolume[0];
  // - t_0 < x_2 - x_0 > dV_2
  ref_soln_x -= 7736514797.0571524*cellVolume[0];
  ref_soln_y -= 15402697641.41379*cellVolume[0];
  ref_soln_z -= -11956431959.088315*cellVolume[0];
  // t_2 < x_1 - x_2 >
  ref_soln_x += -4649613866.523322*cellVolume[1];
  ref_soln_y += 1139155397.2982139*cellVolume[1];
  ref_soln_z += -5579536639.827986*cellVolume[1];
  // - t_1 < x_2 - x_1 >
  ref_soln_x -= 5110873051.7924152*cellVolume[1];
  ref_soln_y -= -1252163897.689138*cellVolume[1];
  ref_soln_z -= 6133047662.15088*cellVolume[1];
  // assert the values of net force on cell 2
  TEST_FLOATING_EQUALITY(force[6], ref_soln_x, 1.0e-11);
  TEST_FLOATING_EQUALITY(force[7], ref_soln_y, 1.0e-11);
  TEST_FLOATING_EQUALITY(force[8], ref_soln_z, 1.0e-11);

  delete[] ownedIDs;
  delete[] neighborhoodList;
}

//! Tests the automatic-differentiation Jacobian for a two-point system.
TEUCHOS_UNIT_TEST(MultiphysicsElasticMaterial, twoPointADJacobian) {

  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("Horizon", 10.0);
	params.set("Apply Automatic Differentiation Jacobian", true);
  params.set("Permeability", 1.e-12);
	params.set("Fluid density", 1000.0);
	params.set("Fluid dynamic viscosity", 1.0);
	params.set("Fluid compressibility", 1.e-4);
	params.set("Fluid Reynolds viscosity temperature effect", 0.0);
	params.set("Fluid linear thermal expansion", 0.0);
	params.set("Permeability curve inflection damage", .50);
	params.set("Permeability alpha", .25);
	params.set("Max permeability", 1.e-4);

  MultiphysicsElasticMaterial mat(params);
	ElasticMaterial solidMat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
	/*  
  Epetra_Map nodeMap(2, 0, comm);
  Epetra_Map unknownMap(6, 0, comm);
  Epetra_Map bondMap(1, 0, comm); 
	*/
	// Maps are needed to form the jacobian matrices which are incompatible with the
	// blockmaps needed by the data manager.
	Epetra_Map jacobianRowMap(8, 0, comm);
	Epetra_Map jacobianRowMapSolids(6, 0, comm);

  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;
	double dt = 1.0;
  // \todo check field specs

  // set up discretization
  numOwnedPoints = 2;
  ownedIDs = new int[numOwnedPoints];
	for(int i=0; i<numOwnedPoints; ++i){
		ownedIDs[i] = i;
	}
	  
  neighborhoodList = new int[4];
  neighborhoodList[0] = 1;
  neighborhoodList[1] = 1;
  neighborhoodList[2] = 1;
  neighborhoodList[3] = 0;

	double rowOEight[8];
	int colIndices[8];
	for(int col=0; col<8; ++col){
		rowOEight[col] = 0.0;
		colIndices[col] = col;
	}

  Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(2, 2, &ownedIDs[0], 1, 0, comm));
  Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(2, 2, &ownedIDs[0], 3, 0, comm));
  Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(2, 2, &ownedIDs[0], 1, 0, comm));

  // create the data manager
  // in serial, the overlap and non-overlap maps are the same
  PeridigmNS::DataManager dataManager;
  PeridigmNS::DataManager dataManagerSolids;
	dataManager.setMaps(tempOneDimensionalMap,
											tempOneDimensionalMap,
											tempThreeDimensionalMap,
											tempThreeDimensionalMap,
											tempBondMap);
	dataManagerSolids.setMaps(tempOneDimensionalMap,
														tempOneDimensionalMap,
														tempThreeDimensionalMap,
														tempThreeDimensionalMap,
														tempBondMap);

  dataManager.allocateData(mat.FieldIds());
  dataManagerSolids.allocateData(solidMat.FieldIds());

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  PeridigmNS::FieldManager& fieldManagerSolids = PeridigmNS::FieldManager::self();
  int modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  int fluidPressureYFieldId = fieldManager.getFieldId("Fluid_Pressure_Y");
  int volumeFieldId = fieldManager.getFieldId("Volume");
//   int weightedVolumeFieldId = fieldManager.getFieldId("Weighted_Volume");
//   int dilatationFieldId = fieldManager.getFieldId("Dilatation");
//   int bondDamageFieldId = fieldManager.getFieldId("Bond_Damage");
//   int forceDensityFieldId = fieldManager.getFieldId("Force_Density");
//   int flowDensityFieldId = fieldManager.getFieldId("Flux_Density");

  int modelCoordinatesFieldIdSolids = fieldManagerSolids.getFieldId("Model_Coordinates");
  int coordinatesFieldIdSolids = fieldManagerSolids.getFieldId("Coordinates");
  int volumeFieldIdSolids = fieldManagerSolids.getFieldId("Volume");
//   int weightedVolumeFieldIdSolids = fieldManagerSolids.getFieldId("Weighted_Volume");
//   int dilatationFieldIdSolids = fieldManagerSolids.getFieldId("Dilatation");
//   int bondDamageFieldIdSolids = fieldManagerSolids.getFieldId("Bond_Damage");
//   int forceDensityFieldIdSolids = fieldManagerSolids.getFieldId("Force_Density");

  // create the Jacobian

	Teuchos::RCP<Epetra_FECrsMatrix> jacobianFEC = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, jacobianRowMap /* rowmap */, 8 /*entries per row*/));
	for(int row=0; row <8; ++row){
		jacobianFEC->InsertGlobalValues(row, 
																		8 /* num entries */,
																		rowOEight /* values */,
																		colIndices /* column indices */);
	}
	jacobianFEC->GlobalAssemble();
  PeridigmNS::SerialMatrix jacobian(jacobianFEC);

	Teuchos::RCP<Epetra_FECrsMatrix> jacobianSolidsFEC = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, jacobianRowMapSolids/* rowmap */, 6 /*entries per row*/));
	for(int row=0; row <6; ++row){
		jacobianSolidsFEC->InsertGlobalValues(row, 
																		6 /* num entries */,
																		rowOEight /* values */,
																		colIndices /* column indices */);
	}
	jacobianSolidsFEC->GlobalAssemble();
  PeridigmNS::SerialMatrix jacobianSolids(jacobianSolidsFEC);

  Epetra_Vector& x = *dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE);
  Epetra_Vector& y = *dataManager.getData(coordinatesFieldId, PeridigmField::STEP_NP1);
	Epetra_Vector& fluidPressureY = *dataManager.getData(fluidPressureYFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& cellVolume = *dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE);

	Epetra_Vector& xSolids = *dataManagerSolids.getData(modelCoordinatesFieldIdSolids, PeridigmField::STEP_NONE);
  Epetra_Vector& ySolids = *dataManagerSolids.getData(coordinatesFieldIdSolids, PeridigmField::STEP_NP1);
  Epetra_Vector& cellVolumeSolids = *dataManagerSolids.getData(volumeFieldIdSolids, PeridigmField::STEP_NONE);

  x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
  x[3] = 1.0; x[4] = 0.0; x[5] = 0.0;
  y[0] = 0.0; y[1] = 0.0; y[2] = 0.0;
  y[3] = 2.0; y[4] = 0.0; y[5] = 0.0;

	fluidPressureY[0] = 0.0;
	fluidPressureY[1] = 0.0;

	xSolids[0] = 0.0; xSolids[1] = 0.0; xSolids[2] = 0.0;
  xSolids[3] = 1.0; xSolids[4] = 0.0; xSolids[5] = 0.0;
  ySolids[0] = 0.0; ySolids[1] = 0.0; ySolids[2] = 0.0;
  ySolids[3] = 2.0; ySolids[4] = 0.0; ySolids[5] = 0.0;

  for(int i=0; i<cellVolume.MyLength(); ++i){
 		cellVolume[i] = 1.0;
		cellVolumeSolids[i] = 1.0;
  }

	PeridigmNS::Material::JacobianType jacotype = PeridigmNS::Material::FULL_MATRIX;

  mat.initialize(dt, 
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 dataManager);
	  
  mat.computeJacobian(dt, 
                      numOwnedPoints,
                      ownedIDs,
                      neighborhoodList,
                      dataManager,
                      jacobian,
											jacotype);

  solidMat.initialize(dt, 
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 dataManagerSolids);

  solidMat.computeJacobian(dt, 
                      numOwnedPoints,
                      ownedIDs,
                      neighborhoodList,
                      dataManagerSolids,
                      jacobianSolids,
											jacotype);
  
	std::cout << "For the multiphysics material, the jacobian: " << std::endl;
	jacobian.getFECrsMatrix()->Print(std::cout);
	std::cout << "For the solids only material, the jacobian: " << std::endl;
	jacobianSolids.getFECrsMatrix()->Print(std::cout);

}

//! Tests the automatic-differentiation Jacobian for a two-point system.
/*  
TEUCHOS_UNIT_TEST(MultiphysicsElasticMaterial, twoPointProbeJacobianJAM) {

  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);  // \todo Fix these units, if needed.
  params.set("Bulk Modulus", 130000.0);
  params.set("Shear Modulus", 78000.0);
  params.set("Horizon", 10.0);
	params.set("Apply Automatic Differentiation Jacobian", true);
  params.set("Permeability", 1.0);
	params.set("Fluid density", 1000.0);
	params.set("Fluid dynamic viscosity", 1.0);
	params.set("Fluid compressibility", 1.0);
	params.set("Fluid Reynolds viscosity temperature effect", 0.0);
	params.set("Fluid linear thermal expansion", 0.0);
	params.set("Permeability curve inflection damage", .50);

  MultiphysicsElasticMaterial mat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(2, 0, comm);
  Epetra_Map unknownMap(6, 0, comm);
  Epetra_Map bondMap(2, 0, comm);
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;
  // \todo check field specs

  // set up discretization
  numOwnedPoints = 2;
  ownedIDs = new int[numOwnedPoints];
  ownedIDs[0] = 0;
  ownedIDs[1] = 1;
  neighborhoodList = new int[4];
  neighborhoodList[0] = 1;
  neighborhoodList[1] = 1;
  neighborhoodList[2] = 1;
  neighborhoodList[3] = 0;

  // create the data manager
  // in serial, the overlap and non-overlap maps are the same
  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&bondMap, false));
  dataManager.allocateData(mat.FieldIds());

  // create the Jacobian
//   PeridigmNS::SerialMatrix jacobian(2*3);

//   Epetra_Vector& x = *dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE);
//   Epetra_Vector& y = *dataManager.getData(coordinatesFieldId, PeridigmField::STEP_NP1);
//   Epetra_Vector& cellVolume = *dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE);

//   x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
//   x[3] = 0.5; x[4] = 0.0; x[5] = 0.0;
//   y[0] = 0.0; y[1] = 0.0; y[2] = 0.0;
//   y[3] = 0.5 - 0.058327; y[4] = 0.255; y[5] = 0.0;
//   for(int i=0; i<cellVolume.MyLength(); ++i){
// 	cellVolume[i] = 0.5;
//   }

//   mat.initialize(dt, 
//                  numOwnedPoints,
//                  ownedIDs,
//                  neighborhoodList,
//                  dataManager);

//   mat.computeJacobian(dt, 
//                       numOwnedPoints,
//                       ownedIDs,
//                       neighborhoodList,
//                       dataManager,
//                       jacobian);

//   cout << "\nJacobian for twoPointProbeJacobianJAM:" << endl;
//   jacobian.print(cout);
}
*/

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}

