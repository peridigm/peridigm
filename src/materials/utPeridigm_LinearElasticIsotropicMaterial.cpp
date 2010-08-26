/*! \file utPeridigm_LinearElasticIsotropicMaterial.cpp */

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
#include "Peridigm_LinearElasticIsotropicMaterial.hpp"
#include <Epetra_SerialComm.h>

using namespace boost::unit_test;
using namespace std;
using namespace PeridigmNS;
using namespace Teuchos;

//! Tests state variable count and name accessor functions.
void testStateVariableAccessors()
{
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  LinearElasticIsotropicMaterial mat(params);
  BOOST_REQUIRE(mat.NumScalarConstitutiveVariables() == 3);
  BOOST_CHECK(mat.ScalarConstitutiveVariableName(0) == "Weighted Volume");
  BOOST_CHECK(mat.ScalarConstitutiveVariableName(1) == "Dilatation");
  BOOST_CHECK(mat.ScalarConstitutiveVariableName(2) == "Damage");
  BOOST_CHECK_THROW(mat.ScalarConstitutiveVariableName(-1), std::range_error);
  BOOST_CHECK_THROW(mat.ScalarConstitutiveVariableName(3), std::range_error);
  BOOST_REQUIRE(mat.NumVectorConstitutiveVariables() == 1);
  BOOST_CHECK(mat.VectorConstitutiveVariableName(0) == "Current Position");
  BOOST_CHECK_THROW(mat.VectorConstitutiveVariableName(-1), std::range_error);
  BOOST_CHECK_THROW(mat.VectorConstitutiveVariableName(1), std::range_error);
}

//! Tests two-point system under compression against hand calculations.
void testTwoPts()
{
  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  LinearElasticIsotropicMaterial mat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(2, 0, comm);
  Epetra_Map unknownMap(6, 0, comm);
  Epetra_Vector x(unknownMap);
  Epetra_Vector u(unknownMap);
  Epetra_Vector v(unknownMap);
  double dt = 1.0;
  Epetra_Vector cellVolume(nodeMap);
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;
  double* bondState;
  int numScalarStateVariables = mat.NumScalarConstitutiveVariables();
  BOOST_CHECK(mat.NumScalarConstitutiveVariables() == 3);
  Epetra_MultiVector scalarConstitutiveData(nodeMap, numScalarStateVariables);
  int numVectorStateVariables = mat.NumVectorConstitutiveVariables();
  BOOST_CHECK(mat.NumVectorConstitutiveVariables() == 1);
  Epetra_MultiVector vectorConstitutiveData(unknownMap, numVectorStateVariables);
  Epetra_Vector force(unknownMap);

  x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
  x[3] = 1.0; x[4] = 0.0; x[5] = 0.0;
  u[0] = 0.0; u[1] = 0.0; u[2] = 0.0;
  u[3] = 0.0; u[4] = 0.0; u[5] = 0.0;
  v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
  v[3] = 1.0; v[4] = 0.0; v[5] = 0.0;
  for(int i=0; i<cellVolume.MyLength(); ++i){
	cellVolume[i] = 1.0;
  }
  scalarConstitutiveData.PutScalar(0.0);
  vectorConstitutiveData.PutScalar(0.0);

  numOwnedPoints = 2;
  ownedIDs = new int[numOwnedPoints];
  ownedIDs[0] = 0;
  ownedIDs[1] = 1;
  neighborhoodList = new int[4];
  neighborhoodList[0] = 1;
  neighborhoodList[1] = 1;
  neighborhoodList[2] = 1;
  neighborhoodList[3] = 0;
  int numBonds = 2;
  bondState = new double[numBonds];
  bondState[0] = 0.0;
  bondState[1] = 0.0;
  // bondMap
  // used for storing constitutive data on bonds
  int numGlobalElements = -1;
  int numMyElements = numBonds;
  int indexBase = 0;
  Epetra_Map bondMap(numGlobalElements, numMyElements, indexBase, comm);
  int numBondConstitutiveVariables = mat.NumBondConstitutiveVariables();
  BOOST_CHECK(mat.NumBondConstitutiveVariables() == 0);
  if(numBondConstitutiveVariables < 1)
    numBondConstitutiveVariables = 1;
  Epetra_MultiVector bondConstitutiveData(bondMap, numBondConstitutiveVariables);

  mat.initialize(x, 
                 u, 
                 v, 
                 dt, 
                 cellVolume,
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 bondState,
                 scalarConstitutiveData,
                 vectorConstitutiveData,
                 bondConstitutiveData,
                 force);

  mat.updateConstitutiveData(x, 
							 u, 
							 v, 
							 dt, 
							 cellVolume,
							 numOwnedPoints,
							 ownedIDs,
							 neighborhoodList,
							 bondState,
							 scalarConstitutiveData,
							 vectorConstitutiveData,
                             bondConstitutiveData,
							 force);

  double currentPositionX1 = vectorConstitutiveData[0][0];
  BOOST_CHECK_SMALL(currentPositionX1, 1.0e-14);
  double currentPositionY1 = vectorConstitutiveData[0][1];
  BOOST_CHECK_SMALL(currentPositionY1, 1.0e-14);
  double currentPositionZ1 = vectorConstitutiveData[0][2];
  BOOST_CHECK_SMALL(currentPositionZ1, 1.0e-14);
  double currentPositionX2 = vectorConstitutiveData[0][3];
  BOOST_CHECK_CLOSE(currentPositionX2, 2.0, 1.0e-12);
  double currentPositionY2 = vectorConstitutiveData[0][4];
  BOOST_CHECK_SMALL(currentPositionY2, 1.0e-14);
  double currentPositionZ2 = vectorConstitutiveData[0][5];
  BOOST_CHECK_SMALL(currentPositionZ2, 1.0e-14);
  double weightedVolume = scalarConstitutiveData[0][0];
  BOOST_CHECK_CLOSE(weightedVolume, 1.0, 1.0e-15);
  weightedVolume = scalarConstitutiveData[0][1];
  BOOST_CHECK_CLOSE(weightedVolume, 1.0, 1.0e-15);
  double dilatation = scalarConstitutiveData[1][0];
  BOOST_CHECK_CLOSE(dilatation, 3.0, 1.0e-15);
  dilatation = scalarConstitutiveData[1][1];
  BOOST_CHECK_CLOSE(dilatation, 3.0, 1.0e-15);
  double bondDatum = bondConstitutiveData[0][0];
  BOOST_CHECK_SMALL(bondDatum, 1.0e-15);
  bondDatum = bondConstitutiveData[0][1];
  BOOST_CHECK_SMALL(bondDatum, 1.0e-15);

  mat.computeForce(x, 
				   u, 
				   v, 
				   dt, 
				   cellVolume,
				   numOwnedPoints,
				   ownedIDs,
				   neighborhoodList,
				   bondState,
				   scalarConstitutiveData,
				   vectorConstitutiveData,
                   bondConstitutiveData,
				   force);

  BOOST_CHECK_CLOSE(force[0], 2.34e+12, 1.0e-2);
  BOOST_CHECK_SMALL(force[1], 1.0e-14);
  BOOST_CHECK_SMALL(force[2], 1.0e-14);
  BOOST_CHECK_CLOSE(force[3], -2.34e+12, 1.0e-2);
  BOOST_CHECK_SMALL(force[4], 1.0e-14);
  BOOST_CHECK_SMALL(force[5], 1.0e-14);

  delete[] ownedIDs;
  delete[] neighborhoodList;
  delete[] bondState;
}

//! Tests eight-cell block under compression against hand calculations.
void testEightPts()
{
  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  LinearElasticIsotropicMaterial mat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(8, 0, comm);
  Epetra_Map unknownMap(24, 0, comm);
  Epetra_Vector x(unknownMap);
  Epetra_Vector u(unknownMap);
  Epetra_Vector v(unknownMap);
  double dt = 1.0;
  Epetra_Vector cellVolume(nodeMap);
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;
  double* bondState;
  int numScalarStateVariables = mat.NumScalarConstitutiveVariables();
  BOOST_CHECK(mat.NumScalarConstitutiveVariables() == 3);
  Epetra_MultiVector scalarConstitutiveData(nodeMap, numScalarStateVariables);
  int numVectorStateVariables = mat.NumVectorConstitutiveVariables();
  BOOST_CHECK(mat.NumVectorConstitutiveVariables() == 1);
  Epetra_MultiVector vectorConstitutiveData(unknownMap, numVectorStateVariables);
  Epetra_Vector force(unknownMap);

  // initial positions
  x[0]  = 0.0; x[1]  = 0.0; x[2]  = 0.0;
  x[3]  = 0.0; x[4]  = 0.0; x[5]  = 1.0;
  x[6]  = 0.0; x[7]  = 1.0; x[8]  = 0.0;
  x[9]  = 0.0; x[10] = 1.0; x[11] = 1.0;
  x[12] = 1.0; x[13] = 0.0; x[14] = 0.0;
  x[15] = 1.0; x[16] = 0.0; x[17] = 1.0;
  x[18] = 1.0; x[19] = 1.0; x[20] = 0.0;
  x[21] = 1.0; x[22] = 1.0; x[23] = 1.0;

  // displacements (compression in z direction)
  u[0]  = 0.0; u[1]  = 0.0; u[2]  = 0.0;
  u[3]  = 0.0; u[4]  = 0.0; u[5]  = -0.02;
  u[6]  = 0.0; u[7]  = 0.0; u[8]  = 0.0;
  u[9]  = 0.0; u[10] = 0.0; u[11] = -0.02;
  u[12] = 0.0; u[13] = 0.0; u[14] = 0.0;
  u[15] = 0.0; u[16] = 0.0; u[17] = -0.02;
  u[18] = 0.0; u[19] = 0.0; u[20] = 0.0;
  u[21] = 0.0; u[22] = 0.0; u[23] = -0.02;

  // velocities (not used by material model)
  v[0]  = 0.0; v[1]  = 0.0; v[2]  = 0.0;
  v[3]  = 0.0; v[4]  = 0.0; v[5]  = 0.0;
  v[6]  = 0.0; v[7]  = 0.0; v[8]  = 0.0;
  v[9]  = 0.0; v[10] = 0.0; v[11] = 0.0;
  v[12] = 0.0; v[13] = 0.0; v[14] = 0.0;
  v[15] = 0.0; v[16] = 0.0; v[17] = 0.0;
  v[18] = 0.0; v[19] = 0.0; v[20] = 0.0;
  v[21] = 0.0; v[22] = 0.0; v[23] = 0.0;

  // cell volumes
  for(int i=0; i<cellVolume.MyLength(); ++i){
	cellVolume[i] = 1.0;
  }

  // zero out constitutive data
  scalarConstitutiveData.PutScalar(0.0);
  vectorConstitutiveData.PutScalar(0.0);

  // all cellss are neighbors of each other
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
  // total number of bonds = 8(7) = 56
  int numBonds = 56;
  bondState = new double[numBonds];
  for(int i=0 ; i<numBonds ; ++i){
	bondState[0] = 0.0;
  }
  // bondMap
  // used for storing constitutive data on bonds
  int numGlobalElements = -1;
  int numMyElements = numBonds;
  int indexBase = 0;
  Epetra_Map bondMap(numGlobalElements, numMyElements, indexBase, comm);
  int numBondConstitutiveVariables = mat.NumBondConstitutiveVariables();
  BOOST_CHECK(mat.NumBondConstitutiveVariables() == 0);
  if(numBondConstitutiveVariables < 1)
    numBondConstitutiveVariables = 1;
  Epetra_MultiVector bondConstitutiveData(bondMap, numBondConstitutiveVariables);

  mat.initialize(x, 
                 u, 
                 v, 
                 dt, 
                 cellVolume,
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 bondState,
                 scalarConstitutiveData,
                 vectorConstitutiveData,
                 bondConstitutiveData,
                 force);

  mat.updateConstitutiveData(x, 
							 u, 
							 v, 
							 dt, 
							 cellVolume,
							 numOwnedPoints,
							 ownedIDs,
							 neighborhoodList,
							 bondState,
							 scalarConstitutiveData,
							 vectorConstitutiveData,
                             bondConstitutiveData,
							 force);

  double currentPosition;
  currentPosition = vectorConstitutiveData[0][0];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][1];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][2];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][3];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][4];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][5];
  BOOST_CHECK_CLOSE(currentPosition, 0.98, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][6];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][7];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][8];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][9];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][10];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][11];
  BOOST_CHECK_CLOSE(currentPosition, 0.98, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][12];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][13];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][14];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][15];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][16];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][17];
  BOOST_CHECK_CLOSE(currentPosition, 0.98, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][18];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][19];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][20];
  BOOST_CHECK_SMALL(currentPosition, 1.0e-14);
  currentPosition = vectorConstitutiveData[0][21];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][22];
  BOOST_CHECK_CLOSE(currentPosition, 1.0, 1.0e-15);
  currentPosition = vectorConstitutiveData[0][23];
  BOOST_CHECK_CLOSE(currentPosition, 0.98, 1.0e-15);

  // the weighted volumes and dilatations are the
  // same for all points in this test problem
  for(int i=0; i<8; ++i){
	double weightedVolume = scalarConstitutiveData[0][i];
	double dilatation = scalarConstitutiveData[1][i];
	BOOST_CHECK_CLOSE(weightedVolume, 12.0, 1.0e-12);
	BOOST_CHECK_CLOSE(dilatation, -0.01991593994643333, 1.0e-12);
  }

  // the bond constitutive data should be all zeros (not
  // used by this material model)
  for(int i=0 ; i<numBonds ; ++i){
      BOOST_CHECK_SMALL(bondConstitutiveData[0][i], 1.0e-15);
  }

  mat.computeForce(x, 
				   u, 
				   v, 
				   dt, 
				   cellVolume,
				   numOwnedPoints,
				   ownedIDs,
				   neighborhoodList,
				   bondState,
				   scalarConstitutiveData,
				   vectorConstitutiveData,
                   bondConstitutiveData,
				   force);

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

  // assert the net force on cell 0
  BOOST_CHECK_CLOSE(force[0], f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(force[1], f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(force[2], f[2], 1.0e-11);

  // assert the net forces on the other cells
  // these forces differ only by direction

  // force on cell 1
  BOOST_CHECK_CLOSE(force[3],       f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(force[4],       f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[5],  f[2], 1.0e-11);

  // force on cell 2
  BOOST_CHECK_CLOSE(force[6],       f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[7],  f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(force[8],       f[2], 1.0e-11);

  // force on cell 3
  BOOST_CHECK_CLOSE(force[9],       f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[10], f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[11], f[2], 1.0e-11);

  // force on cell 4
  BOOST_CHECK_CLOSE(-1.0*force[12], f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(force[13],      f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(force[14],      f[2], 1.0e-11);

  // force on cell 5
  BOOST_CHECK_CLOSE(-1.0*force[15], f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(force[16],      f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[17], f[2], 1.0e-11);

  // force on cell 6
  BOOST_CHECK_CLOSE(-1.0*force[18], f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[19], f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(force[20],      f[2], 1.0e-11);

  // force on cell 7
  BOOST_CHECK_CLOSE(-1.0*force[21], f[0], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[22], f[1], 1.0e-11);
  BOOST_CHECK_CLOSE(-1.0*force[23], f[2], 1.0e-11);

  delete[] ownedIDs;
  delete[] neighborhoodList;
  delete[] bondState;
}

//! Tests arbitrary three-cell system under arbitrary deformation against hand calculations.
void testThreePts()
{
  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  LinearElasticIsotropicMaterial mat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(3, 0, comm);
  Epetra_Map unknownMap(9, 0, comm);
  Epetra_Vector x(unknownMap);
  Epetra_Vector u(unknownMap);
  Epetra_Vector v(unknownMap);
  double dt = 1.0;
  Epetra_Vector cellVolume(nodeMap);
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;
  double* bondState;
  int numScalarStateVariables = mat.NumScalarConstitutiveVariables();
  BOOST_CHECK(mat.NumScalarConstitutiveVariables() == 3);
  Epetra_MultiVector scalarConstitutiveData(nodeMap, numScalarStateVariables);
  int numVectorStateVariables = mat.NumVectorConstitutiveVariables();
  BOOST_CHECK(mat.NumVectorConstitutiveVariables() == 1);
  Epetra_MultiVector vectorConstitutiveData(unknownMap, numVectorStateVariables);
  Epetra_Vector force(unknownMap);

  // initial positions
  x[0] =  1.1; x[1] = 2.6;  x[2] = -0.1;
  x[3] = -2.0; x[4] = 0.9;  x[5] = -0.3;
  x[6] =  0.0; x[7] = 0.01; x[8] =  1.8;

  // displacements
  u[0] = 0.1; u[1] = -0.2; u[2] =  0.0;
  u[3] = 0.1; u[4] = -0.2; u[5] = -0.5;
  u[6] = 0.1; u[7] =  0.2; u[8] = -0.2;

  // velocities (not used by material model)
  v[0]  = 0.0; v[1]  = 0.0; v[2]  = 0.0;
  v[3]  = 0.0; v[4]  = 0.0; v[5]  = 0.0;
  v[6]  = 0.0; v[7]  = 0.0; v[8]  = 0.0;

  // cell volumes
  cellVolume[0] = 0.9;
  cellVolume[1] = 1.1;
  cellVolume[2] = 0.8;

  // zero out constitutive data
  scalarConstitutiveData.PutScalar(0.0);
  vectorConstitutiveData.PutScalar(0.0);

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
  // total number of bonds = 3(2) = 6
  int numBonds = 6;
  bondState = new double[numBonds];
  for(int i=0 ; i<numBonds ; ++i){
	bondState[0] = 0.0;
  }
  // bondMap
  // used for storing constitutive data on bonds
  int numGlobalElements = -1;
  int numMyElements = numBonds;
  int indexBase = 0;
  Epetra_Map bondMap(numGlobalElements, numMyElements, indexBase, comm);
  int numBondConstitutiveVariables = mat.NumBondConstitutiveVariables();
  BOOST_CHECK(mat.NumBondConstitutiveVariables() == 0);
  if(numBondConstitutiveVariables < 1)
    numBondConstitutiveVariables = 1;
  Epetra_MultiVector bondConstitutiveData(bondMap, numBondConstitutiveVariables);

  mat.initialize(x, 
                 u, 
                 v, 
                 dt, 
                 cellVolume,
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 bondState,
                 scalarConstitutiveData,
                 vectorConstitutiveData,
                 bondConstitutiveData,
                 force);

  mat.updateConstitutiveData(x, 
							 u, 
							 v, 
							 dt, 
							 cellVolume,
							 numOwnedPoints,
							 ownedIDs,
							 neighborhoodList,
							 bondState,
							 scalarConstitutiveData,
							 vectorConstitutiveData,
                             bondConstitutiveData,
							 force);

  double currentPosition;
  currentPosition = vectorConstitutiveData[0][0];
  BOOST_CHECK_CLOSE(currentPosition, 1.2, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][1];
  BOOST_CHECK_CLOSE(currentPosition, 2.4, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][2];
  BOOST_CHECK_CLOSE(currentPosition, -0.1, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][3];
  BOOST_CHECK_CLOSE(currentPosition, -1.9, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][4];
  BOOST_CHECK_CLOSE(currentPosition, 0.7, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][5];
  BOOST_CHECK_CLOSE(currentPosition, -0.8, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][6];
  BOOST_CHECK_CLOSE(currentPosition, 0.1, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][7];
  BOOST_CHECK_CLOSE(currentPosition, 0.21, 1.0e-13);
  currentPosition = vectorConstitutiveData[0][8];
  BOOST_CHECK_CLOSE(currentPosition, 1.6, 1.0e-13);

  // check the weighted volume and dilatation
  // against hand calculations
  double weightedVolume, dilatation;
  weightedVolume = scalarConstitutiveData[0][0];
  BOOST_CHECK_CLOSE(weightedVolume, 23.016479999999931, 1.0e-12);
  dilatation = scalarConstitutiveData[1][0];
  BOOST_CHECK_CLOSE(dilatation, -0.114127034572639, 1.0e-11);
  weightedVolume = scalarConstitutiveData[0][1];
  BOOST_CHECK_CLOSE(weightedVolume,18.64767999999997, 1.0e-12);
  dilatation = scalarConstitutiveData[1][1];
  BOOST_CHECK_CLOSE(dilatation, 0.08257537372985, 1.0e-11);
  weightedVolume = scalarConstitutiveData[0][2];
  BOOST_CHECK_CLOSE(weightedVolume, 20.497599999999963, 1.0e-12);
  dilatation = scalarConstitutiveData[1][2];
  BOOST_CHECK_CLOSE(dilatation, -0.12166177890553, 1.0e-11);

  // the bond constitutive data should be all zeros (not
  // used by this material model)
  for(int i=0 ; i<numBonds ; ++i){
      BOOST_CHECK_SMALL(bondConstitutiveData[0][i], 1.0e-15);
  }

  mat.computeForce(x, 
				   u, 
				   v, 
				   dt, 
				   cellVolume,
				   numOwnedPoints,
				   ownedIDs,
				   neighborhoodList,
				   bondState,
				   scalarConstitutiveData,
				   vectorConstitutiveData,
                   bondConstitutiveData,
				   force);

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
  BOOST_CHECK_CLOSE(force[0], ref_soln_x, 1.0e-11);
  BOOST_CHECK_CLOSE(force[1], ref_soln_y, 1.0e-11);
  BOOST_CHECK_CLOSE(force[2], ref_soln_z, 1.0e-11);

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
  BOOST_CHECK_CLOSE(force[3], ref_soln_x, 1.0e-11);
  BOOST_CHECK_CLOSE(force[4], ref_soln_y, 1.0e-10);
  BOOST_CHECK_CLOSE(force[5], ref_soln_z, 1.0e-11);

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
  BOOST_CHECK_CLOSE(force[6], ref_soln_x, 1.0e-11);
  BOOST_CHECK_CLOSE(force[7], ref_soln_y, 1.0e-11);
  BOOST_CHECK_CLOSE(force[8], ref_soln_z, 1.0e-11);

  delete[] ownedIDs;
  delete[] neighborhoodList;
  delete[] bondState;
}

bool init_unit_test_suite()
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPeridigm_LinearElasticIsotropicMaterial");
  proc->add(BOOST_TEST_CASE(&testStateVariableAccessors));
  proc->add(BOOST_TEST_CASE(&testTwoPts));
  proc->add(BOOST_TEST_CASE(&testEightPts));
  proc->add(BOOST_TEST_CASE(&testThreePts));
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
  // Initialize UTF
  return unit_test_main(init_unit_test, argc, argv);
}
