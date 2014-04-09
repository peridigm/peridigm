/*! \file utPeridigm_GeometryUtils.cpp */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2014) Sandia Corporation
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

#include "Peridigm_GeometryUtils.hpp"
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#include <Epetra_SerialComm.h>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif

using namespace Teuchos;
using namespace PeridigmNS;
using namespace std;

//! Exercise tetVolume() and tetCentroid functions
TEUCHOS_UNIT_TEST(GeometryUtils, TetGeometry) {

  // Volume of a tetrahedron is A h / 3, where A is the area of the base and h is the height
  // Centroid is at the intersection of lines from each node to the centroid of the opposite face (a quarter of the distance along one of these lines)

  vector<double*> nodes(4);
  vector<double> n1(3), n2(3), n3(3), n4(3);
  nodes[0] = &n1[0] ; nodes[1] = &n2[0] ; nodes[2] = &n3[0] ; nodes[3] = &n4[0] ;
  vector<double> u(3), v(3), cross(3), normal(3), centroid(3), baseCentroid(3), trueCentroid(3);
  double volume(0.0), baseArea(0.0), trueVolume(0.0), mag(0.0), height(0.0);
  double relTolerance = 1.0e-14;
  
  // Simple case
  nodes[0][0] = 0.0 ; nodes[0][1] = 0.0 ; nodes[0][2] = 0.0 ;
  nodes[1][0] = 1.0 ; nodes[1][1] = 0.0 ; nodes[1][2] = 0.0 ;
  nodes[2][0] = 0.0 ; nodes[2][1] = 1.0 ; nodes[2][2] = 0.0 ;
  nodes[3][0] = 0.0 ; nodes[3][1] = 0.0 ; nodes[3][2] = 1.0 ;
  baseArea = 0.5;
  trueVolume = baseArea * 1.0 / 3.0;
  volume = tetVolume(nodes);
  TEST_FLOATING_EQUALITY(volume, trueVolume, relTolerance);
  // Compute the centroid of the base
  baseCentroid[0] = 1.0/3.0; 
  baseCentroid[1] = 1.0/3.0; 
  baseCentroid[2] = 0.0;
  // Find the centroid of the tet
  trueCentroid[0] = 0.25*(nodes[3][0] - baseCentroid[0]) + baseCentroid[0];
  trueCentroid[1] = 0.25*(nodes[3][1] - baseCentroid[1]) + baseCentroid[1];
  trueCentroid[2] = 0.25*(nodes[3][2] - baseCentroid[2]) + baseCentroid[2];
  tetCentroid(nodes, centroid);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);

  // Slightly more complex case
  nodes[0][0] =  2.0 ; nodes[0][1] =  2.0 ; nodes[0][2] =  0.0 ;
  nodes[1][0] = -2.0 ; nodes[1][1] = -2.0 ; nodes[1][2] =  0.0 ;
  nodes[2][0] = -1.0 ; nodes[2][1] =  1.0 ; nodes[2][2] =  0.0 ;
  nodes[3][0] = -0.1 ; nodes[3][1] = 50.0 ; nodes[3][2] = -1.2 ;
  baseArea = 0.5 * 2.0 * std::sqrt(8.0) * std::sqrt(2.0);
  trueVolume = baseArea * 1.2 / 3.0;
  volume = tetVolume(nodes);
  TEST_FLOATING_EQUALITY(volume, trueVolume, relTolerance);
  // Compute the centroid of the base
  baseCentroid[0] = (nodes[0][0] + nodes[1][0] + nodes[2][0])/3.0; 
  baseCentroid[1] = (nodes[0][1] + nodes[1][1] + nodes[2][1])/3.0; 
  baseCentroid[2] = (nodes[0][2] + nodes[1][2] + nodes[2][2])/3.0;
  // Find the centroid of the tet
  trueCentroid[0] = 0.25*(nodes[3][0] - baseCentroid[0]) + baseCentroid[0];
  trueCentroid[1] = 0.25*(nodes[3][1] - baseCentroid[1]) + baseCentroid[1];
  trueCentroid[2] = 0.25*(nodes[3][2] - baseCentroid[2]) + baseCentroid[2];
  tetCentroid(nodes, centroid);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);

  // Made-up numbers for node locations
  nodes[0][0] =  12.1   ; nodes[0][1] = -0.01 ; nodes[0][2] =    0.0 ;
  nodes[1][0] =  1.1e5  ; nodes[1][1] = -2.0  ; nodes[1][2] =   10.0 ;
  nodes[2][0] = -2.3e-8 ; nodes[2][1] =  1.0  ; nodes[2][2] =    0.0 ;
  nodes[3][0] = -0.1    ; nodes[3][1] = -50.0 ; nodes[3][2] = -0.002 ;
  // Base area = 0.5 det( (x_3 - x_1) X (x_3 - x_2) )
  u[0] = nodes[2][0] - nodes[0][0];
  u[1] = nodes[2][1] - nodes[0][1];
  u[2] = nodes[2][2] - nodes[0][2];
  v[0] = nodes[2][0] - nodes[1][0];
  v[1] = nodes[2][1] - nodes[1][1];
  v[2] = nodes[2][2] - nodes[1][2];
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];
  mag = std::sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2] );
  baseArea = 0.5 * mag;
  // Find normal to base
  u[0] = nodes[1][0] - nodes[0][0];
  u[1] = nodes[1][1] - nodes[0][1];
  u[2] = nodes[1][2] - nodes[0][2];
  v[0] = nodes[2][0] - nodes[0][0];
  v[1] = nodes[2][1] - nodes[0][1];
  v[2] = nodes[2][2] - nodes[0][2];
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];
  mag = std::sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2] );
  normal[0] = cross[0] / mag;
  normal[1] = cross[1] / mag;
  normal[2] = cross[2] / mag;
  // Compute distance from plane to 4th node
  height = normal[0]*(nodes[3][0] - nodes[0][0]) + normal[1]*(nodes[3][1] - nodes[0][1]) + normal[2]*(nodes[3][2] - nodes[0][2]);
  trueVolume = baseArea * height / 3.0;
  volume = tetVolume(nodes);
  TEST_FLOATING_EQUALITY(volume, trueVolume, relTolerance);
  // Compute the centroid of the base
  baseCentroid[0] = (nodes[0][0] + nodes[1][0] + nodes[2][0])/3.0; 
  baseCentroid[1] = (nodes[0][1] + nodes[1][1] + nodes[2][1])/3.0; 
  baseCentroid[2] = (nodes[0][2] + nodes[1][2] + nodes[2][2])/3.0;
  // Find the centroid of the tet
  trueCentroid[0] = 0.25*(nodes[3][0] - baseCentroid[0]) + baseCentroid[0];
  trueCentroid[1] = 0.25*(nodes[3][1] - baseCentroid[1]) + baseCentroid[1];
  trueCentroid[2] = 0.25*(nodes[3][2] - baseCentroid[2]) + baseCentroid[2];
  tetCentroid(nodes, centroid);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);

  // How 'bout some random numbers from random.org
  nodes[0][0] = 273025420.0 ; nodes[0][1] = 589350625.0 ; nodes[0][2] = 607719774.0 ;
  nodes[1][0] = 183383954.0 ; nodes[1][1] = 816159710.0 ; nodes[1][2] = 897521393.0 ;
  nodes[2][0] = 113150325.0 ; nodes[2][1] = 418500907.0 ; nodes[2][2] =   4068851.0 ;
  nodes[3][0] = -24119980.0 ; nodes[3][1] = 783619877.0 ; nodes[3][2] = 783619877.0 ;
  // Base area = 0.5 det( (x_3 - x_1) X (x_3 - x_2) )
  u[0] = nodes[2][0] - nodes[0][0];
  u[1] = nodes[2][1] - nodes[0][1];
  u[2] = nodes[2][2] - nodes[0][2];
  v[0] = nodes[2][0] - nodes[1][0];
  v[1] = nodes[2][1] - nodes[1][1];
  v[2] = nodes[2][2] - nodes[1][2];
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];
  mag = std::sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2] );
  baseArea = 0.5 * mag;
  // Find normal to base
  u[0] = nodes[1][0] - nodes[0][0];
  u[1] = nodes[1][1] - nodes[0][1];
  u[2] = nodes[1][2] - nodes[0][2];
  v[0] = nodes[2][0] - nodes[0][0];
  v[1] = nodes[2][1] - nodes[0][1];
  v[2] = nodes[2][2] - nodes[0][2];
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];
  mag = std::sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2] );
  normal[0] = cross[0] / mag;
  normal[1] = cross[1] / mag;
  normal[2] = cross[2] / mag;
  // Compute distance from plane to 4th node
  height = normal[0]*(nodes[3][0] - nodes[0][0]) + normal[1]*(nodes[3][1] - nodes[0][1]) + normal[2]*(nodes[3][2] - nodes[0][2]);
  trueVolume = baseArea * height / 3.0;
  volume = tetVolume(nodes);
  TEST_FLOATING_EQUALITY(volume, trueVolume, relTolerance);
  // Compute the centroid of the base
  baseCentroid[0] = (nodes[0][0] + nodes[1][0] + nodes[2][0])/3.0; 
  baseCentroid[1] = (nodes[0][1] + nodes[1][1] + nodes[2][1])/3.0; 
  baseCentroid[2] = (nodes[0][2] + nodes[1][2] + nodes[2][2])/3.0;
  // Find the centroid of the tet
  trueCentroid[0] = 0.25*(nodes[3][0] - baseCentroid[0]) + baseCentroid[0];
  trueCentroid[1] = 0.25*(nodes[3][1] - baseCentroid[1]) + baseCentroid[1];
  trueCentroid[2] = 0.25*(nodes[3][2] - baseCentroid[2]) + baseCentroid[2];
  tetCentroid(nodes, centroid);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
  
  // And for the grand finale, put exponents on the random numbers
  nodes[0][0] = 273025420.0e5   ; nodes[0][1] = 589350625.0e3 ; nodes[0][2] =  607719774.0e3  ;
  nodes[1][0] = 183383954.0e5   ; nodes[1][1] = 816159710.0e2 ; nodes[1][2] =  897521393.0e-6 ;
  nodes[2][0] = 113150325.0e-10 ; nodes[2][1] = 418500907.0e3 ; nodes[2][2] =    4068851.0e-2 ;
  nodes[3][0] = -24119980.0e3   ; nodes[3][1] = 783619877.0e4 ; nodes[3][2] = -783619877.0e1  ;
  // Base area = 0.5 det( (x_3 - x_1) X (x_3 - x_2) )
  u[0] = nodes[2][0] - nodes[0][0];
  u[1] = nodes[2][1] - nodes[0][1];
  u[2] = nodes[2][2] - nodes[0][2];
  v[0] = nodes[2][0] - nodes[1][0];
  v[1] = nodes[2][1] - nodes[1][1];
  v[2] = nodes[2][2] - nodes[1][2];
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];
  mag = std::sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2] );
  baseArea = 0.5 * mag;
  // Find normal to base
  u[0] = nodes[1][0] - nodes[0][0];
  u[1] = nodes[1][1] - nodes[0][1];
  u[2] = nodes[1][2] - nodes[0][2];
  v[0] = nodes[2][0] - nodes[0][0];
  v[1] = nodes[2][1] - nodes[0][1];
  v[2] = nodes[2][2] - nodes[0][2];
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];
  mag = std::sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2] );
  normal[0] = cross[0] / mag;
  normal[1] = cross[1] / mag;
  normal[2] = cross[2] / mag;
  // Compute distance from plane to 4th node
  height = normal[0]*(nodes[3][0] - nodes[0][0]) + normal[1]*(nodes[3][1] - nodes[0][1]) + normal[2]*(nodes[3][2] - nodes[0][2]);
  trueVolume = baseArea * height / 3.0;
  volume = tetVolume(nodes);
  TEST_FLOATING_EQUALITY(volume, trueVolume, relTolerance);
  // Compute the centroid of the base
  baseCentroid[0] = (nodes[0][0] + nodes[1][0] + nodes[2][0])/3.0; 
  baseCentroid[1] = (nodes[0][1] + nodes[1][1] + nodes[2][1])/3.0; 
  baseCentroid[2] = (nodes[0][2] + nodes[1][2] + nodes[2][2])/3.0;
  // Find the centroid of the tet
  trueCentroid[0] = 0.25*(nodes[3][0] - baseCentroid[0]) + baseCentroid[0];
  trueCentroid[1] = 0.25*(nodes[3][1] - baseCentroid[1]) + baseCentroid[1];
  trueCentroid[2] = 0.25*(nodes[3][2] - baseCentroid[2]) + baseCentroid[2];
  tetCentroid(nodes, centroid);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
}

//! Exercise tetCentroidAndVolume()
TEUCHOS_UNIT_TEST(GeometryUtils, TetCentroidAndVolume) {

  vector<double> nodes(12), centroid(3), trueCentroid(3), volume(1);
  double trueVolume;
  double relTolerance = 1.0e-14;

  // Simple case
  nodes[0] = 0.0  ; nodes[1] = 0.0  ; nodes[2] = 0.0 ;
  nodes[3] = 1.0  ; nodes[4] = 0.0  ; nodes[5] = 0.0 ;
  nodes[6] = 0.0  ; nodes[7] = 1.0  ; nodes[8] = 0.0 ;
  nodes[9] = 0.0  ; nodes[10] = 0.0 ; nodes[11] = 1.0 ;
  trueCentroid[0] = 1.0/4.0;
  trueCentroid[1] = 1.0/4.0;
  trueCentroid[2] = 1.0/4.0;
  trueVolume = 1.0/6.0;
  tetCentroidAndVolume(&nodes[0], &centroid[0], &volume[0]);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
  TEST_FLOATING_EQUALITY(volume[0], trueVolume, relTolerance);

  // Another fairly simple case
  nodes[0] = 0.25  ; nodes[1] = 0.25  ; nodes[2] = 0.521 ;
  nodes[3] = -0.25 ; nodes[4] = 0.25  ; nodes[5] = 0.521 ;
  nodes[6] = 0.0   ; nodes[7] = -0.01 ; nodes[8] = 0.521 ;
  nodes[9] = 0.2   ; nodes[10] = 5.0  ; nodes[11] = 0.621;
  trueCentroid[0] = 0.2/4.0;
  trueCentroid[1] = (0.25 - 0.26/3.0) + 0.25*(5.0 - (0.25 - 0.26/3.0));
  trueCentroid[2] = 0.521 + 0.1/4.0;
  trueVolume = 0.0065/3.0;
  tetCentroidAndVolume(&nodes[0], &centroid[0], &volume[0]);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
  TEST_FLOATING_EQUALITY(volume[0], trueVolume, relTolerance);

  // Case with random numbers
  nodes[0] = 268.0  ; nodes[1] = -593.0 ; nodes[2] = 280.0  ;
  nodes[3] = -405.0 ; nodes[4] = -586.0 ; nodes[5] = 73.0   ;
  nodes[6] = -475.0 ; nodes[7] = -611.0 ; nodes[8] = 540.0  ;
  nodes[9] = -740.0 ; nodes[10] = 282.0 ; nodes[11] = -900.0;
  trueCentroid[0] = 0.25*(nodes[0] + nodes[3] + nodes[6] + nodes[9]);
  trueCentroid[1] = 0.25*(nodes[1] + nodes[4] + nodes[7] + nodes[10]);
  trueCentroid[2] = 0.25*(nodes[2] + nodes[5] + nodes[8] + nodes[11]);

  double DA[3] = {nodes[0]-nodes[9], nodes[1]-nodes[10], nodes[2]-nodes[11]};
  double DB[3] = {nodes[3]-nodes[9], nodes[4]-nodes[10], nodes[5]-nodes[11]};
  double DC[3] = {nodes[6]-nodes[9], nodes[7]-nodes[10], nodes[8]-nodes[11]};
  double DB_Cross_DC[3] = {DB[1]*DC[2] - DB[2]*DC[1], DB[2]*DC[0] - DB[0]*DC[2], DB[0]*DC[1] - DB[1]*DC[0]};
  double DA_Dot_DB_Cross_DC = DA[0]*DB_Cross_DC[0] + DA[1]*DB_Cross_DC[1] + DA[2]*DB_Cross_DC[2];
  trueVolume = std::abs(DA_Dot_DB_Cross_DC)/6.0;
  tetCentroidAndVolume(&nodes[0], &centroid[0], &volume[0]);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
  TEST_FLOATING_EQUALITY(volume[0], trueVolume, relTolerance);
}

//! Exercise hexCentroidAndVolume()
TEUCHOS_UNIT_TEST(GeometryUtils, HexCentroidAndVolume) {

  vector<double> nodes(24), centroid(3), trueCentroid(3), volume(1);
  double trueVolume;
  double relTolerance = 1.0e-14;

  // Simple case
  nodes[0] = 0.0  ; nodes[1] = 0.0  ; nodes[2] = 0.0 ;
  nodes[3] = 1.0  ; nodes[4] = 0.0  ; nodes[5] = 0.0 ;
  nodes[6] = 1.0  ; nodes[7] = 1.0  ; nodes[8] = 0.0 ;
  nodes[9] = 0.0  ; nodes[10] = 1.0 ; nodes[11] = 0.0 ;
  nodes[12] = 0.0 ; nodes[13] = 0.0 ; nodes[14] = 1.0 ;
  nodes[15] = 1.0 ; nodes[16] = 0.0 ; nodes[17] = 1.0 ;
  nodes[18] = 1.0 ; nodes[19] = 1.0 ; nodes[20] = 1.0 ;
  nodes[21] = 0.0 ; nodes[22] = 1.0 ; nodes[23] = 1.0 ;
  trueCentroid[0] = 0.5;
  trueCentroid[1] = 0.5;
  trueCentroid[2] = 0.5;
  trueVolume = 1.0;
  hexCentroidAndVolume(&nodes[0], &centroid[0], &volume[0]);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
  TEST_FLOATING_EQUALITY(volume[0], trueVolume, relTolerance);

  // Composite of a cube and a wedge
  nodes[0] = 0.0  ; nodes[1] = 0.0  ; nodes[2] = 0.0 ;
  nodes[3] = 1.0  ; nodes[4] = 0.0  ; nodes[5] = 0.0 ;
  nodes[6] = 1.0  ; nodes[7] = 1.0  ; nodes[8] = 0.0 ;
  nodes[9] = 0.0  ; nodes[10] = 1.0 ; nodes[11] = 0.0 ;
  nodes[12] = 0.0 ; nodes[13] = 0.0 ; nodes[14] = 2.0 ;
  nodes[15] = 1.0 ; nodes[16] = 0.0 ; nodes[17] = 1.0 ;
  nodes[18] = 1.0 ; nodes[19] = 1.0 ; nodes[20] = 1.0 ;
  nodes[21] = 0.0 ; nodes[22] = 1.0 ; nodes[23] = 2.0 ;
  trueCentroid[0] = (0.5*1.0/3.0 + 0.5)/1.5;
  trueCentroid[1] = (0.5*0.5     + 0.5)/1.5;
  trueCentroid[2] = (0.5*4.0/3.0 + 0.5)/1.5;
  trueVolume = 1.5;
  hexCentroidAndVolume(&nodes[0], &centroid[0], &volume[0]);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
  TEST_FLOATING_EQUALITY(volume[0], trueVolume, relTolerance);

  // Composite of a cube and two wedges
  nodes[0] =  0.0 ; nodes[1] = 0.0  ; nodes[2] = -1.0 ;
  nodes[3] =  1.0 ; nodes[4] = 0.0  ; nodes[5] =  0.0 ;
  nodes[6] =  1.0 ; nodes[7] = 1.0  ; nodes[8] =  0.0 ;
  nodes[9] =  0.0 ; nodes[10] = 1.0 ; nodes[11] = -1.0 ;
  nodes[12] = 0.0 ; nodes[13] = 0.0 ; nodes[14] = 2.0 ;
  nodes[15] = 1.0 ; nodes[16] = 0.0 ; nodes[17] = 1.0 ;
  nodes[18] = 1.0 ; nodes[19] = 1.0 ; nodes[20] = 1.0 ;
  nodes[21] = 0.0 ; nodes[22] = 1.0 ; nodes[23] = 2.0 ;
  trueCentroid[0] = (0.5*1.0/3.0 + 0.5*1.0/3.0 + 0.5)/2.0;
  trueCentroid[1] = (0.5*0.5 + 0.5*0.5 + 0.5)/2.0;
  trueCentroid[2] = (0.5*4.0/3.0 - 0.5*1.0/3.0 + 0.5)/2.0;
  trueVolume = 2.0;
  hexCentroidAndVolume(&nodes[0], &centroid[0], &volume[0]);

  // Nonplaner face, which is approximated with four triangles
  nodes[0] = 0.0  ; nodes[1] = 0.0  ; nodes[2] = 0.0 ;
  nodes[3] = 1.0  ; nodes[4] = 0.0  ; nodes[5] = 0.0 ;
  nodes[6] = 1.0  ; nodes[7] = 1.0  ; nodes[8] = 0.0 ;
  nodes[9] = 0.0  ; nodes[10] = 1.0 ; nodes[11] = 0.0 ;
  nodes[12] = 0.0 ; nodes[13] = 0.0 ; nodes[14] = 1.0 ;
  nodes[15] = 1.0 ; nodes[16] = 0.0 ; nodes[17] = 1.0 ;
  nodes[18] = 1.0 ; nodes[19] = 1.0 ; nodes[20] = 2.0 ;
  nodes[21] = 0.0 ; nodes[22] = 1.0 ; nodes[23] = 1.0 ;
  // The geometry is a cube + pyramid + two tetrahedra
  double cubeVolume = 1.0;
  double cubeCentroid[3] = {0.5, 0.5, 0.5};
  double pyramidVolume = (1.0/3.0)*1.0*0.25;
  double pyramidCentroid[3] = {0.5, 0.5, 1.0+0.25*0.25};
  double tet1Volume = (1.0/3.0)*0.5*0.5;
  double tet1Centroid[3] = {(1.0+1.0+0.5+1.0)/4.0, (0.0+1.0+0.5+1.0)/4.0, (1.0+1.0+1.25+2.0)/4.0}; 
  double tet2Volume = (1.0/3.0)*0.5*0.5;
  double tet2Centroid[3] = {(0.0+1.0+0.5+1.0)/4.0, (1.0+1.0+0.5+1.0)/4.0, (1.0+1.0+1.25+2.0)/4.0}; 
  trueVolume = cubeVolume + pyramidVolume + tet1Volume + tet2Volume;
  trueCentroid[0] = (cubeCentroid[0]*cubeVolume + pyramidCentroid[0]*pyramidVolume + tet1Centroid[0]*tet1Volume + tet2Centroid[0]*tet2Volume) / trueVolume;
  trueCentroid[1] = (cubeCentroid[1]*cubeVolume + pyramidCentroid[1]*pyramidVolume + tet1Centroid[1]*tet1Volume + tet2Centroid[1]*tet2Volume) / trueVolume;
  trueCentroid[2] = (cubeCentroid[2]*cubeVolume + pyramidCentroid[2]*pyramidVolume + tet1Centroid[2]*tet1Volume + tet2Centroid[2]*tet2Volume) / trueVolume;
  hexCentroidAndVolume(&nodes[0], &centroid[0], &volume[0]);
  TEST_FLOATING_EQUALITY(centroid[0], trueCentroid[0], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[1], trueCentroid[1], relTolerance);
  TEST_FLOATING_EQUALITY(centroid[2], trueCentroid[2], relTolerance);
  TEST_FLOATING_EQUALITY(volume[0], trueVolume, relTolerance);
}

//! Exercise triangleSphereIntersection()
TEUCHOS_UNIT_TEST(GeometryUtils, TriangleSphereIntersection) {

  vector<double*> nodes(3);
  vector<double> n1(3), n2(3), n3(3), sphereCenter(3);
  nodes[0] = &n1[0] ; nodes[1] = &n2[0] ; nodes[2] = &n3[0] ;
  double sphereRadius;
  PeridigmNS::SphereIntersection sphereIntersection;

  nodes[0][0] = 1.0 ; nodes[0][1] = -1.0 ; nodes[0][2] = -1.0;
  nodes[1][0] = 1.0 ; nodes[1][1] =  1.0 ; nodes[1][2] =  0.0;
  nodes[2][0] = 1.0 ; nodes[2][1] = -1.0 ; nodes[2][2] =  1.0;
  sphereCenter[0] = 0.0 ; sphereCenter[1] = 0.0 ; sphereCenter[2] = 0.0;
  sphereRadius = std::sqrt(3) + 1.0e-12;
  
  // -- Cases on simple triangle defined above --

  // Case 1:  All the nodes are within the sphere
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);

  // Case 2:  Some nodes are inside the sphere, and some are outside the sphere
  sphereCenter[0] = 0.0 ; sphereCenter[1] = 1.0 ; sphereCenter[2] = 0.0;
  sphereRadius = 1.0 + 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereCenter[0] = 0.0 ; sphereCenter[1] = -1.0 ; sphereCenter[2] = 0.0;
  sphereRadius = std::sqrt(2) + 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // Case 3:  The sphere does not intersect the place of the triangle
  sphereRadius = 1.0 - 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  sphereCenter[0] = 1.1 ; sphereCenter[1] = 0.5 ; sphereCenter[2] = 0.1;
  sphereRadius = 0.1 - 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);

  // Case 4:  None of the nodes are with the sphere, the sphere intersects the (infinite) plane, and the sphere intersects the triangle
  sphereCenter[0] = 0.0 ; sphereCenter[1] = 0.0 ; sphereCenter[2] = 0.0;
  sphereRadius = 1.0 + 1.0e-3;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereCenter[0] = 1.1 ; sphereCenter[1] = 0.5 ; sphereCenter[2] = 0.1;
  sphereRadius = 0.1 + 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  // In this case, the projection of the sphere center onto the plane of the triangle is not in the triangle
  sphereCenter[0] = 0.0 ; sphereCenter[1] = -2.0 ; sphereCenter[2] = 0.0;
  sphereRadius = std::sqrt(2) + 1.0e-3;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // Case 5:  None of the nodes are with the sphere, the sphere intersects the (infinite) plane, but the sphere does not intersect the triangle
  sphereCenter[0] = 0.0 ; sphereCenter[1] = 0.0 ; sphereCenter[2] = 2.0;
  sphereRadius = 1.0 + 1.0e-3;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  sphereCenter[0] = 1.1 ; sphereCenter[1] = 1.2 ; sphereCenter[2] = 0.1;
  sphereRadius = 0.1 + 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);

  // -- A few trickier cases --
  
  nodes[0][0] = 20.0  ; nodes[0][1] = -54.0 ; nodes[0][2] = -31.0;
  nodes[1][0] = 1.0   ; nodes[1][1] = -43.0 ; nodes[1][2] = -57.0;
  nodes[2][0] = -96.0 ; nodes[2][1] = -41.0 ; nodes[2][2] =  72.0;
  double triangleCentroid[3] = {0.0, 0.0, 0.0};
  for(int i=0 ; i<3 ; ++i)
    for(int j=0 ; j<3 ; ++j)
      triangleCentroid[j] += (1.0/3.0)*nodes[i][j];
  double u[3], v[3], triangleNormal[3];
  for(int i=0 ; i<3 ; ++i){
    u[i] = nodes[1][i] - nodes[0][i];
    v[i] = nodes[2][i] - nodes[0][i];
  }
  triangleNormal[0] = u[1]*v[2] - u[2]*v[1];
  triangleNormal[1] = u[2]*v[0] - u[0]*v[2];
  triangleNormal[2] = u[0]*v[1] - u[1]*v[0];
  double mag = std::sqrt(triangleNormal[0]*triangleNormal[0] + triangleNormal[1]*triangleNormal[1] + triangleNormal[2]*triangleNormal[2]);
  for(int i=0 ; i<3 ; ++i)
    triangleNormal[i] /= mag;
  sphereRadius = 33.0;
  for(int i=0 ; i<3 ; ++i)
    sphereCenter[i] = triangleCentroid[i] + sphereRadius*triangleNormal[i];
  sphereRadius = 33.0 - 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  sphereRadius = 33.0 + 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 1.0e30;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);

  sphereRadius = 33.0;
  for(int i=0 ; i<3 ; ++i)
    sphereCenter[i] = triangleCentroid[i] - sphereRadius*triangleNormal[i];
  sphereRadius = 33.0 - 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  sphereRadius = 33.0 + 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 1.0e30;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);

  sphereRadius = 33.0;
  for(int i=0 ; i<3 ; ++i)
    sphereCenter[i] = nodes[0][i] + sphereRadius*triangleNormal[i];
  sphereRadius = 33.0 - 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  sphereRadius = 33.0 + 1.0e-12;
  sphereIntersection = triangleSphereIntersection(nodes, sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
}

//! Exercise tetrahedronSphereIntersection()
TEUCHOS_UNIT_TEST(GeometryUtils, TetrahedronSphereIntersection) {

  vector<double> nodes(12), sphereCenter(3);
  double sphereRadius;
  PeridigmNS::SphereIntersection sphereIntersection;

  // Simple case
  nodes[0]  = 0.0 ; nodes[1]  = 0.0 ; nodes[2]  = 0.0 ;
  nodes[3]  = 1.0 ; nodes[4]  = 0.0 ; nodes[5]  = 0.0 ;
  nodes[6]  = 0.0 ; nodes[7]  = 1.0 ; nodes[8]  = 0.0 ;
  nodes[9]  = 0.0 ; nodes[10] = 0.0 ; nodes[11] = 1.0 ;
  sphereCenter[0] = -1.0 ; sphereCenter[1] = 0.0 ; sphereCenter[2] = 0.0 ;
  sphereRadius = 1.0 - 1.0e-12;
  
  // All nodes outside sphere, no intersection
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);

  // All nodes inside sphere
  sphereRadius = 123.0;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);

  // One node inside sphere
  sphereRadius = 1.0 + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // Three nodes inside sphere
  sphereCenter[0] = 1.0 ; sphereCenter[1] = 1.0 ; sphereCenter[2] = 1.0 ;
  sphereRadius = std::sqrt(2) + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // All nodes outside sphere, but there is an intersection (or not)
  // Kisses first face
  sphereCenter[0] = -1.0 ; sphereCenter[1] = 0.25 ; sphereCenter[2] = 0.25 ;
  sphereRadius = 1.0 + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 1.0 - 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  // Kisses second face
  sphereCenter[0] = 0.25 ; sphereCenter[1] = -1.0 ; sphereCenter[2] = 0.25 ;
  sphereRadius = 1.0 + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 1.0 - 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  // Kisses third face
  sphereCenter[0] = 0.25 ; sphereCenter[1] = 0.25 ; sphereCenter[2] = -1.0 ;
  sphereRadius = 1.0 + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 1.0 - 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  // Kissed fourth face
  sphereCenter[0] = 1.0/3.0 + std::sqrt(1.0/3.0) ; sphereCenter[1] = 1.0/3.0 + + std::sqrt(1.0/3.0) ; sphereCenter[2] = 1.0/3.0 + + std::sqrt(1.0/3.0) ;
  sphereRadius = 1.0 + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 1.0 - 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);

  // Cases where edge kisses sphere (or not)
  nodes[0]  = 1.1 ; nodes[1]  = 0.0 ; nodes[2]  = 1.1 ;
  nodes[3]  = 0.0 ; nodes[4]  = 1.1 ; nodes[5]  = 1.1 ;
  nodes[6]  = 0.0 ; nodes[7]  = 0.0 ; nodes[8]  = 1.0e8 ;
  nodes[9]  = 1.0e8 ; nodes[10] = 1.0e7 ; nodes[11] = 0.0 ;
  sphereCenter[0] = 0.25 ; sphereCenter[1] = 0.25 ; sphereCenter[2] = 0.0 ;
  double X1X0[3] = {sphereCenter[0] - nodes[0], sphereCenter[1] - nodes[1], sphereCenter[2] - nodes[2]};
  double X2X0[3] = {sphereCenter[0] - nodes[3], sphereCenter[1] - nodes[4], sphereCenter[2] - nodes[5]};
  double X2X1[3] = {nodes[3] - nodes[0], nodes[4] - nodes[1], nodes[5] - nodes[2]};
  double X1X0_Cross_X2X0[3] = {X1X0[1]*X2X0[2] - X1X0[2]*X2X0[1], X1X0[2]*X2X0[0] - X1X0[0]*X2X0[2], X1X0[0]*X2X0[1] - X1X0[1]*X2X0[0]};
  double shortestDistanceNumerator = std::sqrt(X1X0_Cross_X2X0[0]*X1X0_Cross_X2X0[0] + X1X0_Cross_X2X0[1]*X1X0_Cross_X2X0[1] + X1X0_Cross_X2X0[2]*X1X0_Cross_X2X0[2]);
  double shortestDistanceDenominator = std::sqrt(X2X1[0]*X2X1[0] + X2X1[1]*X2X1[1] + X2X1[2]*X2X1[2]);
  double shortestDistance = shortestDistanceNumerator/shortestDistanceDenominator;
  sphereRadius = shortestDistance + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = shortestDistance - 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  nodes[0]  = 8.1 ; nodes[1]  = 0.0 ; nodes[2]  = 1.1 ;
  nodes[3]  = 0.0 ; nodes[4]  = 1.1 ; nodes[5]  = -0.1 ;
  nodes[6]  = 100.0 ; nodes[7]  = 100.0 ; nodes[8]  = 0.0 ;
  nodes[9]  = 100.0 ; nodes[10] = 100.0 ; nodes[11] = 0.01;
  sphereCenter[0] = 0.25 ; sphereCenter[1] = 0.25 ; sphereCenter[2] = 0.0 ;
  X1X0[0] = sphereCenter[0] - nodes[0];
  X1X0[1] = sphereCenter[1] - nodes[1];
  X1X0[2] = sphereCenter[2] - nodes[2];
  X2X0[0] = sphereCenter[0] - nodes[3];
  X2X0[1] = sphereCenter[1] - nodes[4];
  X2X0[2] = sphereCenter[2] - nodes[5];
  X2X1[0] = nodes[3] - nodes[0];
  X2X1[1] = nodes[4] - nodes[1];
  X2X1[2] = nodes[5] - nodes[2];
  X1X0_Cross_X2X0[0] = X1X0[1]*X2X0[2] - X1X0[2]*X2X0[1];
  X1X0_Cross_X2X0[1] = X1X0[2]*X2X0[0] - X1X0[0]*X2X0[2];
  X1X0_Cross_X2X0[2] = X1X0[0]*X2X0[1] - X1X0[1]*X2X0[0];
  shortestDistanceNumerator = std::sqrt(X1X0_Cross_X2X0[0]*X1X0_Cross_X2X0[0] + X1X0_Cross_X2X0[1]*X1X0_Cross_X2X0[1] + X1X0_Cross_X2X0[2]*X1X0_Cross_X2X0[2]);
  shortestDistanceDenominator = std::sqrt(X2X1[0]*X2X1[0] + X2X1[1]*X2X1[1] + X2X1[2]*X2X1[2]);
  shortestDistance = shortestDistanceNumerator/shortestDistanceDenominator;
  sphereRadius = shortestDistance + 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = shortestDistance - 1.0e-12;
  sphereIntersection = tetrahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
}

//! Exercise hexahedronSphereIntersection()
TEUCHOS_UNIT_TEST(GeometryUtils, HexahedronSphereIntersection) {

  vector<double> nodes(24), sphereCenter(3);
  double sphereRadius;
  PeridigmNS::SphereIntersection sphereIntersection;

  // Simple case
  nodes[0]  = 0.0 ; nodes[1]  = 0.0 ; nodes[2]  = 0.0 ;
  nodes[3]  = 1.0 ; nodes[4]  = 0.0 ; nodes[5]  = 0.0 ;
  nodes[6]  = 1.0 ; nodes[7]  = 1.0 ; nodes[8]  = 0.0 ;
  nodes[9]  = 0.0 ; nodes[10] = 1.0 ; nodes[11] = 0.0 ;
  nodes[12] = 0.0 ; nodes[13] = 0.0 ; nodes[14] = 1.0 ;
  nodes[15] = 1.0 ; nodes[16] = 0.0 ; nodes[17] = 1.0 ;
  nodes[18] = 1.0 ; nodes[19] = 1.0 ; nodes[20] = 1.0 ;
  nodes[21] = 0.0 ; nodes[22] = 1.0 ; nodes[23] = 1.0 ;
  sphereCenter[0] = -1.0 ; sphereCenter[1] = 0.0 ; sphereCenter[2] = 0.0 ;
  sphereRadius = 1.0 - 1.0e-12;
  
  // All nodes outside sphere, no intersection
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);

  // All nodes inside sphere
  sphereRadius = 123.0;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);

  // One node inside sphere
  sphereRadius = 1.0 + 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // All nodes outside sphere, but there is an intersection
  sphereCenter[0] = -1.0 ; sphereCenter[1] = 0.5 ; sphereCenter[2] = 0.5 ;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // Four nodes inside sphere
  sphereRadius = 2.0 - 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // Slightly trickier case #1
  nodes[0]  = 0.0 ; nodes[1]  = 0.0 ; nodes[2]  = 0.0 ;
  nodes[3]  = 1.0 ; nodes[4]  = 0.0 ; nodes[5]  = 0.0 ;
  nodes[6]  = 1.0 ; nodes[7]  = 1.0 ; nodes[8]  = 0.0 ;
  nodes[9]  = 0.0 ; nodes[10] = 1.0 ; nodes[11] = 0.0 ;
  nodes[12] = 0.0 ; nodes[13] = 0.0 ; nodes[14] = 2.0 ;
  nodes[15] = 1.0 ; nodes[16] = 0.0 ; nodes[17] = 1.0 ;
  nodes[18] = 1.0 ; nodes[19] = 1.0 ; nodes[20] = 1.0 ;
  nodes[21] = 0.0 ; nodes[22] = 1.0 ; nodes[23] = 2.0 ;
  sphereCenter[0] = 1.0 ; sphereCenter[1] = 0.5 ; sphereCenter[2] = 2.0 ;
  sphereRadius = 0.5*std::sqrt(2.0) + 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 0.5*std::sqrt(2.0) - 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  sphereRadius = std::sqrt(5.25) - 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = std::sqrt(5.25) + 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);
  // Case with center of sphere inside the element, sphere touches element faces but no nodes are within sphere
  sphereCenter[0] = 0.5 ; sphereCenter[1] = 0.5 ; sphereCenter[2] = 0.5 ;
  sphereRadius = 0.5 + 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);

  // Slightly trickier case #2
  nodes[0] = -6.0 ; nodes[1] =  8.0 ; nodes[2] =  1.0 ;
  nodes[3] = -1.0 ; nodes[4] =  2.0 ; nodes[5] = -8.0 ;
  nodes[6] =  0.0 ; nodes[7] = -4.0 ; nodes[8] = -4.0 ;
  double edge1[3], edge2[3];
  for(int i=0 ; i<3 ; ++i){
    edge1[i] = nodes[3+i] - nodes[i];
    edge2[i] = nodes[6+i] - nodes[i];
  }
  nodes[9]  = nodes[6] - 2.0*edge1[0] ; nodes[10] = nodes[7] - 2.0*edge1[1]; nodes[11] = nodes[8] - 2.0*edge1[2];
  nodes[12] = nodes[0] + 5.0 ; nodes[13] = nodes[1] ; nodes[14] = nodes[2] ;
  nodes[15] = nodes[3] + 5.0 ; nodes[16] = nodes[4] ; nodes[17] = nodes[5] ;
  nodes[18] = nodes[4] + 1.1 ; nodes[19] = nodes[5] ; nodes[20] = nodes[6] ;
  nodes[21] = nodes[7] + 1.1 ; nodes[22] = nodes[8] ; nodes[23] = nodes[9] ;
  sphereCenter[0] = 0.0 ; sphereCenter[1] = 0.0 ; sphereCenter[2] = 0.0 ;
  sphereRadius = 1234.0;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);
  sphereRadius = sqrt(32.0) + 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  for(int i=0 ; i<3 ; ++i)
    sphereCenter[i] = 0.5*(nodes[i] + nodes[3+i]);
  sphereRadius = 0.1;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  double normal[3];
  normal[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1];
  normal[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2];
  normal[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0];
  double mag = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  for(int i=0 ; i<3 ; ++i)
    normal[i] /= mag;
  for(int i=0 ; i<3 ; ++i)
    sphereCenter[i] += 1.123*normal[i];
  sphereRadius = 1.123 + 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INTERSECTS_SPHERE);
  sphereRadius = 1.123 - 1.0e-12;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::OUTSIDE_SPHERE);
  sphereRadius = 100.0;
  sphereIntersection = hexahedronSphereIntersection(&nodes[0], sphereCenter, sphereRadius);
  TEST_EQUALITY(sphereIntersection, PeridigmNS::INSIDE_SPHERE);
}

int main( int argc, char* argv[] ) {

    int numProcs = 1;
    int returnCode = -1;
   
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
   
    if(numProcs == 1){
       returnCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
    }
    else{
       std::cerr << "Unit test runtime ERROR: utPeridigm_GeometryUtils only makes sense on 1 processor." << std::endl;
    }

    return returnCode;
}

