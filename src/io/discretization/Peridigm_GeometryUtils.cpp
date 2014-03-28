/*! \file Peridigm_GeometryUtils.cpp */

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

#include "Peridigm_GeometryUtils.hpp"
#include <Teuchos_Assert.hpp>
using namespace std;

void PeridigmNS::hexCentroidAndVolume(double* const nodeCoordinates,
                                      double* centroid,
                                      double* volume)
{
  // Divide into tets
  
  // Exodus node connectivity for tets (1-based)
  // 4, 5, 7, 8
  // 3, 4, 5, 7
  // 1, 3, 4, 5
  // 3, 5, 6, 7
  // 1, 3, 5, 6
  // 1, 2, 3, 6

  // 0-based
  // 3, 4, 6, 7
  // 2, 3, 4, 6
  // 0, 2, 3, 4
  // 2, 4, 5, 6
  // 0, 2, 4, 5
  // 0, 1, 2, 5

  // 0-based with stride of 3
  // 9, 12, 18, 21
  // 6, 9, 12, 18
  // 0, 6, 9, 12
  // 6, 12, 15, 18
  // 0, 6, 12, 15
  // 0, 3, 6, 15

  std::vector<double*> tetrahedronNodeCoordinates(4);
  vector<double> tetrahedronCentroid(3), hexahedronCentroid(3);
  double tetrahedronVolume(0.0);

  *volume = 0.0;

  // 9, 12, 18, 21
  tetrahedronNodeCoordinates[0] = nodeCoordinates + 9;
  tetrahedronNodeCoordinates[1] = nodeCoordinates + 12;
  tetrahedronNodeCoordinates[2] = nodeCoordinates + 18;
  tetrahedronNodeCoordinates[3] = nodeCoordinates + 21;
  tetrahedronVolume = tetVolume(tetrahedronNodeCoordinates);
  *volume += tetrahedronVolume;
  if(centroid != 0){
    tetCentroid(tetrahedronNodeCoordinates, tetrahedronCentroid);
    for(int i=0 ; i<3 ; ++i)
      hexahedronCentroid[i] += tetrahedronCentroid[i]*tetrahedronVolume;
  }

  // 6, 9, 12, 18
  tetrahedronNodeCoordinates[0] = nodeCoordinates + 6;
  tetrahedronNodeCoordinates[1] = nodeCoordinates + 9;
  tetrahedronNodeCoordinates[2] = nodeCoordinates + 12;
  tetrahedronNodeCoordinates[3] = nodeCoordinates + 18;
  tetrahedronVolume = tetVolume(tetrahedronNodeCoordinates);
  *volume += tetrahedronVolume;
  if(centroid != 0){
    tetCentroid(tetrahedronNodeCoordinates, tetrahedronCentroid);
    for(int i=0 ; i<3 ; ++i)
      hexahedronCentroid[i] += tetrahedronCentroid[i]*tetrahedronVolume;
  }

  // 0, 6, 9, 12
  tetrahedronNodeCoordinates[0] = nodeCoordinates + 0;
  tetrahedronNodeCoordinates[1] = nodeCoordinates + 6;
  tetrahedronNodeCoordinates[2] = nodeCoordinates + 9;
  tetrahedronNodeCoordinates[3] = nodeCoordinates + 12;
  tetrahedronVolume = tetVolume(tetrahedronNodeCoordinates);
  *volume += tetrahedronVolume;
  if(centroid != 0){
    tetCentroid(tetrahedronNodeCoordinates, tetrahedronCentroid);
    for(int i=0 ; i<3 ; ++i)
      hexahedronCentroid[i] += tetrahedronCentroid[i]*tetrahedronVolume;
  }

  // 6, 12, 15, 18
  tetrahedronNodeCoordinates[0] = nodeCoordinates + 6;
  tetrahedronNodeCoordinates[1] = nodeCoordinates + 12;
  tetrahedronNodeCoordinates[2] = nodeCoordinates + 15;
  tetrahedronNodeCoordinates[3] = nodeCoordinates + 18;
  tetrahedronVolume = tetVolume(tetrahedronNodeCoordinates);
  *volume += tetrahedronVolume;
  if(centroid != 0){
    tetCentroid(tetrahedronNodeCoordinates, tetrahedronCentroid);
    for(int i=0 ; i<3 ; ++i)
      hexahedronCentroid[i] += tetrahedronCentroid[i]*tetrahedronVolume;
  }

  // 0, 6, 12, 15
  tetrahedronNodeCoordinates[0] = nodeCoordinates + 0;
  tetrahedronNodeCoordinates[1] = nodeCoordinates + 6;
  tetrahedronNodeCoordinates[2] = nodeCoordinates + 12;
  tetrahedronNodeCoordinates[3] = nodeCoordinates + 15;
  tetrahedronVolume = tetVolume(tetrahedronNodeCoordinates);
  *volume += tetrahedronVolume;
  if(centroid != 0){
    tetCentroid(tetrahedronNodeCoordinates, tetrahedronCentroid);
    for(int i=0 ; i<3 ; ++i)
      hexahedronCentroid[i] += tetrahedronCentroid[i]*tetrahedronVolume;
  }

  // 0, 3, 6, 15
  tetrahedronNodeCoordinates[0] = nodeCoordinates + 0;
  tetrahedronNodeCoordinates[1] = nodeCoordinates + 3;
  tetrahedronNodeCoordinates[2] = nodeCoordinates + 6;
  tetrahedronNodeCoordinates[3] = nodeCoordinates + 15;
  tetrahedronVolume = tetVolume(tetrahedronNodeCoordinates);
  *volume += tetrahedronVolume;
  if(centroid != 0){
    tetCentroid(tetrahedronNodeCoordinates, tetrahedronCentroid);
    for(int i=0 ; i<3 ; ++i)
      hexahedronCentroid[i] += tetrahedronCentroid[i]*tetrahedronVolume;
  }

  if(centroid != 0){
    for(int i=0 ; i<3 ; ++i)
      centroid[i] = hexahedronCentroid[i]/(*volume);
  }
}

void PeridigmNS::hexVolume(double* const nodeCoordinates, double* volume){
  hexCentroidAndVolume(nodeCoordinates, 0, volume);
}

void PeridigmNS::tetCentroid(const std::vector<double*>& nodeCoordinates,
                             std::vector<double>& centroid)
{
  centroid[0] = (nodeCoordinates[0][0] + nodeCoordinates[1][0] + nodeCoordinates[2][0] + nodeCoordinates[3][0])*0.25;
  centroid[1] = (nodeCoordinates[0][1] + nodeCoordinates[1][1] + nodeCoordinates[2][1] + nodeCoordinates[3][1])*0.25;
  centroid[2] = (nodeCoordinates[0][2] + nodeCoordinates[1][2] + nodeCoordinates[2][2] + nodeCoordinates[3][2])*0.25;
}

double PeridigmNS::tetVolume(const std::vector<double*>& nodeCoordinates)
{
  // Change the coordinate system such that the first point is at the origin.
  // The other three points are labeled a, b, and c.
  std::vector<double> a(3), b(3), c(3);
  for(int dof=0 ; dof<3 ; ++dof){
    a[dof] = nodeCoordinates[1][dof] - nodeCoordinates[0][dof];
    b[dof] = nodeCoordinates[2][dof] - nodeCoordinates[0][dof];
    c[dof] = nodeCoordinates[3][dof] - nodeCoordinates[0][dof];
  }

  // The volume is then | a . (b x c) | / 6
  double volume = scalarTripleProduct(a, b, c) / 6.0;

  TEUCHOS_TEST_FOR_EXCEPT_MSG(volume < 0.0, "\n**** Error:  tetVolume() computed negative volume, possible element inversion or incorrect node numbering.\n");

  return volume;
}

PeridigmNS::SphereIntersection PeridigmNS::triangleSphereIntersection(const std::vector<double*>& nodeCoordinates,
                                                                      const std::vector<double>& sphereCenter,
                                                                      double sphereRadius)
{
  // This algorithm was found on Christer Ericson's blog
  // He's the author of Real-Time Collision Detection (Morgan Kaufmann, 2005).
  
  // The algorithm is based on the separating axis test.
  // A separating axis is a line, parallel to which a plane
  // exists that separates two entities.  In this case,
  // the entities are the sphere and the triangle.  Separating
  // axes may exist between the sphere center and the triangle plane,
  // the sphere center and the triangle vertices, or the sphere center
  // and the triangle edges.  Each of these possibilities is examined,
  // and if no separating axis is found, it's concluded that the
  // sphere and the triangle intersect.

  // Translate problem so that the sphere center is the origin
  // The three nodes of the (translated) triangle are A, B, C
  // The sphere origin is P
  double A[3], B[3], C[3];
  subtract(nodeCoordinates[0], &sphereCenter[0], A);
  subtract(nodeCoordinates[1], &sphereCenter[0], B);
  subtract(nodeCoordinates[2], &sphereCenter[0], C);
  double P[3] = {0.0, 0.0, 0.0};

  // Raduis squared
  double rr = sphereRadius*sphereRadius;
  
  // Squared distances of each node from the sphere center
  double aa, bb, cc;
  dot(A, A, &aa);
  dot(B, B, &bb);
  dot(C, C, &cc);

  // Check how many nodes are within the sphere
  int numWithinSphere(0);
  if(aa < rr)
    numWithinSphere += 1;
  if(bb < rr)
    numWithinSphere += 1;
  if(cc < rr)
    numWithinSphere += 1;

  // Check for all-in or intersection conditions
  if(numWithinSphere == 1 || numWithinSphere == 2)
    return INTERSECTS_SPHERE;
  if(numWithinSphere == 3)
    return INSIDE_SPHERE;

  // If the function didn't just exit, all the nodes must be outside the sphere.

  // Store the triangle edges
  double AB[3], BC[3], CA[3];
  subtract(B, A, AB);
  subtract(C, B, BC);
  subtract(A, C, CA);

  // N is the plane normal (not normalized to avoid sqrt calculation)
  double N[3], temp[3];
  subtract(C, A, temp); // temp is AC
  cross(AB, temp, N);
  
  // d is the shortest distance from P to the plane of the triangle
  // N is not normalized, so to get the true distance d would have to be divided by the magnitude of N
  double d;
  dot(A, N, &d);
  
  // e is the squared magnitude of plane normal
  double e;
  dot(N, N, &e);

  // If the shortest distance from the sphere center to the plane is larger than the sphere radius, there is a separation.
  // Note that d*d must be divided by e because a normalization step was skipped above (to avoid sqrt calculation)
  if( d*d/e > rr )
    return OUTSIDE_SPHERE;

  // The triangle is outside the sphere if the following two conditions are met:
  //   1) Node A is outside the sphere (we already know it is from a check above)
  //   2) Node B and the center of the sphere are on opposite sides of the plane
  //      defined by the point A and the vector between A and the sphere center
  //      (because we moved the sphere center to the origin, this vector is just A)
  //   3) Node C and the center of the sphere are on opposite sides of the plane
  //      defined by the point A and the vector between A and the sphere center
  double ab, ac, bc;
  dot(A, B, &ab);
  dot(A, C, &ac);
  dot(B, C, &bc);
  // Check to see if both points are on the opposite side of the plane from the origin
  // Because the plane normal is A, the distance to the origin will be negative.
  // So B and C will be on opposite sides of the plane if their distances are positive.
  //
  // We could perform this check:
  //   dot(AB, A, check1);
  //   dot(AC, A, check2);
  //   if(check1 > 0 and check2 > 0)
  //     return OUTSIDE_SPHERE;
  //
  // Some calculations can be eliminated by distributing the dot product (AB is defined as B-A)
  //   dot(AB,A) -> dot(AB) - dot(AA)
  // Rearranging yields
  if(ab > aa && ac > aa)
    return OUTSIDE_SPHERE;
  //
  // Repeat the check for planes defined at points B and C
  if(ab > bb && bc > bb)
    return OUTSIDE_SPHERE;
  if(ac > cc && bc > cc)
    return OUTSIDE_SPHERE;

  // TODO:  RE-USE VARAIBLES WHERE POSSIBLE BELOW

  // The final checks for possible separating axes involve the edges.
  //
  // Construct a plane that contains the edge AB and has a normal that
  // passes through the sphere center.
  //
  // Consider the point Q on the (infinite) line AB to which P projects.
  double AP[3], Q[3];
  double ap_ab, ab_ab, t;
  subtract(P, A, AP);
  dot(AP, AB, &ap_ab);
  dot(AB, AB, &ab_ab);
  t = ap_ab/ab_ab;
  for(int i=0 ; i<3 ; ++i)
    Q[i] = A[i] + t*AB[i];
  // There is a separation if the following conditions are both true:
  //   1) Q is a distance greater than the radius from the sphere center
  //   2) The points C and the sphere center lie on opposite sides of the plane
  double QP[3], QC[3];
  double qp_qp, qp_qc;
  // The first condition is true if the squared length of QP is greater than the radius squared, (qp_qp > rr)
  subtract(P, Q, QP);
  dot(QP, QP, &qp_qp);
  // The second condition is true if the vectors QP and QC point in opposite directions (dot product qp_qc is negative)
  subtract(C, Q, QC);
  dot(QP, QC, &qp_qc);
  if( (qp_qp > rr) && (qp_qc < 0.0) )
    return OUTSIDE_SPHERE;
  //
  // Now consider edge BC
  double BP[3];
  double bp_bc, bc_bc;
  subtract(P, B, BP);
  dot(BP, BC, &bp_bc);
  dot(BC, BC, &bc_bc);
  t = bp_bc/bc_bc;
  for(int i=0 ; i<3 ; ++i)
    Q[i] = B[i] + t*BC[i];
  double QA[3];
  double qp_qa;
  subtract(P, Q, QP);
  dot(QP, QP, &qp_qp);
  subtract(A, Q, QA);
  dot(QP, QA, &qp_qa);
  if( (qp_qp > rr) && (qp_qa < 0.0) )
    return OUTSIDE_SPHERE;
  //
  // Now consider edge CA
  double CP[3];
  double cp_ca, ca_ca;
  subtract(P, C, CP);
  dot(CP, CA, &cp_ca);
  dot(CA, CA, &ca_ca);
  t = cp_ca/ca_ca;
  for(int i=0 ; i<3 ; ++i)
    Q[i] = C[i] + t*CA[i];
  double QB[3];
  double qp_qb;
  subtract(P, Q, QP);
  dot(QP, QP, &qp_qp);
  subtract(B, Q, QB);
  dot(QP, QB, &qp_qb);
  if( (qp_qp > rr) && (qp_qb < 0.0) )
    return OUTSIDE_SPHERE;

  // All possible separating axes have been eliminated, so there must be an intersection
  return INTERSECTS_SPHERE;
}

PeridigmNS::SphereIntersection PeridigmNS::hexahedronSphereIntersection(double* const nodeCoordinates,
                                                                        const std::vector<double>& sphereCenter,
                                                                        double sphereRadius)
{
  // Perform check on the number of nodes within the sphere
  int numNodesInSphere(0);
  double radiusSquared, distanceSquared;
  radiusSquared = sphereRadius*sphereRadius;
  for(int nodeIndex=0 ; nodeIndex<24 ; nodeIndex+=3){
    distanceSquared =
      (nodeCoordinates[nodeIndex]   - sphereCenter[0])*(nodeCoordinates[nodeIndex]   - sphereCenter[0]) +
      (nodeCoordinates[nodeIndex+1] - sphereCenter[1])*(nodeCoordinates[nodeIndex+1] - sphereCenter[1]) +
      (nodeCoordinates[nodeIndex+2] - sphereCenter[2])*(nodeCoordinates[nodeIndex+2] - sphereCenter[2]);
    if(distanceSquared < radiusSquared)
      numNodesInSphere += 1;
  }

  if(numNodesInSphere == 8)
    return INSIDE_SPHERE;
  else if(numNodesInSphere > 0)
    return INTERSECTS_SPHERE;

  // If all the nodes are outside the sphere, we need to look carefully at each face

  // Divide each face into two triangles
  // Check for intersection of triangles and sphere  

  SphereIntersection sphereIntersection;
  bool allInside(true), allOutside(true);
  std::vector<double*> triangleNodeCoordinates(3);

  for(int face=0 ; face<6 ; face++){
    for(int tri=0 ; tri<2 ; ++tri){

      if(face == 0 and tri == 0){
        // Exodus nodes 1, 2, 3
        // 0-based      0, 1, 2
        // stride 3     0, 3, 6
        triangleNodeCoordinates[0] = nodeCoordinates;
        triangleNodeCoordinates[1] = nodeCoordinates + 3;
        triangleNodeCoordinates[2] = nodeCoordinates + 6;
      }
      else if(face == 0 and tri == 1){
        // Exodus nodes 1, 3, 4
        // 0-based      0, 2, 3
        // stride 3     0, 6, 9
        triangleNodeCoordinates[0] = nodeCoordinates;
        triangleNodeCoordinates[1] = nodeCoordinates + 6;
        triangleNodeCoordinates[2] = nodeCoordinates + 9;
      }
      else if(face == 1 and tri == 0){
        // Exodus nodes 5, 6, 7
        // 0-based      4, 5, 6
        // stride 3     12, 15, 18
        triangleNodeCoordinates[0] = nodeCoordinates + 12;
        triangleNodeCoordinates[1] = nodeCoordinates + 15;
        triangleNodeCoordinates[2] = nodeCoordinates + 18;
      }
      else if(face == 1 and tri == 1){
        // Exodus nodes 5, 7, 8
        // 0-based      4, 6, 7
        // stride 3     12, 18, 21
        triangleNodeCoordinates[0] = nodeCoordinates + 12;
        triangleNodeCoordinates[1] = nodeCoordinates + 18;
        triangleNodeCoordinates[2] = nodeCoordinates + 21;
      }
      else if(face == 2 and tri == 0){
        // Exodus nodes 4, 5, 8
        // 0-based      3, 4, 7
        // stride 3     9, 12, 21
        triangleNodeCoordinates[0] = nodeCoordinates + 9;
        triangleNodeCoordinates[1] = nodeCoordinates + 12;
        triangleNodeCoordinates[2] = nodeCoordinates + 21;
      }
      else if(face == 2 and tri == 1){
        // Exodus nodes 1, 4, 5
        // 0-based      0, 3, 4
        // stride 3     0, 9, 12
        triangleNodeCoordinates[0] = nodeCoordinates;
        triangleNodeCoordinates[1] = nodeCoordinates + 9;
        triangleNodeCoordinates[2] = nodeCoordinates + 12;
      }
      else if(face == 3 and tri == 0){
        // Exodus nodes 3, 6, 7
        // 0-based      2, 5, 6
        // stride 3     6, 15, 18
        triangleNodeCoordinates[0] = nodeCoordinates + 6;
        triangleNodeCoordinates[1] = nodeCoordinates + 15;
        triangleNodeCoordinates[2] = nodeCoordinates + 18;
      }
      else if(face == 3 and tri == 1){
        // Exodus nodes 2, 3, 6
        // 0-based      1, 2, 5
        // stride 3     3, 6, 15
        triangleNodeCoordinates[0] = nodeCoordinates + 3;
        triangleNodeCoordinates[1] = nodeCoordinates + 6;
        triangleNodeCoordinates[2] = nodeCoordinates + 15;
      }
      else if(face == 4 and tri == 0){
        // Exodus nodes 1, 2, 6
        // 0-based      0, 1, 5
        // stride 3     0, 3, 15
        triangleNodeCoordinates[0] = nodeCoordinates;
        triangleNodeCoordinates[1] = nodeCoordinates + 3;
        triangleNodeCoordinates[2] = nodeCoordinates + 15;
      }
      else if(face == 4 and tri == 1){
        // Exodus nodes 1, 5, 6
        // 0-based      0, 4, 5
        // stride 3     0, 12, 15
        triangleNodeCoordinates[0] = nodeCoordinates;
        triangleNodeCoordinates[1] = nodeCoordinates + 12;
        triangleNodeCoordinates[2] = nodeCoordinates + 15;
      }
      else if(face == 5 and tri == 0){
        // Exodus nodes 3, 4, 7
        // 0-based      2, 3, 6
        // stride 3     6, 9, 18
        triangleNodeCoordinates[0] = nodeCoordinates + 6;
        triangleNodeCoordinates[1] = nodeCoordinates + 9;
        triangleNodeCoordinates[2] = nodeCoordinates + 18;
      }
      else if(face == 5 and tri == 1){
        // Exodus nodes 4, 7, 8
        // 0-based      3, 6, 7
        // stride 3     9, 18, 21
        triangleNodeCoordinates[0] = nodeCoordinates + 9;
        triangleNodeCoordinates[1] = nodeCoordinates + 18;
        triangleNodeCoordinates[2] = nodeCoordinates + 21;
      }

      sphereIntersection = triangleSphereIntersection(triangleNodeCoordinates, sphereCenter, sphereRadius);
      if(sphereIntersection == INTERSECTS_SPHERE)
        return INTERSECTS_SPHERE;
      else if(sphereIntersection == INSIDE_SPHERE)
        allOutside = false;
      else if(sphereIntersection == OUTSIDE_SPHERE)
        allInside = false;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG(allInside && allOutside, "\n**** Error:  Nonsense result in hexahedronSphereIntersection().\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!allInside && !allOutside, "\n**** Error:  Nonsense result in hexahedronSphereIntersection().\n");

  if(allOutside)
    return OUTSIDE_SPHERE;
  return INSIDE_SPHERE;
}

// void circumscribeElementWithSphere(const std::vector<double*>& nodeCoordinates,
//                                    double* sphereCenter,
//                                    double* radius)
// {
//   if(nodeCoordinates.size
// }

void PeridigmNS::subtract(const double* const a,
                          const double* const b, 
                          double* c)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
  return;
}

void PeridigmNS::dot(const double* const a,
                     const double* const b, 
                     double* c)
{
  *c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return;
}

void PeridigmNS::cross(const double* const a,
                       const double* const b, 
                       double* c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return;
}

double PeridigmNS::scalarTripleProduct(const std::vector<double>& a,
                                       const std::vector<double>& b,
                                       const std::vector<double>& c)
{
  double tripleProduct = 
    a[0]*(b[1]*c[2] - b[2]*c[1]) + a[1]*(b[2]*c[0] - b[0]*c[2]) + a[2]*(b[0]*c[1] - b[1]*c[0]);

  return tripleProduct;
}

double PeridigmNS::maxDistanceToNode(int numNodes,
                                     const double* const nodeCoordinates,
                                     const double* point)
{
  double maxDistance(0.0), squaredDistance(0.0), x(0.0), y(0.0), z(0.0);
  for(int i=0 ; i<numNodes ; ++i){
    x = nodeCoordinates[i*3]   - point[0];
    y = nodeCoordinates[i*3+1] - point[1];
    z = nodeCoordinates[i*3+2] - point[2];
    squaredDistance = x*x + y*y + z*z;
    if(squaredDistance > maxDistance)
      maxDistance = squaredDistance;
  }
  maxDistance = std::sqrt(maxDistance);
  return maxDistance;
}


// TODO
//
// unit test hexahedronSphereIntersection()
// recursive algorithm for computing partial volume
// revist function to find largest element dimension (circumscribing sphere?)
