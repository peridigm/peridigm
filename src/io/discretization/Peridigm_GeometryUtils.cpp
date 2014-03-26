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
  double radiusSquared(0.0), distanceSquared(0.0);
  int numNodesInSphere(0);

  radiusSquared = sphereRadius*sphereRadius;

  for(int node=0 ; node<3 ; ++node){
    distanceSquared =
      (nodeCoordinates[node][0] - sphereCenter[0])*(nodeCoordinates[node][0] - sphereCenter[0]) +
      (nodeCoordinates[node][1] - sphereCenter[1])*(nodeCoordinates[node][1] - sphereCenter[1]) +
      (nodeCoordinates[node][2] - sphereCenter[2])*(nodeCoordinates[node][2] - sphereCenter[2]);
    if(distanceSquared < radiusSquared)
      numNodesInSphere += 1;
  }

  // If all the nodes are in the sphere, then the entire triangle is inside the sphere
  if(numNodesInSphere == 3)
    return INSIDE_SPHERE;

  // If some node are inside the sphere and some are outside the sphere, then the triangle intersects the sphere
  if(numNodesInSphere == 1 || numNodesInSphere == 2)
    return INTERSECTS_SPHERE;

  // Compute the normal to the triangle
  double u[3], v[3], magnitude, normal[3];
  for(int i=0 ; i<3 ; ++i){
    u[i] = nodeCoordinates[1][i] - nodeCoordinates[0][i];
    v[i] = nodeCoordinates[2][i] - nodeCoordinates[0][i];
  }
  normal[0] = u[1]*v[2] - u[2]*v[1];
  normal[1] = u[2]*v[0] - u[0]*v[2];
  normal[2] = u[0]*v[1] - u[1]*v[0];
  magnitude = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  for(int i=0 ; i<3 ; ++i)
    normal[i] /= magnitude;

  // Compute the closest distance from the plane of the triangle to the center of the sphere
  double distance =
    normal[0]*(sphereCenter[0] - nodeCoordinates[0][0]) +
    normal[1]*(sphereCenter[1] - nodeCoordinates[0][1]) +
    normal[2]*(sphereCenter[2] - nodeCoordinates[0][2]);

  if(std::abs(distance) > sphereRadius)
    return OUTSIDE_SPHERE;

  // Find the closest point on the plane to the center of the sphere
  double pt[3];
  for(int i=0 ; i<3 ; ++i)
    pt[i] = sphereCenter[i] - distance*normal[i];

  // Determine if the point is within the triangle
  // Technique is based on barycentric coordinates (http://www.blackpawn.com/texts/pointinpoly/)
  double w[3];
  for(int i=0 ; i<3 ; ++i)
    w[i] = pt[i] - nodeCoordinates[0][i];

  double dot00 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
  double dot01 = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
  double dot02 = u[0]*w[0] + u[1]*w[1] + u[2]*w[2];
  double dot11 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  double dot12 = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];

  double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
  double bc_1 = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double bc_2 = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // Check if point is in triangle
  if( (bc_1 >= 0.0) && (bc_2 >= 0.0) && (bc_1 + bc_2 <= 1.0) )
    return INTERSECTS_SPHERE;
 
  return OUTSIDE_SPHERE;
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
