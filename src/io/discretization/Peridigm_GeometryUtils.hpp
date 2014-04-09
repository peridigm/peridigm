/*! \file Peridigm_GeometryUtils.hpp */

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

#ifndef PERIDIGM_GEOMETRYUTILS_HPP
#define PERIDIGM_GEOMETRYUTILS_HPP

#include <vector>

namespace PeridigmNS {

  enum SphereIntersection { INSIDE_SPHERE = 0, INTERSECTS_SPHERE = 1, OUTSIDE_SPHERE = 2 };

  //! Approximate the volume and centroid of a hexahedron.
  void hexCentroidAndVolume(double* const nodeCoordinates, double* centroid, double* volume);

  //! Approximate the volume of a hexahedron.
  void hexVolume(double* const nodeCoordinates, double* volume);

  //! Compute the centroid of a tetrahedron.
  void tetCentroid(const std::vector<double*>& nodeCoordinates, std::vector<double>& centroid);

  //! Compute the volume of a tetrahedron.
  double tetVolume(const std::vector<double*>& nodeCoordinates);

  //! Determine if a triangle intersects a sphere
  SphereIntersection triangleSphereIntersection(const std::vector<double*>& nodeCoordinates,
                                                const std::vector<double>& sphereCenter,
                                                double sphereRadius);

  //! Determine if a hexahedron intersects a sphere
  SphereIntersection hexahedronSphereIntersection(double* const nodeCoordinates,
                                                  const std::vector<double>& sphereCenter,
                                                  double sphereRadius);

  //! Compute the difference of two three-dimensional vectors
  void subtract(const double* const a, const double* const b, double* c);

  //! Compute dot product of two three-dimensional vectors
  void dot(const double* const a, const double* const b, double* c);

  //! Compute cross product of two three-dimensional vectors
  void cross(const double* const a, const double* const b, double* c);

  //! Compute a scalar triple product.
  double scalarTripleProduct(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);

  //! Compute the maxmimum distance from a given point to a node in an element.
  double maxDistanceToNode(int numNodes, const double* const nodeCoordinates, const double* point);
}

#endif // PERIDIGM_GEOMETRYUTILS_HPP
