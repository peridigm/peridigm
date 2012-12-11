/*! \file Peridigm_SearchTree.hpp */
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
#ifndef PERIDIGM_SEARCHTREE_HPP
#define PERIDIGM_SEARCHTREE_HPP

namespace PeridigmNS {

  class SearchTree {

  public:

    /** \brief Constructor.
     *
     *  \param numPoint     The number of points within the tree.
     *  \param coordinates  The coordinates of all the points in the tree, stored as (X0, Y0, Z0, X1, Y1, Z1, ..., XN, YN, ZN).
     **/
    SearchTree(int numPoints, const double* const coordinates);

    //! Destructor.
    virtual ~SearchTree(){}

    /** \brief Finds the set of points within a given radius of a given point.
     *
     *  \param point         The coordinates of the point at the center of the search sphere; this is an array of length three, (X, Y, Z).
     *  \param searchRadius  The radius defining the search sphere.
     *  \param neighborList  The list of ids for all points found within the search sphere; input as an empty list and filled by this function.
     *
     *  This function searches all the points provided to the constructor and returns the ids of those point that are within a
     *  sphere defined by the arguments point and searchRadius.  The ids refer to the positions of the points in the array supplied
     *  to the constructor.  Ids start at zero and increase as (X0, Y0, Z0, X1, Y1, Z1, ..., XN, YN, ZN).
     *
     *  For efficiency, the neighborList argument should be sized to approximately the size of the final neighbor list.
     **/
    virtual FindPointsWithinRadius(const double* const point, double searchRadius, vector<int>& neighborList) = 0;

  private:

    //! Default constructor is private to prevent use
    SearchTree(){}
  }

}

#endif // PERIDIGM_SEARCHTREE_HPP
