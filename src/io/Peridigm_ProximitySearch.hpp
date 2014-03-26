/*! \file Peridigm_ProximitySearch.hpp */
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
#ifndef PERIDIGM_PROXIMITYSEARCH_HPP
#define PERIDIGM_PROXIMITYSEARCH_HPP

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <vector>
#include "pdneigh/BondFilter.h"

namespace PeridigmNS {
namespace ProximitySearch {

  template<class T>
  struct NonDeleter{
	void operator()(T* d) {}
  };

  /** \brief Rebalance a neighborhood list.
   *
   *  \param currentOwnedMap          [input]
   *  \param currentOverlapMap        [input]
   *  \param currentNeighborListSize  [input]
   *  \param currentNeighborList      [input]
   *  \param targetOwnedMap           [input]
   *  \param targetOverlapMap         [output]
   *  \param targetNeighborListsize   [output]
   *  \param targetNeighborList       [output] Allocated within this function.
   **/
  void RebalanceNeighborhoodList(Teuchos::RCP<const Epetra_BlockMap> currentOwnedMap,
                                 Teuchos::RCP<const Epetra_BlockMap> currentOverlapMap,
                                 int currentNeighborListSize,
                                 const int* currentNeighborList,
                                 Teuchos::RCP<const Epetra_BlockMap> targetOwnedMap,
                                 Teuchos::RCP<Epetra_BlockMap>& targetOverlapMap,
                                 int& targetNeighborListSize,
                                 int*& targetNeighborList);

    /** \brief Global proximity search.
     *
     *  \param x                 [input]           Set of points; the neighbors of each point will be found from among the other points in the vector.
     *  \param searchRadii       [input]           The radii list defining the search sphere for each point.
     *  \param overlapMap        [output]          Epetra_BlockMap for importing (ghosting) off-processor neighbors detected by search.
     *  \param neighborListSize  [output]          The length of the neighbor list vector.
     *  \param neighborList      [output]          Pointer to the neighbor list containing the number of neighbors for each point and the list of neighbors for each point (indexes into x).
     *  \param bondFilters       [optional input]  Set of bond filters to employ during the proximity search.
     *  \param radiusAddition    [optional input]  An additional length added to each radius defining the search sphere for each point.
     *
     *  The global proximity search finds, for each point in x, all the points that are within the specified search radius.  The search radius is defined separately for
     *  each point.  The neighborList is allocated within this function and becomes the responsibility of the calling routine (i.e., the calling routine is responsible for deallocation).
     *
     *  The vector x is stored as (x_0, y_0, z_0, x_1, x_2, x_3, ..., x_n, y_n, z_n).
     *
     *  The neighborList is stored as (num_neighbors_0, index1_0, index2_0, ..., indexFinal_0, num_neighbors_1, index1_1, index2_1, ..., indexFinal_1, ..., num_neighbors_N, index1_N, index2_N, ..., indexFinal_N)
     **/
  void GlobalProximitySearch(Teuchos::RCP<Epetra_Vector> x,
                             Teuchos::RCP<Epetra_Vector> searchRadii,
                             Teuchos::RCP<Epetra_BlockMap>& overlapMap,
                             int& neighborListSize,
                             int*& neighborList,
                             std::vector< std::tr1::shared_ptr<PdBondFilter::BondFilter> > bondFilters = std::vector< std::tr1::shared_ptr<PdBondFilter::BondFilter> >(),
                             double radiusAddition = 0.0);

}
}

#endif // PERIDIGM_PROXIMITYSEARCH_HPP


