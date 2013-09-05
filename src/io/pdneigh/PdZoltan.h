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

#ifndef PD_ZOLTAN_H_
#define PD_ZOLTAN_H_

#include "zoltan.h"
#include "QuickGridData.h"

namespace PDNEIGH {


/*
 * Re-usable component that creates initializes, and returns a "Zoltan" object
 */
struct Zoltan_Struct * createAndInitializeZoltan(QUICKGRID::QuickGridData& pdGridData);

/*
 * Load balancing given a pre-computed neighborhood list -- eg for use with PdQuickGrid where
 * the neighborhood is generally pre-computed
 */
QUICKGRID::QuickGridData& getLoadBalancedDiscretization(QUICKGRID::QuickGridData& pdGridData);

/*
 * Zoltan call back functions
 */
int zoltanQuery_numObjectsOnProc
(
		void *pdGridData,
		int *ierr
);

void zoltanQuery_objectList
(
		void *pdGridData,
		int numGids,
		int numLids,
		ZOLTAN_ID_PTR zoltangIds,
		ZOLTAN_ID_PTR zoltanlIds,
		int numWeights,
		float *objectWts,
		int *ierr
);

int zoltanQuery_dimension
(
		void *pdGridData,
		int *ierr
);

void zoltanQuery_gridData
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numObjects,
		ZOLTAN_ID_PTR zoltangIds,
		ZOLTAN_ID_PTR zoltanlIds,
		int dimension,
		double *gridData,
		int *ierr
);

void zoltanQuery_pointSizeInBytes
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numPoints,
		ZOLTAN_ID_PTR zoltangIds,
		ZOLTAN_ID_PTR zoltanlIds,
		int *sizes,
		int *ierr
);

void zoltanQuery_packPointsMultiFunction
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numPoints,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int *dest,
		int *sizes,
		int *idx,
		char *buf,
		int *ierr
);

void zoltanQuery_unPackPointsMultiFunction
(
		void *pdGridData,
		int numGids,
		int numPoints,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		int *sizes,
		int *idx,
		char *buf,
		int *ierr
);

}  // namespace PDNEIGH

#endif // PD_ZOLTAN_H_
