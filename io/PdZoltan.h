
#ifndef PD_ZOLTAN_H_
#define PD_ZOLTAN_H_

#include "zoltan.h"
#include "PdGridData.h"
#include<tr1/memory>

using std::tr1::shared_ptr;

/*
 * Re-usable component that creates initializes, and returns a "Zoltan" object
 */
struct Zoltan_Struct * createAndInitializeZoltan(PdGridData& pdGridData);

/*
 * Load balancing given a pre-computed neighborhood list -- eg for use with PdQuickGrid where
 * the neighborhood is generally pre-computed
 */
PdGridData& getLoadBalancedDiscretization(PdGridData& pdGridData);

/*
 * Use this function to perform parallel search and create a new neighborhood list;
 * This function will construct neighborhood lists including accross processor boundaries;
 * The resulting list will then be added to the incoming pdGridData -- note that any pre-existing
 * list on pdGridData will be overwritten/lost.
 * @param horizon -- this is the distance that should be used to form the neighborhood list
 * Use CASE Scenario:
 * 1) Create a mesh
 * 2) Load balance mesh
 * 3) Call this function
 * 4) 4th argument -- withSelf : set true if each point should be part of its own neighborhood
 */
PdGridData& createAndAddNeighborhood(PdGridData& decomp, double horizon, bool withSelf=false);

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

#endif // PD_ZOLTAN_H_
