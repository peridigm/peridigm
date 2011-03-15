
#ifndef PD_ZOLTAN_H_
#define PD_ZOLTAN_H_

#include "zoltan.h"
#include "quick_grid/QuickGridData.h"

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
