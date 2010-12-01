#include <tr1/memory>
#include <vector>
#include "zoltan.h"
#include "PdZoltan.h"
#include "PdQuickGrid.h"
#include "PdVTK.h"
#include "PdQuickGridParallel.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "PdNeighborhood.h"
#include <iostream>

using namespace PdQuickGrid;
using namespace PdNeighborhood;
using std::tr1::shared_ptr;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;

struct ZoltanDestroyer{
	void operator()(struct Zoltan_Struct *zoltan) {
		Zoltan_Destroy(&zoltan);
	}
};

/*
 * This function is private to the library and should not be called otherwise
 * This function is called from the zoltan callback: zoltanQuery_unPackPointsMultiFuntion
 * Inputs to this function are those coming from the zoltan call back: zoltanQuery_unPackPointsMultiFuntion
 */
int computeSizeNewNeighborhoodList(int initialValue, int numImport, int *idx, char *buf, int dimension);


struct Zoltan_Struct * createAndInitializeZoltan(PdGridData& pdGridData){

	/*
	 * The Zoltan_Initialize function initializes MPI for Zoltan.
	 * If the application uses MPI, this function should be called after calling MPI_Init.
	 * If the application does not use MPI, this function calls MPI_Init for use by Zoltan.
	 * This function is called with the argc and argv command-line arguments from the main program,
	 * which are used if Zoltan_Initialize calls MPI_Init. From C,  if MPI_Init has already been called,
	 * the argc and argv arguments may have any value because their values will be ignored.
	 * From Fortran, if one of argc or argv is omitted, they must both be omitted.
	 * If they are omitted, ver does NOT have to be passed as a keyword argument.
	 */
	int numArgs=0;
	char **argv=0;
	float version=0;
	int zoltanErr = Zoltan_Initialize(numArgs,argv,&version);

	struct Zoltan_Struct *zoltan;

	/******************************************************************
	 ** Create a Zoltan library structure for this instance of load
	 ** balancing.  Set the parameters and query functions that will
	 ** govern the library's calculation.  See the Zoltan User's
	 ** Guide for the definition of these and many other parameters.
	 ******************************************************************/
	/*
	 * Passing MPI_COMM_WORLD -- all processors should participate
	 */
	zoltan = Zoltan_Create(MPI_COMM_WORLD);

	/* Zoltan general parameters */

	/*
	 *  DEBUG-LEVEL=0 Quiet mode; no output unless an error or warning is produced.
	 */
	Zoltan_Set_Param(zoltan, "DEBUG_LEVEL", "0");
	/*
	 * Load balance method
	 */
	Zoltan_Set_Param(zoltan, "LB_METHOD", "RCB");

	/*
	 * Tell Zoltan to automatically migrate data for the application.  Zoltan_Migrate is called below
	 */
	//	Zoltan_Set_Param(zoltan, "AUTO_MIGRATE", "TRUE");

	/*
	 * The number of unsigned integers that should be used to represent a global identifier (ID).
	 */
	Zoltan_Set_Param(zoltan, "NUM_GID_ENTRIES", "1");
	/*
	 * The number of unsigned integers that should be used to represent a local identifier (ID).
	 */
	Zoltan_Set_Param(zoltan, "NUM_LID_ENTRIES", "1");
	/*
	 * The number of weights (to be supplied by the user in a query function) associated with an object.
	 * If this parameter is zero, all objects have equal weight.
	 */
	Zoltan_Set_Param(zoltan, "OBJ_WEIGHT_DIM", "0");

	/*
	 * Must set this so that we can later call Zoltan_LB_Box_PP_Assign
	 * Zoltan_LB_Box_PP_Assign: determines which processors and parts own geometry that intersects
	 * an input bounding box; In this case, we look at points on the boundary of each processor, put a
	 * bounding box around them (not that BB will have some portion off processor), then call this
	 * function with the "zoltan'.  We will get back the processors that need to share the input point
	 */
	Zoltan_Set_Param(zoltan,"KEEP_CUTS", "1");

	/*
	 * "IMPORT", to return only information about objects to be imported to a processor
	 * "EXPORT", to return only information about objects to be exported from a processor
	 * "ALL", or "IMPORT AND EXPORT" (or any string with both "IMPORT" and "EXPORT" in it) to return both import and export information
	 * "PARTITION ASSIGNMENTS" (or any string with "PARTITION" in it) to return the new process and part assignment
	 * of every local object, including those not being exported.
	 */
	Zoltan_Set_Param(zoltan, "RETURN_LISTS", "ALL");

	/* RCB parameters */

	/*
	 * Flag controlling the amount of timing and diagnostic output the routine produces.
	 * 0 = no output; 1 = print summary; 2 = print data for each processor.
	 */
	Zoltan_Set_Param(zoltan, "RCB_OUTPUT_LEVEL", "0");

	/*
	 * Flag controlling the shape of the resulting regions. If this option is specified, then when a cut is made,
	 * all of the dots located on the cut are moved to the same side of the cut. The resulting regions are
	 * then rectilinear. When these dots are treated as a group, then the resulting load balance may not
	 * be as good as when the group of dots is split by the cut.
	 * 0 = move dots individually; 1 = move dots in groups.
	 */
	Zoltan_Set_Param(zoltan, "RCB_RECTILINEAR_BLOCKS", "0");

	/*
	 * query function returns the number of objects that are currently assigned to the processor
	 */
	Zoltan_Set_Num_Obj_Fn(zoltan, zoltanQuery_numObjectsOnProc, &pdGridData);

	/*
	 * query function fills two (three if weights are used) arrays with information about the objects
	 * currently assigned to the processor. Both arrays are allocated (and subsequently freed) by Zoltan;
	 * their size is determined by a call to a ZOLTAN_NUM_OBJ_FN query function to get the array size.
	 */
	Zoltan_Set_Obj_List_Fn(zoltan, zoltanQuery_objectList, &pdGridData);

	/*
	 * query function returns the number of values needed to express the geometry of an object.
	 * For example, for a two-dimensional mesh-based application, (x,y) coordinates are needed
	 * to describe an object's geometry; thus the ZOLTAN_NUM_GEOM_FN query function should return
	 * the value of two. For a similar three-dimensional application, the return value should be three.
	 */
	Zoltan_Set_Num_Geom_Fn(zoltan, zoltanQuery_dimension, &pdGridData);

	/*
	 * query function returns a vector of geometry values for a list of given objects. The geometry
	 * vector is allocated by Zoltan to be of size num_obj * num_dim.
	 */
	Zoltan_Set_Geom_Multi_Fn(zoltan, zoltanQuery_gridData, &pdGridData);


	return zoltan;
}

PdGridData& getLoadBalancedDiscretization(PdGridData& pdGridData){
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) Start"  << std::endl; std::cout.flush();

	struct Zoltan_Struct *zoltan = createAndInitializeZoltan(pdGridData);

//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) A"  << std::endl; std::cout.flush();
	pdGridData.zoltanPtr = shared_ptr<struct Zoltan_Struct>(zoltan,ZoltanDestroyer());

	/* Query functions, to provide geometry to Zoltan */


	/*
	 * Set number of bytes per node -- this is used for migration of data based after load balancing
	 */
	Zoltan_Set_Obj_Size_Multi_Fn(zoltan, zoltanQuery_pointSizeInBytes, &pdGridData);
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) B"  << std::endl; std::cout.flush();
	/*
	 * A ZOLTAN_PACK_OBJ_FN query function allows the application to tell Zoltan how to copy all
	 * needed data for a given object into a communication buffer. The object's data can then be
	 * sent to another processor as part of data migration. It may also perform other operations,
	 * such as removing the object from the processor's data structure. This routine is called by
	 * Zoltan_Migrate for each object to be sent to another processor.
	 */
	Zoltan_Set_Pack_Obj_Multi_Fn(zoltan,zoltanQuery_packPointsMultiFunction,&pdGridData);
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) C"  << std::endl; std::cout.flush();
	/*
	 * A ZOLTAN_UNPACK_OBJ_FN query function allows the application to tell Zoltan how to copy
	 * all needed data for a given object from a communication buffer into the application's data
	 * structure. This operation is needed as the final step of importing objects during data
	 * migration. The query function may also perform other computation, such as building request
	 * lists for related data.
	 */
	Zoltan_Set_Unpack_Obj_Multi_Fn(zoltan,zoltanQuery_unPackPointsMultiFunction,&pdGridData);
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) D"  << std::endl; std::cout.flush();

	/******************************************************************
	 ** Zoltan partition mesh: the number of partitions is
	 ** equal to the number of processes.  Process rank 0 will own
	 ** partition 0, process rank 1 will own partition 1, and so on.
	 ******************************************************************/
	int changes, numGidEntries, numLidEntries, numImport, numExport;
	ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
	int *importProcs, *importToPart, *exportProcs, *exportToPart;

	int zoltanErr = Zoltan_LB_Partition
			(
					zoltan,             /* input (all remaining fields are output) */
					&changes,           /* 1 if partitioning was changed, 0 otherwise */
					&numGidEntries,     /* Number of integers used for a global ID */
					&numLidEntries,     /* Number of integers used for a local ID */
					&numImport,         /* Number of vertices to be sent to me */
					&importGlobalGids,  /* Global IDs of vertices to be sent to me */
					&importLocalGids,   /* Local IDs of vertices to be sent to me */
					&importProcs,       /* Process rank for source of each incoming vertex */
					&importToPart,      /* New partition for each incoming vertex */
					&numExport,         /* Number of vertices I must send to other processes*/
					&exportGlobalGids,  /* Global IDs of the vertices I must send */
					&exportLocalGids,   /* Local IDs of the vertices I must send */
					&exportProcs,       /* Process to which I send each of the vertices */
					&exportToPart       /* Partition to which each vertex will belong */
			);
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) E"  << std::endl; std::cout.flush();
	Zoltan_Migrate
	(
			zoltan,
			numImport,
			importGlobalGids,
			importLocalGids,
			importProcs,
			importToPart,
			numExport,
			exportGlobalGids,
			exportLocalGids,
			exportProcs,
			exportToPart
	);
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) F"  << std::endl; std::cout.flush();
	/*
	 * Now insure that all processors were unpacked
	 * NOTE:
	 * This is a very special case in which a processor decomposition is such that no points are imported; In this case,
	 * the unpack function must be called for that processor -- otherwise PdGridData on that processor will be inconsistent if
	 * that processor exported points.
	 */
	if(pdGridData.unPack){
		ZOLTAN_ID_PTR gIds = 0;
		int numImport = 0;
		int *sizes=0;
		int *idx=0;
		char *buf = 0;
		zoltanQuery_unPackPointsMultiFunction(&pdGridData,numGidEntries,numImport,gIds,sizes,idx,buf,&zoltanErr);
	}
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) G"  << std::endl; std::cout.flush();
	/* Free memory allocated for load-balancing results by Zoltan */
	Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
	Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
	if (zoltanErr != ZOLTAN_OK){
		std::cerr << "PdQuickGrid::getLoadBalancedDiscretization -- Zoltan_LB_Partition Failure.  Abort " << std::endl;
		MPI_Finalize();
//		Zoltan_Destroy(&zoltan);
		exit(0);
	}

//	Zoltan_Destroy(&zoltan);
//	std::cout << "getLoadBalancedDiscretization(PdGridData& pdGridData) Finish" << std::endl;
	return pdGridData;
}


int zoltanQuery_numObjectsOnProc
(
		void *pdGridData,
		int *ierr
)
{
	*ierr = ZOLTAN_OK;
	PdGridData *gridData = (PdGridData *)pdGridData;
	return gridData->numPoints;
}

void zoltanQuery_objectList
(
		void *pdGridData,
		int numGids,
		int numLids,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int numWeights,
		float *objectWts,
		int *ierr
)
{

	*ierr = ZOLTAN_OK;
	PdGridData *gridData = (PdGridData *)pdGridData;
	int *gIds = gridData->myGlobalIDs.get();
	for(int i=0; i<gridData->numPoints; i++){
		zoltanGlobalIds[i] = gIds[i];
		zoltanLocalIds[i] = i;
	}
}

int zoltanQuery_dimension
(
		void *pdGridData,
		int *ierr
)
{
	*ierr = ZOLTAN_OK;
	PdGridData *gridData = (PdGridData *)pdGridData;
	return gridData->dimension;
}

/*
 * A ZOLTAN_GEOM_MULTI FN query function returns a vector of geometry values for a list of given objects.
 * The geometry vector is allocated by Zoltan to be of size num_obj * num_dim; its format is described below.
 * Function Type: 	ZOLTAN_GEOM_MULTI_FN_TYPE
 * Arguments:
    data 	Pointer to user-defined data.
    num_gid_entries 	The number of array entries used to describe a single global ID.  This value is the
		maximum value over all processors of the parameter NUM_GID_ENTRIES.
    num_lid_entries 	The number of array entries used to describe a single local ID.  This value is the
		maximum value over all processors of the parameter NUM_LID_ENTRIES. (It should be zero if local ids are not used.)
    num_obj 	The number of object IDs in arrays global_ids and local_ids.
    global_ids 	Array of global IDs of objects whose geometry values should be returned.
    local_ids 	Array of local IDs of objects whose geometry values should be returned. (Optional.)
    num_dim 	Number of coordinate entries per object (typically 1, 2, or 3).
    geom_vec 	Upon return, an array containing geometry values.
		For object i (specified by global_ids[i*num_gid_entries] and local_ids[i*num_lid_entries], i=0,1,...,num_obj-1),
		coordinate values should be stored in geom_vec[i*num_dim:(i+1)*num_dim-1].
    ierr 	Error code to be set by function.
 */
void zoltanQuery_gridData
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numPoints,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int dimension,
		double *zoltan_gridData,
		int *ierr
)
{


	PdGridData *gridData = (PdGridData *)pdGridData;

	/*
	 * In this app -- numGids and numLids should be "1"; also assert dimension is correct
	 */
	if ( (numGids != 1) || (numLids != 1) || (dimension != gridData->dimension)){
		*ierr = ZOLTAN_FATAL;
		return;
	}

	/*
	 * This is also bad!
	 */
	if ( numPoints != gridData->numPoints ){
		*ierr = ZOLTAN_FATAL;
		return;
	}

	*ierr = ZOLTAN_OK;
	double *x = gridData->myX.get();
	int c=0;
	for(int point=0;point<numPoints;point++){
		int localId = zoltanLocalIds[point];
		for(int d=0;d<dimension;d++,c++)
			zoltan_gridData[c] = x[localId*dimension+d];
	}
}

void zoltanQuery_pointSizeInBytes
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numPoints,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int *sizes,
		int *ierr
)
{

	PdGridData *gridData = (PdGridData *)pdGridData;
	/*
	 * In this app -- numGids and numLids should be "1"; also assert dimension is correct
	 */
	if ( (numGids != 1) || (numLids != 1) ){
		*ierr = ZOLTAN_FATAL;
		cout << "pointSize FATAL error" << endl;
		return;
	}

	/*
	 * What are we packing up?  For each point:
	 * 1)  coordinates:  size = dimension*sizeof(double)
	 * 1a) volume:       size = sizeof(double)
	 * 2)  numNeighbors: size = sizeof(int)
	 * 3)  neighbors:    size = numNeighbors*sizeof(int)
	 */
	*ierr = ZOLTAN_OK;
	int bytesPerDouble = sizeof(double);
	int bytesPerInt = sizeof(int);
	int dimension = gridData->dimension;


	//neighbor lists
	shared_ptr<int> neighborList = gridData->neighborhood;
	int *nPtr = neighborList.get();
	shared_ptr<int> neighborListPtr = gridData->neighborhoodPtr;
	int *ptrIntoList = neighborListPtr.get();


	// loop over incoming number of points
	for(int point=0;point<numPoints;point++){

		// get local id
		ZOLTAN_ID_TYPE localId = zoltanLocalIds[point];
		int intPtr = ptrIntoList[localId];

		// numNeighbors for this point
		int numNeigh = nPtr[intPtr];

		// number of bytes per point
		// coordinates
		int numBytesPerPoint = dimension*bytesPerDouble;
		// volume
		numBytesPerPoint += bytesPerDouble;
		// numNeighbors
		numBytesPerPoint += bytesPerInt;
		// neighbor list
		numBytesPerPoint += numNeigh*bytesPerInt;
		sizes[point] = numBytesPerPoint;
	}
}

void zoltanQuery_packPointsMultiFunction
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numExport,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int *dest,
		int *sizes,
		int *idx,
		char *buf,
		int *ierr
)

{

	PdGridData *gridData = (PdGridData *)pdGridData;
	/*
	 * In this app -- numGids and numLids should be "1";
	 */
	if ( (numGids != 1) || (numLids != 1) ){
		*ierr = ZOLTAN_FATAL;
		return;
	}


	*ierr = ZOLTAN_OK;

	// mark points that are for export
	char *exportFlagPtr = gridData->exportFlag.get();
	for(int p=0;p<gridData->numPoints;p++){
		exportFlagPtr[p]=0;
	}
	gridData->numExport=numExport;

//	std::cout << "zoltanQuery_packPointsMultiFunction " << std::endl; std::cout.flush();
//	std::cout << "\tNumber of points owned = " << gridData->numPoints
//			  << "\n\tNumber of points to be exported = " << numExport << std::endl; std::cout.flush();
	/*
	 * Pack Data in three steps:
	 * 1) coordinates and using one memcpy
	 * 2) volume using one memcpy
	 * 3) numNeigh and neighbors using one memcpy
	 */
	int dimension = gridData->dimension;
	double *X = gridData->myX.get();
	double *V = gridData->cellVolume.get();
	int *neighborList = gridData->neighborhood.get();
	int *neighborListPtr = gridData->neighborhoodPtr.get();

	// iterate over addresses and copy bytes
	int *idxEnd = idx+numExport;

	ZOLTAN_ID_PTR localIdsPtr = zoltanLocalIds;

	for(int *idxPtr = idx, *sizesPtr = sizes; idxPtr != idxEnd; idxPtr++, sizesPtr++, localIdsPtr++ ){

		ZOLTAN_ID_TYPE id = *localIdsPtr;

		// Mark this point as exported
		exportFlagPtr[id]=1;

		char *tmp = &buf[*idxPtr];
		// coordinates
		void *xPtr = (void*)(&X[dimension*id]);
		int numBytes = dimension*sizeof(double);
		memcpy((void*)tmp,xPtr,numBytes);

		// advance buffer pointer
		tmp += numBytes;

		// cell volume
		numBytes = sizeof(double);
		void *volPtr = (void*)(&V[id]);
		memcpy((void*)tmp,volPtr,numBytes);

		// advance buffer pointer
		tmp += numBytes;

		// num neighbors plus neighbor list
		// "i" points into the list at the correct location for this localId
		int i = neighborListPtr[id];
		int numNeigh = neighborList[i];

		numBytes = (1+numNeigh)*sizeof(int);
		memcpy((void*)tmp,(void*)(&neighborList[i]),numBytes);

	}


}



void zoltanQuery_unPackPointsMultiFunction
(
		void *pdGridData,
		int numGids,
		int numImport,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		int *sizes,
		int *idx,
		char *buf,
		int *ierr
)
{

	PdGridData *gridData = (PdGridData *)pdGridData;

	/*
	 * Mark this processor as being unpacked
	 */
	gridData->unPack = false;

	/*
	 * In this app -- numGids should be "1"; also assert dimension is correct
	 */
	if ( numGids != 1  ){
		*ierr = ZOLTAN_FATAL;
		return;
	}
	/*
	 * Reconstruct a new PdGridData
	 * 0) numPoints = numPoints - numExport + numImport
	 * 1) Create a new PdGridData
	 * 2) Loop over points in existing
	 * a) if not exported that pack into new PdGridData
	 * b) else ignore
	 * Points that were not originally exported will be at the top
	 * Imported points will be at the bottom
	 */
	int dimension = gridData->dimension;
	int newNumPoints = gridData->numPoints - gridData->numExport + numImport;
	PdGridData newGridData = PdQuickGrid::allocatePdGridData(newNumPoints,dimension);

//	std::cout << "zoltanQuery_unPackPointsMultiFunction " << std::endl; std::cout.flush();
//	std::cout << "\tNumber of points owned = " << gridData->numPoints
//			  << "\n\tNumber of points exported = " << gridData->numExport
//			  << "\n\tNumber of points to import "<< numImport << std::endl; std::cout.flush();

	/*
	 * The above allocation handles:
	 * 1) X
	 * 2) V
	 * 3) globalIds
	 * 4) neighborhoodPtr
	 */
	std::tr1::shared_ptr<double> newX = newGridData.myX;                         double *newXPtr   = newX.get();
	std::tr1::shared_ptr<double> newV = newGridData.cellVolume;                  double *newVPtr   = newV.get();
	std::tr1::shared_ptr<int> newGlobalIds = newGridData.myGlobalIDs;            int    *newIdsPtr = newGlobalIds.get();
	std::tr1::shared_ptr<int> newNeighborhoodPtr = newGridData.neighborhoodPtr;  int    *newNeighPtrPtr = newNeighborhoodPtr.get();

	std::tr1::shared_ptr<double> X = gridData->myX;                              double *xPtr   = X.get();
	std::tr1::shared_ptr<double> V = gridData->cellVolume;                       double *vPtr   = V.get();
	std::tr1::shared_ptr<int> globalIds = gridData->myGlobalIDs;                 int    *idsPtr = globalIds.get();
	std::tr1::shared_ptr<int> neighborhoodPtr = gridData->neighborhoodPtr;       int    *neighPtrPtr = neighborhoodPtr.get();
	std::tr1::shared_ptr<int> neighborhood = gridData->neighborhood;             int    *neighPtr = neighborhood.get();

	std::tr1::shared_ptr<char> exportFlag = gridData->exportFlag;                char *exportPtr = exportFlag.get();

	// Sum over points that stay on processor for determining size of new list
	int newSizeNeighborhoodList = 0;

	// Copy over points from old gridData that have not been exported
	for(int p=0;p<gridData->numPoints;p++, exportPtr++, neighPtrPtr++, vPtr++, idsPtr++){
		// this means we keep this point
		if(0==*exportPtr){

			// coordinates
			for(int d=0;d<dimension;d++)
				newXPtr[d] = xPtr[d];
			newXPtr+=dimension;

			// volume
			*newVPtr = *vPtr;
			newVPtr++;

			// global id
			*newIdsPtr = *idsPtr;
			newIdsPtr++;

			// new neighborhood pointer
			*newNeighPtrPtr = newSizeNeighborhoodList;
			newNeighPtrPtr++;

			// Accumulate length of new neighborhoodList
			int numNeigh = neighPtr[*neighPtrPtr];
			newSizeNeighborhoodList += (1+numNeigh);

		}
		xPtr+=dimension;

	}

	// Allocate neigbhorhood list
	// This function call determines the additional length of the neighborhood needed for the incoming points
	// Note that the initial value of "newSizeNeighborhoodList" is an input
	int allocateSizeNeighborhoodList = computeSizeNewNeighborhoodList(newSizeNeighborhoodList, numImport, idx, buf, dimension);
	std::tr1::shared_ptr<int> newNeighborhood(new int[allocateSizeNeighborhoodList],Deleter<int>());
	int *newNeighPtr = newNeighborhood.get();
	// Loop over points not exported and copy existing neighborhood into new neighborhood
	// Re-initialize
	exportPtr = exportFlag.get();
	neighPtrPtr = neighborhoodPtr.get();
	for(int p=0;p<gridData->numPoints;p++, exportPtr++, neighPtrPtr++){
		// this means we keep this point
		if(0==*exportPtr){

			// Copy existing neighborhood list into new list
			int ptr = *neighPtrPtr;
			int numNeigh = neighPtr[ptr];
			*newNeighPtr = numNeigh; newNeighPtr++;
			for(int n=0;n<numNeigh;n++, newNeighPtr++)
				*newNeighPtr = neighPtr[ptr+1+n];

		}

	}

	// After the above loops,
	//  ** newNeighPtr is positioned exactly at the start of memory for "imported" points (loop directly above)
	//  ** Same is true for newXPtr
	//  ** Same is true for newVPtr

	// iterate over addresses and copy bytes
	int *idxEnd = idx+numImport;
	ZOLTAN_ID_PTR globalIdsPtr = zoltanGlobalIds;

	for(int *idxPtr = idx, *sizesPtr = sizes; idxPtr != idxEnd; idxPtr++, sizesPtr++ ){

		int totalNumBytes = *sizesPtr;

		// Use this temporary pointer for copying data of single point
		char *tmp = &buf[*idxPtr];

		// coordinates
		int numBytes = dimension*sizeof(double);
		memcpy((void*)newXPtr,(void*)tmp,numBytes);
		// 1) advance buffer pointer and point coordinates as well
		// 2) decrement number of bytes
		tmp += numBytes;
		newXPtr+=dimension;
		totalNumBytes -= numBytes;

		// cell volume
		numBytes = sizeof(double);
		memcpy((void*)newVPtr,(void*)tmp,numBytes);
		// 1) advance buffer pointer and volume pointer
		// 2) decrement number of bytes
		tmp += numBytes;
		newVPtr++;
		/*
		 * This is the remaining number of bytes in buffer for point
		 * NOTE: totalNumBytes != (1+numNeigh)*sizeof(int) -- most of the
		 * time this is true but on some decompositions its not -- its padded
		 * with extra bytes that ends up cause the memcpy below to do a bad
		 * write
		 */
		totalNumBytes -= numBytes;

		// tmp buffer now points to start of neighorhood list; extract number of neighbors
		int numNeigh = *((int*)tmp);

		// get neighbor list
		numBytes = (1+numNeigh)*sizeof(int);
		memcpy((void*)newNeighPtr,(void*)tmp,numBytes);
		newNeighPtr+=(1+numNeigh);

		// Compute and save neighborhood ptr
		*newNeighPtrPtr = newSizeNeighborhoodList;
		newSizeNeighborhoodList += (1+numNeigh);
		newNeighPtrPtr++;

		// cell id
		ZOLTAN_ID_TYPE id = *globalIdsPtr;
		*newIdsPtr = id;
		globalIdsPtr++; newIdsPtr++;


	}

	// Replace old data with new data
	// dimension (does no change)
	// globalNumPoints (does not change)
	gridData->numPoints = newNumPoints;
	gridData->sizeNeighborhoodList = allocateSizeNeighborhoodList;
	// now set number of export points to 0
	gridData->numExport = 0;
	gridData->myGlobalIDs = newGlobalIds;
	gridData->myX = newX;
	gridData->cellVolume = newV;
	gridData->neighborhood = newNeighborhood;
	gridData->neighborhoodPtr = newNeighborhoodPtr;
	gridData->exportFlag = newGridData.exportFlag;
//	std::cout << "zoltanQuery_unPackPointsMultiFunction: Finish" << std::endl;
}



PdGridData& createAndAddNeighborhood(PdGridData& decomp, double horizon, bool withSelf){
	shared_ptr< std::set<int> > frameSet = constructParallelDecompositionFrameSet(decomp,horizon);
	struct Zoltan_Struct *zoltan = decomp.zoltanPtr.get();

//	std::cout << "createAndAddNeighborhood:A" << std::endl;
	int rank, numProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
//	std::cout << "createAndAddNeighborhood:Aa" << std::endl;
	/*
	 * This is the maximum possible length of communication -- every frame point goes to every processor
	 */
	Pd_shared_ptr_Array<int> procsArrayPoint(numProcs);
	Pd_shared_ptr_Array<int> sendProcsArray(numProcs*frameSet->size());
	Pd_shared_ptr_Array<int> pointLocalIdsArray(numProcs*frameSet->size());
	/*
	 * Total number of points to be sent
	 */
	int nSend(0), nReceive(0);

	/*
	 * Size of each point
	 */
	int dimension = decomp.dimension;

	/*
	 * Communication plan
	 */
	struct Zoltan_Comm_Obj *plan;
	int error;

//	std::cout << "rank, nSend, numProcs*frameSet->size() = " << rank << ": " << nSend << ", " << numProcs*frameSet->size() << std::endl;
	{
		int* sendProcsPtr = sendProcsArray.get();
		int* procsPointPtr = procsArrayPoint.get();
		int* localIdsPtr = pointLocalIdsArray.get();


		int numProcsPoint;
		set<int>::iterator myPointsIter = frameSet->begin();
		double *x = decomp.myX.get();

		for(;myPointsIter!=frameSet->end();myPointsIter++){

			/*
			 * This is a local id
			 */
			int id = *myPointsIter;
			const double *xP = x+dimension*id;
			PointBB bb(xP,horizon);
//			if(1==rank){
//				std::cout << "localId, x,y,z = " << id << ", " << *xP << ", " << *(xP+1) << ", " << *(xP+2) << std::endl;
//				std::cout << "xMin, xMax " << bb.get_xMin() << ", " << bb.get_xMax() << std::endl;
//				std::cout << "yMin, yMax " << bb.get_yMin() << ", " << bb.get_yMax() << std::endl;
//				std::cout << "zMin, zMax " << bb.get_zMin() << ", " << bb.get_zMax() << std::endl;
//			}


			Zoltan_LB_Box_Assign(zoltan,bb.get_xMin(),bb.get_yMin(),bb.get_zMin(),bb.get_xMax(),bb.get_yMax(),bb.get_zMax(),procsPointPtr,&numProcsPoint);
			/*
			 * Get local id and copy it nSend times into idsPtr
			 */
			int *ptrToSendProcs = sendProcsPtr+nSend;
			int *i = localIdsPtr+nSend;
//			cout << "Zoltan_LB_Box_Assign -- numProcsPoint = " << numProcsPoint << endl;
			for(int j=0;j<numProcsPoint;j++){
				/*
				 * Skip this point and processor if myRank = proc
				 * We do not send points to ourself
				 */
				if(rank==procsPointPtr[j]) continue;
				*i = id;
				/*
				 * Increment pointers
				 */
				*ptrToSendProcs = procsPointPtr[j];
				ptrToSendProcs++; i++;

				/*
				 * Increment total number of points to be sent
				 */
				nSend++;

			}
		}


		/*
		 * Create "communication" plan
		 */
		error = Zoltan_Comm_Create(&plan,nSend,sendProcsPtr,MPI_COMM_WORLD,PdQuickGrid::commCreate,&nReceive);
	}


	/*
	 * Calculate size of each point sent
	 */
	int nBytes(0);
	/*
	 * Global id
	 */
	nBytes += sizeof(int);
	/*
	 * Coordinates
	 */
	nBytes += dimension * sizeof(double);

	/*
	 * Buffer of data to be sent and received
	 */
	shared_ptr<char> sendBuffPtr(new char[nSend*nBytes],PdQuickGrid::Deleter<char>());
	shared_ptr<char> receBuffPtr(new char[nReceive*nBytes],PdQuickGrid::Deleter<char>());

	{
		/*
		 * Pack send buffer
		 */
		char *b = sendBuffPtr.get();
		int* idsPtr = pointLocalIdsArray.get();
		double *X = decomp.myX.get();
		int* myGIds = decomp.myGlobalIDs.get();

		for(int i=0;i<nSend;i++,b+=nBytes,idsPtr++){

			char *tmp = b;
			/*
			 * Copy global id
			 * NOTE: idsPtr is a "localId"
			 */
			int localId = *idsPtr;

			std::size_t numBytes = sizeof(int);
			void* gIdPtr = (void*)(myGIds+localId);
			memcpy((void*)tmp,gIdPtr,numBytes);
			tmp += numBytes;

			// coordinates
			numBytes = dimension * sizeof(double);
			void *xPtr = (void*)(X+dimension*localId);
			memcpy((void*)tmp,xPtr,numBytes);

		}
	}

	/*
	 * Send / Receive Data
	 */

	char *receiveBuff = receBuffPtr.get();
	Zoltan_Comm_Do(plan,PdQuickGrid::commDo,sendBuffPtr.get(),nBytes,receiveBuff);

	/*
	 * Now create neighborhood list
	 */
	int newNumPoints = decomp.numPoints + nReceive;

	/*
	 * First, create set of overlap points
	 */
	Pd_shared_ptr_Array<double> xOverlapArray(newNumPoints*dimension);
	Pd_shared_ptr_Array<int> gIdsOverlapArray(newNumPoints);

	{
		double *xOverlap = xOverlapArray.get();
		int *gIdsOverlap = gIdsOverlapArray.get();
		double *xOwned = decomp.myX.get();
		int *gIdsOwned = decomp.myGlobalIDs.get();
		int numPoints = decomp.numPoints;

		/*
		 * Copy existing coordinates and gIds
		 */
		for(int n=0;n<numPoints;n++,xOverlap+=3,xOwned+=3,gIdsOverlap++,gIdsOwned++){
			*gIdsOverlap = *gIdsOwned;
			*(xOverlap+0)=*(xOwned+0);
			*(xOverlap+1)=*(xOwned+1);
			*(xOverlap+2)=*(xOwned+2);
		}

		/*
		 * Now append incoming data
		 */
		char *b = receBuffPtr.get();
		for(int i=0;i<nReceive;i++,b+=nBytes,xOverlap+=3,gIdsOverlap++){
			char *tmp = b;

			/*
			 * gId
			 */
			std::size_t numBytes = sizeof(int);
			memcpy((void*)gIdsOverlap,(void*)tmp,numBytes);
			tmp += numBytes;

			// coordinates
			numBytes = dimension * sizeof(double);
			memcpy((void*)xOverlap,(void*)tmp,numBytes);

			/*
			 * DEBUG
			 */
//			if(rank==3){
//				std::cout << "new number of points = " << newNumPoints << std::endl;
//				std::cout << "receiving gId = " << *gIdsOverlap << ": x,y,z = " << *(xOverlap+0) << ", " << *(xOverlap+1) << ", " << *(xOverlap+2) << std::endl;
//			}
			/*
			 * END DEBUG
			 */

		}

	}

	/*
	 * create neighborhood and store on decomp
	 */
	vtkSmartPointer<vtkUnstructuredGrid> overlapGrid = PdVTK::getGrid(xOverlapArray.get_shared_ptr(),newNumPoints);
	NeighborhoodList list = PdNeighborhood::getNeighborhoodList(horizon,decomp.numPoints,decomp.myX,overlapGrid,withSelf);
	decomp.sizeNeighborhoodList = list.getSizeNeighborhoodList();
	decomp.neighborhood = list.getNeighborhood();
	decomp.neighborhoodPtr = list.getNeighborhoodPtr();
	// now set number of export points to 0
	decomp.numExport = 0;

	/*
	 * Now, we have to map the neighborhood list from "localIds" to "globalIds"
	 * Copy gIds from above gIdsOverlapArray into gidNeigh
	 */
	int numPoints = decomp.numPoints;
	int *gidNeigh = decomp.neighborhood.get();
	int *gIdsOverlap = gIdsOverlapArray.get();
	for(int p=0;p<numPoints;p++){
		int numNeigh = *gidNeigh; gidNeigh++;
		for(int n=0;n<numNeigh;n++,gidNeigh++){
			*gidNeigh = *(gIdsOverlap+*gidNeigh);
		}
	}


	Zoltan_Comm_Destroy(&plan);
	return decomp;
}


/*
 * This function is private to the library and should not be called otherwise
 * This function is called from the zoltan callback: zoltanQuery_unPackPointsMultiFuntion
 * Inputs to this function are those coming from the zoltan call back: zoltanQuery_unPackPointsMultiFuntion
 */

int computeSizeNewNeighborhoodList(int initialValue, int numImport, int *idx, char *buf, int dimension)
{

	int newSizeNeighborhoodList = initialValue;

	// Loop over sizes and compute new length of neighborhood list
	// iterate over addresses and copy bytes
	int *idxEnd = idx+numImport;

	for(int *idxPtr = idx; idxPtr != idxEnd; idxPtr++ ){

		// Use this temporary pointer for copying data of single point
		char *tmp = &buf[*idxPtr];

		// Need to position buffPtr at start of neighborhood for each point
		// increment pointer by coordinates, and volume

		// coordinates
		int numBytes = dimension*sizeof(double);
		// cell volume
		numBytes += sizeof(double);
		// move pointer
		tmp += numBytes;

		int numNeigh = *((int*)tmp);
		newSizeNeighborhoodList += (1+numNeigh);

	}

	return newSizeNeighborhoodList;

}



