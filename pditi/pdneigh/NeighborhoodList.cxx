/*
 * Neighborhood.cxx
 *
 *  Created on: Feb 25, 2011
 *      Author: jamitch
 */
#include "NeighborhoodList.h"
#include "Epetra_Comm.h"
#include "zoltan.h"
#include "PdVTK.h"
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include <stdexcept>
/*
 * @TODO
 * Need to remove this header dependency
 */
#include "PdNeighborhood.h"


namespace PDNEIGH {

/*
 * Prototype for private function
 */
const Epetra_BlockMap getOverlap(int ndf, int numShared, const int* shared, int numOwned,const int* owned, const Epetra_Comm& comm);

class PointBB {
public:
	PointBB (const double* x, double horizon):
	xMin(*(x+0)-horizon), xMax(*(x+0)+horizon),
	yMin(*(x+1)-horizon), yMax(*(x+1)+horizon),
	zMin(*(x+2)-horizon), zMax(*(x+2)+horizon) {}
	inline double get_xMin() const { return xMin; }
	inline double get_yMin() const { return yMin; }
	inline double get_zMin() const { return zMin; }
	inline double get_xMax() const { return xMax; }
	inline double get_yMax() const { return yMax; }
	inline double get_zMax() const { return zMax; }

private:
	double xMin, xMax, yMin, yMax, zMin, zMax;
};

double NeighborhoodList::get_horizon() const {
	return horizon;
}

int NeighborhoodList::get_num_owned_points() const {
	return num_owned_points;
}

shared_ptr<double> NeighborhoodList::get_owned_x() const {
	return owned_x;
}

shared_ptr<int> NeighborhoodList::get_neighborhood_ptr() const {
	return neighborhood_ptr;
}

shared_ptr<int> NeighborhoodList::get_neighborhood() const {
	return neighborhood;
}

const int* NeighborhoodList::get_neighborhood (int localId) const {
	int ptr = *(neighborhood_ptr.get()+localId);
	return neighborhood.get()+ptr;
}

int NeighborhoodList::get_size_neighborhood_list () const {
	return size_neighborhood_list;
}

int NeighborhoodList::get_num_neigh (int localId) const {
	return *(num_neighbors.get()+localId);
}

NeighborhoodList::NeighborhoodList
(
		struct Zoltan_Struct* zz,
		size_t numOwnedPoints,
		shared_ptr<int>& ownedGIDs,
		shared_ptr<double>& owned_coordinates,
		double horizon,
		shared_ptr<PdBondFilter::BondFilter> bondFilterPtr
)
:
		num_owned_points(numOwnedPoints),
		size_neighborhood_list(1),
		horizon(horizon),
		owned_gids(ownedGIDs),
		owned_x(owned_coordinates),
		neighborhood(new int[1],ArrayDeleter<int>()),
		neighborhood_ptr(new int[numOwnedPoints],ArrayDeleter<int>()),
		num_neighbors(new int[numOwnedPoints],ArrayDeleter<int>()),
		zoltan(zz),
		filter_ptr(bondFilterPtr)
{
	/*
	 * This single call performs the parallel search
	 */
	createAndAddNeighborhood();
}

pair<int, shared_ptr<int> > NeighborhoodList::getSharedGlobalIds() const {
	std::set<int> ownedIds(owned_gids.get(),owned_gids.get()+num_owned_points);
	std::set<int> shared;
	const int *neighPtr = neighborhood_ptr.get();
	const int *neigh = neighborhood.get();
	std::set<int>::const_iterator ownedIdsEnd = ownedIds.end();
	for(int p=0;p<num_owned_points;p++){
		int ptr = neighPtr[p];
		int numNeigh = neigh[ptr];
		for(int n=1;n<=numNeigh;n++){
			int id = neigh[ptr+n];
			/*
			 * look for id in owned points
			 */
			if(ownedIdsEnd == ownedIds.find(id)){
				/*
				 * add this point to shared
				 */
				shared.insert(id);
			}
		}
	}

	// Copy set into shared ptr
	shared_ptr<int> sharedGlobalIds(new int[shared.size()],ArrayDeleter<int>());
	int *sharedPtr = sharedGlobalIds.get();
	std::set<int>::iterator it;
	for ( it=shared.begin() ; it != shared.end(); it++, sharedPtr++ )
		*sharedPtr = *it;

	return pair<int, shared_ptr<int> >(shared.size(),sharedGlobalIds);
}


const Epetra_BlockMap NeighborhoodList::getOwnedMap(const Epetra_Comm& comm, int ndf) const {
	int numShared=0;
	const int *sharedPtr=NULL;
	int numOwned = num_owned_points;
	const int *ownedPtr = owned_gids.get();
	return getOverlap(ndf, numShared,sharedPtr,numOwned,ownedPtr,comm);
}

const Epetra_BlockMap NeighborhoodList::getOverlapMap(const Epetra_Comm& comm, int ndf) const {
	std::pair<int, std::tr1::shared_ptr<int> > sharedPair = this->getSharedGlobalIds();
	shared_ptr<int> sharedPtr = sharedPair.second;
	int numShared = sharedPair.first;
	const int *shared = sharedPtr.get();
	const int *owned = owned_gids.get();
	int numOwned = num_owned_points;
	return getOverlap(ndf,numShared,shared,numOwned,owned,comm);
}

void NeighborhoodList::createAndAddNeighborhood(){

	enum Pd_MPI_TAGS {commCreate=9,commDo=10};
	shared_ptr< std::set<int> > frameSet = constructParallelDecompositionFrameSet();


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
	 * dimension for each point must be '3'
	 */
	int dimension = 3;

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
		std::set<int>::iterator myPointsIter = frameSet->begin();
		double *x = owned_x.get();

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
		error = Zoltan_Comm_Create(&plan,nSend,sendProcsPtr,MPI_COMM_WORLD,commCreate,&nReceive);
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
	shared_ptr<char> sendBuffPtr(new char[nSend*nBytes],ArrayDeleter<char>());
	shared_ptr<char> receBuffPtr(new char[nReceive*nBytes],ArrayDeleter<char>());

	{
		/*
		 * Pack send buffer
		 */
		char *b = sendBuffPtr.get();
		int* idsPtr = pointLocalIdsArray.get();
		double *X = owned_x.get();
		int* myGIds = owned_gids.get();

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
	Zoltan_Comm_Do(plan,commDo,sendBuffPtr.get(),nBytes,receiveBuff);

	/*
	 * Now create neighborhood list
	 */
	int newNumPoints = num_owned_points + nReceive;

	/*
	 * First, create set of overlap points
	 */
	Pd_shared_ptr_Array<double> xOverlapArray(newNumPoints*dimension);
	Pd_shared_ptr_Array<int> gIdsOverlapArray(newNumPoints);

	{
		double *xOverlap = xOverlapArray.get();
		int *gIdsOverlap = gIdsOverlapArray.get();
		double *xOwned = owned_x.get();
		int *gIdsOwned = owned_gids.get();
		int numPoints = num_owned_points;

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
	 * create neighborhood
	 */
	buildNeighborhoodList(newNumPoints,xOverlapArray.get_shared_ptr());

	/*
	 * Now, we have to map the neighborhood list from "localIds" to "globalIds"
	 * Copy gIds from above gIdsOverlapArray into gidNeigh
	 * Also compute overall size of neighborhood list
	 */
	size_neighborhood_list=num_owned_points;
	int *gidNeigh = neighborhood.get();
	int *num_neigh = num_neighbors.get();
	int *gIdsOverlap = gIdsOverlapArray.get();
	for(int p=0;p<num_owned_points;p++,num_neigh++){
		*num_neigh = *gidNeigh; gidNeigh++;
		size_neighborhood_list += *num_neigh;
		for(int n=0;n<*num_neigh;n++,gidNeigh++){
			*gidNeigh = *(gIdsOverlap+*gidNeigh);
		}
	}


	Zoltan_Comm_Destroy(&plan);

}

shared_ptr< std::set<int> > NeighborhoodList::constructParallelDecompositionFrameSet() const {
	shared_ptr<double> xPtr = owned_x;
	shared_ptr<int> gIdsPtr = owned_gids;
	int numCells = num_owned_points;
	const PdNeighborhood::Coordinates c(xPtr,numCells);

	int numAxes=3;
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y, PdNeighborhood::Z};
	std::vector<std::tr1::shared_ptr<int> > sortedMaps(numAxes);
	for(int j=0;j<numAxes;j++){

		PdNeighborhood::Coordinates::SortComparator compare = c.getSortComparator(labels[j]);
		sortedMaps[j] = c.getIdentityMap();
		/*
		 * Sort points
		 */
		std::sort(sortedMaps[j].get(),sortedMaps[j].get()+numCells,compare);

	}

	/*
	 * Loop over axes and collect points at min and max ranges
	 * Add Points to frame set
	 */
	shared_ptr< std::set<int> >  frameSetPtr(new std::set<int>);
	PdNeighborhood::CoordinateLabel *label = labels;
	PdNeighborhood::CoordinateLabel *endLabel = labels+numAxes;
	for(;label!=endLabel;label++){

		{
			/*
			 * MINIMUM
			 * Find least upper bound of points for Min+horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			std::tr1::shared_ptr<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();
			/*
			 * First value in map corresponds with minimum value
			 */
			int iMIN = *map;
			double min = x[3*iMIN+axis];
			double value = min + horizon;
			const PdNeighborhood::Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const PdNeighborhood::Coordinates::SearchIterator end=start+numCells;
			PdNeighborhood::Coordinates::SearchIterator lub = std::upper_bound(start,end,value);
			frameSetPtr->insert(lub.mapStart(),lub.mapIterator());

		}

		{
			/*
			 * MAXIMUM
			 * Find greatest lower bound glb for Max-horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			std::tr1::shared_ptr<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();

			/*
			 * Last value in map corresponds with maximum value
			 */
			int iMAX = *(map+numCells-1);
			double max = x[3*iMAX+axis];
			double value = max - horizon;
			const PdNeighborhood::Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const PdNeighborhood::Coordinates::SearchIterator end=start+numCells;
			PdNeighborhood::Coordinates::SearchIterator glb = std::upper_bound(start,end,value);
			frameSetPtr->insert(glb.mapIterator(),glb.mapEnd());
		}
	}

	return frameSetPtr;

}

void NeighborhoodList::buildNeighborhoodList
(
		int numOverlapPoints,
		std::tr1::shared_ptr<double> xOverlapPtr
)
{
	/*
	 * Create KdTree
	 */
	vtkSmartPointer<vtkUnstructuredGrid> overlapGrid = PdVTK::getGrid(xOverlapPtr,numOverlapPoints);
	vtkKdTreePointLocator* kdTree = vtkKdTreePointLocator::New();
	kdTree->SetDataSet(overlapGrid);

	/*
	 * this is used by bond filters
	 */
	const double* xOverlap = xOverlapPtr.get();

	size_t sizeList = 0;
	size_t max=0;
	{
		/*
		 * Loop over owned points and determine number of points in horizon
		 */
		double *x = owned_x.get();
		const double *x_end = x+3*num_owned_points;
		size_t localId=0;
		for(;x!=x_end;x+=3, localId++){

			vtkIdList* kdTreeList = vtkIdList::New();
			/*
			 * Note that list returned includes this point *
			 */
			kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);

			if(0==kdTreeList->GetNumberOfIds()){
				/*
				 * Houston, we have a problem
				 */
				std::stringstream sstr;
				sstr << "\nERROR-->NeighborhoodList::buildNeighborhoodList(..)\n";
				sstr << "\tKdTree search failed to find any points in its neighborhood including itself!\n\tThis is probably a problem.\n";
				sstr << "\tLocal point id = " << localId << "\n"
					 << "\tSearch horizon = " << horizon << "\n"
					 << "\tx,y,z = " << *(x) << ", " << *(x+1) << ", " << *(x+2) << std::endl;
				std::string message=sstr.str();
				throw std::runtime_error(message);
			}

			size_t ptListSize = kdTreeList->GetNumberOfIds()+1;
			sizeList += ptListSize;

			/*
			 * Determine maximum possible number of neighbors over all points
			 */
			{
				size_t numIds = kdTreeList->GetNumberOfIds();
				if(numIds>max) max=numIds;
			}
			kdTreeList->Delete();
		}
	}
	/*
	 * Second pass to populate neighborhood list
	 */
	neighborhood_ptr = shared_ptr<int>(new int[num_owned_points],ArrayDeleter<int>());
	neighborhood=shared_ptr<int>(new int[sizeList],ArrayDeleter<int>());
	shared_ptr<bool> markForExclusion(new bool[max],ArrayDeleter<bool>());

	{
		/*
		 * Loop over owned points and determine number of points in horizon
		 */
		int *ptr = neighborhood_ptr.get();
		int neighPtr = 0;
		int *list = neighborhood.get();
		double *x = owned_x.get();
		for(int p=0;p<num_owned_points;p++,x+=3,ptr++){
			*ptr = neighPtr;
			vtkIdList* kdTreeList = vtkIdList::New();
			/*
			 * Note that list returned includes this point * but at start of list
			 */
			kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);
			bool *bondFlags = markForExclusion.get();
			filter_ptr->filterBonds(kdTreeList, x, p, xOverlap, bondFlags);

			/*
			 * Determine number of neighbors from flags
			 * Save start of list so that number of neighbors can be assigned after it
			 * has been calculated; Then increment pointer to first neighbor
			 */

			/*
			 * Save address for number of neighbors; will assign later
			 */
			int *numNeighPtr = list; list++;
			/*
			 * Loop over flags and save neighbors as appropriate; also accumulate number of neighbors
			 */
			size_t numNeigh=0;
			for(int n=0;n<kdTreeList->GetNumberOfIds();n++,bondFlags++){
				if(1==*bondFlags) continue;
				int uid = kdTreeList->GetId(n);
				*list = uid;
				list++;
				numNeigh++;
			}

			/*
			 * Now save number of neighbors
			 */
			*numNeighPtr = numNeigh;
			/*
			 * increment neighborhood pointer
			 */
			neighPtr += (numNeigh+1);
			/*
			 * Delete list
			 */
			kdTreeList->Delete();
		}
	}

	kdTree->Delete();

}


/*
 * PRIVATE FUNCTION
 */
const Epetra_BlockMap getOverlap(int ndf, int numShared, const int* shared, int numOwned, const int* owned, const Epetra_Comm& comm){

	int numPoints = numShared+numOwned;
	shared_ptr<int> ids(new int[numPoints],ArrayDeleter<int>());
	int *ptr = ids.get();

	for(int j=0;j<numOwned;j++,ptr++)
		*ptr=owned[j];

	for(int j=0;j<numShared;j++,ptr++)
		*ptr=shared[j];

	return Epetra_BlockMap(-1,numPoints, ids.get(),ndf, 0,comm);
}


}
