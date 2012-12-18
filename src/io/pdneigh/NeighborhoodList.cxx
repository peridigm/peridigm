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

#include "NeighborhoodList.h"

#include "Sortable.h"
#include "Array.h"
#include "Vector3D.h"

#include "zoltan.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"

#include "Peridigm_VTKSearchTree.hpp"
#include "Peridigm_JAMSearchTree.hpp"
#include "Peridigm_ZoltanSearchTree.hpp"

#include <stdexcept>

namespace PDNEIGH {


using UTILITIES::Sortable;
using UTILITIES::CartesianComponent;
using UTILITIES::Array;
using UTILITIES::PointCenteredBoundingBox;

/*
 * Prototype for private function
 */
shared_ptr<Epetra_BlockMap> getOverlap(int ndf, int numShared, const int* shared, int numOwned, const int* owned, const Epetra_Comm& comm);

shared_ptr<const Epetra_Comm> NeighborhoodList::get_Epetra_Comm() const {
	return epetraComm;
}

shared_ptr<Epetra_Distributor> NeighborhoodList::create_Epetra_Distributor() const {
	return shared_ptr<Epetra_Distributor>(epetraComm->CreateDistributor());
}

double NeighborhoodList::get_horizon() const {
	return horizon;
}

size_t NeighborhoodList::get_num_owned_points() const {
	return num_owned_points;
}

size_t NeighborhoodList::get_num_shared_points() const {
	return sharedGIDs.get_size();
}

shared_ptr<double> NeighborhoodList::get_owned_x() const {
	return owned_x;
}

shared_ptr<int> NeighborhoodList::get_neighborhood_ptr() const {
	return neighborhood_ptr.get_shared_ptr();
}

shared_ptr<int> NeighborhoodList::get_neighborhood() const {
	return neighborhood.get_shared_ptr();
}

shared_ptr<int> NeighborhoodList::get_local_neighborhood() const {
	return local_neighborhood.get_shared_ptr();
}

shared_ptr<int> NeighborhoodList::get_owned_gids() const {
	return owned_gids;
}

shared_ptr<int> NeighborhoodList::get_shared_gids() const {
	return sharedGIDs.get_shared_ptr();
}

const int* NeighborhoodList::get_neighborhood (int localId) const {
	int ptr = *(neighborhood_ptr.get()+localId);
	return neighborhood.get()+ptr;
}

const int* NeighborhoodList::get_local_neighborhood (int localId) const {
	int ptr = *(neighborhood_ptr.get()+localId);
	return local_neighborhood.get()+ptr;
}


int NeighborhoodList::get_size_neighborhood_list () const {
	return size_neighborhood_list;
}

int NeighborhoodList::get_num_neigh (int localId) const {
	return *(num_neighbors.get()+localId);
}

NeighborhoodList::NeighborhoodList
(
		shared_ptr<const Epetra_Comm> comm,
		struct Zoltan_Struct* zz,
		size_t numOwnedPoints,
		shared_ptr<int>& ownedGIDs,
		shared_ptr<double>& owned_coordinates,
		double horizon,
		shared_ptr<PdBondFilter::BondFilter> bondFilterPtr
)
:
		epetraComm(comm),
		epetra_block_maps(),
		num_owned_points(numOwnedPoints),
		size_neighborhood_list(1),
		horizon(horizon),
		owned_gids(ownedGIDs),
		owned_x(owned_coordinates),
		neighborhood(),
		local_neighborhood(),
		neighborhood_ptr(num_owned_points),
		num_neighbors(num_owned_points),
		sharedGIDs(),
		zoltan(zz),
		filter_ptr(bondFilterPtr)
{
	/*
	 * This single call performs the parallel search
	 */
	createAndAddNeighborhood();

	/*
	 * Create list of shared GIDs; these are IDs that
	 * are appended to 'ownedIDs'; together, ownedIDs
	 * and sharedGIDs form the overlapVectors
	 */
	sharedGIDs = createSharedGlobalIds();

	/*
	 * Create common maps
	 */
	shared_ptr<Epetra_BlockMap> owned1DMap = create_Epetra_BlockMap(Epetra_MapTag(OWNED,1));
	epetra_block_maps[Epetra_MapTag(OWNED,1)] = owned1DMap; 
	shared_ptr<Epetra_BlockMap> owned3DMap = create_Epetra_BlockMap(Epetra_MapTag(OWNED,3));
	epetra_block_maps[Epetra_MapTag(OWNED,3)] = owned3DMap; 
	shared_ptr<Epetra_BlockMap> overlap1DMap = create_Epetra_BlockMap(Epetra_MapTag(OVERLAP,1));
	epetra_block_maps[Epetra_MapTag(OVERLAP,1)] = overlap1DMap;
	shared_ptr<Epetra_BlockMap> overlap3DMap = create_Epetra_BlockMap(Epetra_MapTag(OVERLAP,3));
	epetra_block_maps[Epetra_MapTag(OVERLAP,3)] = overlap3DMap;

	/*
	 * Compute local neighborhood list and
	 */
	local_neighborhood = createLocalNeighborList(*getOverlapMap(1));

}


/**
 * This function creates a 'new neighborhood' based upon an existing
 * set of owned points but with a 'newHorizon'.  This is efficient
 * in that both classes share the same underlying coordinates and
 * global ids and 'zoltan' struct.
 *
 * NOTE: THIS PERFORMS A NEW PARALLEL SEARCH
 * NOTE: Zoltan ptr stored must be valid
 */
NeighborhoodList NeighborhoodList::cloneAndShare(double newHorizon,bool withSelf) {
	NeighborhoodList newList(epetraComm,zoltan,num_owned_points,owned_gids,owned_x,newHorizon,filter_ptr->clone(withSelf));
	return newList;
}

Array<int> NeighborhoodList::createLocalNeighborList(const Epetra_BlockMap& overlapMap){
	Array<int> localNeighborList(size_neighborhood_list);
	int *localNeig = localNeighborList.get();
	int *neighPtr = neighborhood_ptr.get();
	int *neigh = neighborhood.get();
	for(size_t p=0;p<num_owned_points;p++){
		int ptr = neighPtr[p];
		int numNeigh = neigh[ptr];
		localNeig[ptr]=numNeigh;
		for(int n=1;n<=numNeigh;n++){
			int gid = neigh[ptr+n];
			int localId = overlapMap.LID(gid);
			localNeig[ptr+n] = localId;
		}
	}
	return localNeighborList;
}


Array<int> NeighborhoodList::createSharedGlobalIds() const {
	std::set<int> ownedIds(owned_gids.get(),owned_gids.get()+num_owned_points);
	std::set<int> shared;
	const int *neighPtr = neighborhood_ptr.get();
	const int *neigh = neighborhood.get();
	std::set<int>::const_iterator ownedIdsEnd = ownedIds.end();
	for(size_t p=0;p<num_owned_points;p++){
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
	Array<int> sharedGlobalIds(shared.size());
	int *sharedPtr = sharedGlobalIds.get();
	std::set<int>::iterator it;
	for ( it=shared.begin() ; it != shared.end(); it++, sharedPtr++ )
		*sharedPtr = *it;

	return sharedGlobalIds;
}


shared_ptr<Epetra_BlockMap> NeighborhoodList::create_Epetra_BlockMap(Epetra_MapTag key) {
	Epetra_MapType t = key.type;
	shared_ptr<Epetra_BlockMap> m;
	switch(t){
	case OWNED:
		m = getOwnedMap(key.ndf);
		break;
	case OVERLAP:
		m = getOverlapMap(key.ndf);
		break;
	}
	return m;
}


shared_ptr<Epetra_BlockMap> NeighborhoodList::getOwnedMap(int ndf) const {
	std::map<Epetra_MapTag, shared_ptr<Epetra_BlockMap>, MapComparator >::const_iterator i=epetra_block_maps.find(Epetra_MapTag(OWNED,ndf));
	if (epetra_block_maps.end() != i)
		return i->second;
	/*
	 * No map found, so we must create one
	 */
	int numShared=0;
	const int *sharedPtr=NULL;
	int numOwned = num_owned_points;
	const int *ownedPtr = owned_gids.get();
	return getOverlap(ndf, numShared,sharedPtr,numOwned,ownedPtr,*epetraComm);
}


shared_ptr<Epetra_BlockMap> NeighborhoodList::getOverlapMap(int ndf) const {
	std::map<Epetra_MapTag, shared_ptr<Epetra_BlockMap>, MapComparator >::const_iterator i=epetra_block_maps.find(Epetra_MapTag(OVERLAP,ndf));
	if (epetra_block_maps.end() != i)
		return i->second;
	/*
	 * No map found, so we must create one
	 */
	Array<int> sharedPtr = sharedGIDs;
	int numShared = sharedPtr.get_size();
	const int *shared = sharedPtr.get();
	const int *owned = owned_gids.get();
	int numOwned = num_owned_points;
	return getOverlap(ndf,numShared,shared,numOwned,owned,*epetraComm);
}

/*
 * PRIVATE FUNCTION
 */
shared_ptr<Epetra_BlockMap> getOverlap(int ndf, int numShared, const int* shared, int numOwned, const int* owned, const Epetra_Comm& comm){

	int numPoints = numShared+numOwned;
	Array<int> ids(numPoints);
	int *ptr = ids.get();

	for(int j=0;j<numOwned;j++,ptr++)
		*ptr=owned[j];

	for(int j=0;j<numShared;j++,ptr++)
		*ptr=shared[j];

	return shared_ptr<Epetra_BlockMap>(new Epetra_BlockMap(-1,numPoints, ids.get(),ndf, 0,comm));
}

void NeighborhoodList::createAndAddNeighborhood(){

	enum {COMM_CREATE=9,COMM_DO=10};

//	shared_ptr< std::set<int> > frameSet = constructParallelDecompositionFrameSet();
	shared_ptr< std::set<int> > frameSet = UTILITIES::constructFrameSet(num_owned_points,owned_x,horizon);


//	std::cout << "createAndAddNeighborhood:A" << std::endl;
	int rank, numProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
//	std::cout << "createAndAddNeighborhood:Aa" << std::endl;
	/*
	 * This is the maximum possible length of communication -- every frame point goes to every processor
	 */
	Array<int> procsArrayPoint(numProcs);
	Array<int> sendProcsArray(numProcs*frameSet->size());
	Array<int> pointLocalIdsArray(numProcs*frameSet->size());


	/*
	 * dimension for each point must be '3'
	 */
	int dimension = 3;

	/*
	 * Communication plan
	 */
	struct Zoltan_Comm_Obj *plan;

	/*
	 * Total number of points to be sent and received
	 */
	int nSend(0), nReceive(0);
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
			PointCenteredBoundingBox bb(xP,horizon);
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
		error = Zoltan_Comm_Create(&plan,nSend,sendProcsPtr,MPI_COMM_WORLD,COMM_CREATE,&nReceive);
        if(error)
          throw std::runtime_error("****Error in NeighborhoodList::createAndAddNeighborhood(), Zoltan_Comm_Create() returned a nonzero error code.");
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
	Zoltan_Comm_Do(plan,COMM_DO,sendBuffPtr.get(),nBytes,receiveBuff);

	/*
	 * Now create neighborhood list
	 */
	int newNumPoints = num_owned_points + nReceive;

	/*
	 * First, create set of overlap points
	 */
	Array<double> xOverlapArray(newNumPoints*dimension);
	Array<int> gIdsOverlapArray(newNumPoints);

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
	for(size_t p=0;p<num_owned_points;p++,num_neigh++){
		*num_neigh = *gidNeigh; gidNeigh++;
		size_neighborhood_list += *num_neigh;
		for(int n=0;n<*num_neigh;n++,gidNeigh++){
			*gidNeigh = *(gIdsOverlap+*gidNeigh);
		}
	}

	Zoltan_Comm_Destroy(&plan);

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
    PeridigmNS::SearchTree* searchTree = new PeridigmNS::VTKSearchTree(numOverlapPoints, xOverlapPtr.get());
    //PeridigmNS::SearchTree* searchTree = new PeridigmNS::JAMSearchTree(numOverlapPoints, xOverlapPtr.get());
    //PeridigmNS::SearchTree* searchTree = new PeridigmNS::ZoltanSearchTree(numOverlapPoints, xOverlapPtr.get());

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

            std::vector<int> treeList;
			/*
			 * Note that list returned includes this point *
			 */
			searchTree->FindPointsWithinRadius(x, horizon, treeList);

			if(0==treeList.size()){
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

			size_t ptListSize = treeList.size()+1;
			sizeList += ptListSize;

			/*
			 * Determine maximum possible number of neighbors over all points
			 */
			{
                size_t numIds = treeList.size();
				if(numIds>max) max=numIds;
			}
		}
	}
	/*
	 * Second pass to populate neighborhood list
	 */
	neighborhood_ptr = Array<int>(num_owned_points);
	neighborhood     = Array<int>(sizeList);
	Array<bool> markForExclusion(max);

	{
		/*
		 * Loop over owned points and determine number of points in horizon
		 */
		int *ptr = neighborhood_ptr.get();
		int neighPtr = 0;
		int *list = neighborhood.get();
		double *x = owned_x.get();
		for(size_t p=0;p<num_owned_points;p++,x+=3,ptr++){
			*ptr = neighPtr;
            std::vector<int> treeList;
			/*
			 * Note that list returned includes this point * but at start of list
			 */
			searchTree->FindPointsWithinRadius(x, horizon, treeList);

			//sort(treeList.begin(), treeList.end());

			bool *bondFlags = markForExclusion.get();
			filter_ptr->filterBonds(treeList, x, p, xOverlap, bondFlags);

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
			for(unsigned int n=0;n<treeList.size();n++,bondFlags++){
				if(1==*bondFlags) continue;
				int uid = treeList[n];
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
		}
	}

	delete searchTree;
}

}
