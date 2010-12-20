/*
 * PdNeighborhood.cxx
 *
 *  Created on: Mar 10, 2010
 *      Author: jamitch
 */

#include "PdNeighborhood.h"
#include "PdQuickGrid.h"
#include "PdVTK.h"
#include "PdBondFilter.h"
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include <stdexcept>

using std::pair;

namespace PdNeighborhood {

/*
 * Private function
 */
Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CoordinateLabel axis,
		double axisMinimumValue,
		double horizon,
		const Coordinates& c,
		std::tr1::shared_ptr<int>& sortedMap
);

/*
 * Private function
 */
Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CoordinateLabel axis,
		double axisMaximumValue,
		double horizon,
		const Coordinates& c,
		std::tr1::shared_ptr<int>& sortedMap
);

/*
 * Private function
 */
std::tr1::shared_ptr<int> getSortedMap(const Coordinates& c, CoordinateLabel axis);


/**
 * This interface assumes data layout like PdQuickGrid
 * @param neighborsPtr -- pointer into start of
 * neighborhood for localId
 * @param neighbors -- first entry must be "numNeigh"
 */
NeighborhoodList::NeighborhoodList
(
		int numPoints,
		std::tr1::shared_ptr<int> neighborsPtr,
		std::tr1::shared_ptr<int> neighbors
)
:
numPoints(numPoints),
sizeNeighborhoodList(numPoints),
numNeighbors(),
neighborsPtr(neighborsPtr),
neighbors(neighbors)
{
	/*
	 * construct number of neighbors; accumulate sizeNeighborhoodList
	 */

	numNeighbors = std::tr1::shared_ptr<int>(new int[numPoints],PdQuickGrid::Deleter<int>());
	int *nNeigh = numNeighbors.get();
	int *ptr = neighborsPtr.get();
	for(int i=0;i<numPoints;i++,nNeigh++,ptr++){
		int p = *ptr;
		*nNeigh = *(neighbors.get()+p);
		sizeNeighborhoodList += *nNeigh;
	}
}


std::tr1::shared_ptr<int> getSortedMap(const Coordinates& c, CoordinateLabel axis){
	int numPoints = c.getNumPoints();
	Coordinates::SortComparator compare = c.getSortComparator(axis);
	std::tr1::shared_ptr<int> sortedMap = c.getIdentityMap();
	/*
	 * Sort points
	 */
	std::sort(sortedMap.get(),sortedMap.get()+numPoints,compare);
	return sortedMap;
}

/**
 * Mostly for rectangular type meshes; Returns points in plane perpendicular to "axis";
 * All points within a distance of horizon of the axis minimum are returned
 * @param axis -- X || Y || Z
 * @param horizon
 */
Pd_shared_ptr_Array<int> getPointsAxisAlignedMinimum
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon
)
{

	const Coordinates c(xPtr,numPoints);
//	Coordinates::SortComparator compare = c.getSortComparator(axis);
//	std::tr1::shared_ptr<int> sortedMap = c.getIdentityMap();
//	/*
//	 * Sort points
//	 */
//	std::sort(sortedMap.get(),sortedMap.get()+numPoints,compare);

	std::tr1::shared_ptr<int> sortedMap = getSortedMap(c,axis);

	/*
	 * Get sorted points
	 */
	int *mapX = sortedMap.get();
	/*
	 * First value in map corresponds with minimum value of coordinate "axis"
	 */
	int iMIN = *mapX;
	double xMin = *(xPtr.get()+3*iMIN+axis);

	/*
	 * Look for points from "xMin" to "xMin+horizon"
	 */
//	double value = xMin + horizon;
//	const Coordinates::SearchIterator start=c.begin(axis,sortedMap);
//	const Coordinates::SearchIterator end=start+numPoints;
//	Coordinates::SearchIterator upper = std::upper_bound(start,end,value);
//	int numPointsInSet = upper.numPointsFromStart();
//	Pd_shared_ptr_Array<int> pointIds(numPointsInSet);
//	int *ids = pointIds.get();
//	for(const int *i=upper.mapStart();i!=upper.mapIterator();i++,ids++)
//		*ids = *i;
//
//	return pointIds;
	return getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,xMin,horizon,c,sortedMap);
}

Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon,
		double axisMinimumValue
)
{
	const Coordinates c(xPtr,numPoints);
	std::tr1::shared_ptr<int> sortedMap = getSortedMap(c,axis);
	return getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,axisMinimumValue,horizon,c,sortedMap);
}

Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CoordinateLabel axis,
		double axisMinimumValue,
		double horizon,
		const Coordinates& c,
		std::tr1::shared_ptr<int>& sortedMap
){
	int numPoints = c.getNumPoints();
	/*
	 * Look for points from "xMin" to "xMin+horizon"
	 */
	double value = axisMinimumValue + horizon;
	const Coordinates::SearchIterator start=c.begin(axis,sortedMap);
	const Coordinates::SearchIterator end=start+numPoints;
	Coordinates::SearchIterator upper = std::upper_bound(start,end,value);
	int numPointsInSet = upper.numPointsFromStart();
	Pd_shared_ptr_Array<int> pointIds(numPointsInSet);
	int *ids = pointIds.get();
	for(const int *i=upper.mapStart();i!=upper.mapIterator();i++,ids++)
		*ids = *i;

	return pointIds;
}

Pd_shared_ptr_Array<int> getPointsAxisAlignedMaximum
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon
)
{

	const Coordinates c(xPtr,numPoints);
//	Coordinates::SortComparator compare = c.getSortComparator(axis);
//	std::tr1::shared_ptr<int> sortedMap = c.getIdentityMap();
//	/*
//	 * Sort points
//	 */
//	std::sort(sortedMap.get(),sortedMap.get()+numPoints,compare);

	std::tr1::shared_ptr<int> sortedMap = getSortedMap(c,axis);
	/*
	 * Get sorted points
	 */
	int *mapX = sortedMap.get();

	/*
	 * Last value in map corresponds with maximum value
	 */
	int iMAX = *(mapX+numPoints-1);
	double max = *(xPtr.get()+3*iMAX+axis);

	/*
	 * Look for points from "max-horizon" to "max"
	 */
//	double value = max - horizon;
//	const Coordinates::SearchIterator start=c.begin(axis,sortedMap);
//	const Coordinates::SearchIterator end=start+numPoints;
//	Coordinates::SearchIterator upper = std::upper_bound(start,end,value);
//	int numPointsInSet = upper.numPointsToEnd();
//	Pd_shared_ptr_Array<int> pointIds(numPointsInSet);
//	int *ids = pointIds.get();
//	for(const int *i=upper.mapIterator();i!=upper.mapEnd();i++,ids++)
//		*ids = *i;
//	return pointIds;

	return getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,max,horizon,c,sortedMap);
}

Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon,
		double axisMaximumValue
)
{
	const Coordinates c(xPtr,numPoints);
	std::tr1::shared_ptr<int> sortedMap = getSortedMap(c,axis);
	return getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,axisMaximumValue,horizon,c,sortedMap);
}

Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CoordinateLabel axis,
		double axisMaximumValue,
		double horizon,
		const Coordinates& c,
		std::tr1::shared_ptr<int>& sortedMap
){

	int numPoints = c.getNumPoints();
	/*
	 * Look for points from "max-horizon" to "max"
	 */
	double value = axisMaximumValue - horizon;
	const Coordinates::SearchIterator start=c.begin(axis,sortedMap);
	const Coordinates::SearchIterator end=start+numPoints;
	Coordinates::SearchIterator upper = std::upper_bound(start,end,value);
	int numPointsInSet = upper.numPointsToEnd();
	Pd_shared_ptr_Array<int> pointIds(numPointsInSet);
	int *ids = pointIds.get();
	for(const int *i=upper.mapIterator();i!=upper.mapEnd();i++,ids++)
		*ids = *i;

	return pointIds;

}

NeighborhoodList getNeighborhoodList
(
		double horizon,
		int numOwnedPoints,
		int numOverlapPoints,
		std::tr1::shared_ptr<double> xOwnedPtr,
		std::tr1::shared_ptr<double> xOverlapPtr,
		PdBondFilter::BondFilter& bondFilter
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

	std::tr1::shared_ptr<int> neighborsPtr(new int[numOwnedPoints],PdQuickGrid::Deleter<int>());
	size_t sizeList = 0;
	size_t max=0;
	{
		/*
		 * Loop over owned points and determine number of points in horizon
		 */
		int *ptr = neighborsPtr.get();
		const int* const end = neighborsPtr.get()+numOwnedPoints;
		double *x = xOwnedPtr.get();
		int localId=0;
		for(;ptr!=end;ptr++,x+=3, localId++){
			/*
			 * pointer to start of neighborhood list
			 */
			*ptr = sizeList;

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
				sstr << "\nERROR-->PdNeighborhood.cxx::getNeighborhoodList(..)\n";
				sstr << "\tKdTree search failed to find any points in its neighborhood including itself!\n\tThis is probably a problem.\n";
				sstr << "\tLocal point id = " << localId << "\n"
					 << "\tSearch horizon = " << horizon << "\n"
					 << "\tx,y,z = " << *(x) << ", " << *(x+1) << ", " << *(x+2) << std::endl;
				std::string message=sstr.str();
				throw std::runtime_error(message);
			}

			size_t ptListSize = bondFilter.filterListSize(kdTreeList,x,xOverlap);
			if(ptListSize>max) max=ptListSize;
			sizeList += ptListSize;
			kdTreeList->Delete();
		}
	}
	/*
	 * Second pass to populate neighborhood list
	 */
	std::tr1::shared_ptr<int> neighborsList(new int[sizeList],PdQuickGrid::Deleter<int>());
	std::tr1::shared_ptr<bool> markForExclusion(new bool[max],PdQuickGrid::Deleter<bool>());

	{
		/*
		 * Loop over owned points and determine number of points in horizon
		 */
		int *list = neighborsList.get();
		double *x = xOwnedPtr.get();
		for(int p=0;p<numOwnedPoints;p++,x+=3){

			vtkIdList* kdTreeList = vtkIdList::New();
			/*
			 * Note that list returned includes this point * but at start of list
			 */
			kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);
			bool *bondFlags = markForExclusion.get();
			bondFilter.filterBonds(kdTreeList, x, p, xOverlap, bondFlags);

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
			 * Delete list
			 */
			kdTreeList->Delete();
		}
	}

	kdTree->Delete();
	return NeighborhoodList(numOwnedPoints,neighborsPtr,neighborsList);
}

}
