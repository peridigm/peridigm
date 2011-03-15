/*
 * Sortable.cxx
 *
 *  Created on: Mar 15, 2011
 *      Author: jamitch
 */
#include "Sortable.h"
#include "Array.h"
#include <algorithm>

namespace UTILITIES {
/*
 * THIS IS A PRIVATE FUNCTION
 */
Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CartesianComponent axis,
		double axisMinimumValue,
		double horizon,
		const Sortable& c,
		Array<int>& sortedMap
);

/*
 * THIS IS A PRIVATE FUNCTION
 */
Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CartesianComponent axis,
		double axisMaximumValue,
		double horizon,
		const Sortable& c,
		Array<int>& sortedMap
);

Array<int> getSortedMap(const Sortable& c, CartesianComponent axis){
	size_t numPoints = c.get_num_points();
	Sortable::Comparator compare = c.getComparator(axis);
	Array<int> sortedMap = c.getIdentityMap();
	/*
	 * Sort points
	 */
	std::sort(sortedMap.get(),sortedMap.get()+numPoints,compare);
	return sortedMap;
}

Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CartesianComponent axis,
		shared_ptr<double> xPtr,
		size_t numPoints,
		double horizon,
		double axisMinimumValue
)
{
	const Sortable c(numPoints,xPtr);
	Array<int> sortedMap = getSortedMap(c,axis);
	return getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,axisMinimumValue,horizon,c,sortedMap);
}

Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CartesianComponent axis,
		double axisMinimumValue,
		double horizon,
		const Sortable& c,
		Array<int>& sortedMap
){
	int numPoints = c.get_num_points();
	/*
	 * Look for points from "xMin" to "xMin+horizon"
	 */
	double value = axisMinimumValue + horizon;
	const Sortable::SearchIterator start=c.begin(axis,sortedMap.get_shared_ptr());
	const Sortable::SearchIterator end=start+numPoints;
	Sortable::SearchIterator upper = std::upper_bound(start,end,value);
	int numPointsInSet = upper.numPointsFromStart();
	Array<int> pointIds(numPointsInSet);
	int *ids = pointIds.get();
	for(const int *i=upper.mapStart();i!=upper.mapIterator();i++,ids++)
		*ids = *i;

	return pointIds;
}


/**
 * Mostly for rectangular type meshes; Returns points in plane perpendicular to "axis";
 * All points within a distance of horizon of the axis minimum are returned
 * @param axis -- X || Y || Z
 * @param horizon
 */
Array<int> getPointsAxisAlignedMinimum
(
		CartesianComponent axis,
		shared_ptr<double> xPtr,
		size_t numPoints,
		double horizon
)
{

	const Sortable c(numPoints,xPtr);
	Array<int> sortedMap = getSortedMap(c,axis);

	/*
	 * Get sorted points
	 */
	int *mapX = sortedMap.get();
	/*
	 * First value in map corresponds with minimum value of coordinate "axis"
	 */
	int iMIN = *mapX;
	double xMin = *(xPtr.get()+3*iMIN+axis);

	return getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,xMin,horizon,c,sortedMap);
}


Array<int> getPointsAxisAlignedMaximum
(
		CartesianComponent axis,
		shared_ptr<double> xPtr,
		size_t numPoints,
		double horizon
)
{

	const Sortable c(numPoints,xPtr);

	Array<int> sortedMap = getSortedMap(c,axis);
	/*
	 * Get sorted points
	 */
	int *mapX = sortedMap.get();

	/*
	 * Last value in map corresponds with maximum value
	 */
	int iMAX = *(mapX+numPoints-1);
	double max = *(xPtr.get()+3*iMAX+axis);

	return getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,max,horizon,c,sortedMap);
}

Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CartesianComponent axis,
		shared_ptr<double> xPtr,
		size_t numPoints,
		double horizon,
		double axisMaximumValue
)
{
	const Sortable c(numPoints,xPtr);
	Array<int> sortedMap = getSortedMap(c,axis);
	return getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,axisMaximumValue,horizon,c,sortedMap);
}

Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CartesianComponent axis,
		double axisMaximumValue,
		double horizon,
		const Sortable& c,
		Array<int>& sortedMap
){

	int numPoints = c.get_num_points();
	/*
	 * Look for points from "max-horizon" to "max"
	 */
	double value = axisMaximumValue - horizon;
	const Sortable::SearchIterator start=c.begin(axis,sortedMap.get_shared_ptr());
	const Sortable::SearchIterator end=start+numPoints;
	Sortable::SearchIterator upper = std::upper_bound(start,end,value);
	int numPointsInSet = upper.numPointsToEnd();
	Array<int> pointIds(numPointsInSet);
	int *ids = pointIds.get();
	for(const int *i=upper.mapIterator();i!=upper.mapEnd();i++,ids++)
		*ids = *i;

	return pointIds;

}


} /* UTITILITIES */
