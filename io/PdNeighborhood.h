/*
 * PdNeighborhood.h
 *
 *  Created on: Mar 10, 2010
 *      Author: jamitch
 */

#ifndef PDNEIGHBORHOOD_H_
#define PDNEIGHBORHOOD_H_
#include<tr1/memory>
#include<iostream>
#include "PdVTK.h"
#include "PdQuickGrid.h"

namespace PdNeighborhood {
enum CoordinateLabel {X=0,Y,Z};
class NeighborhoodList {
private:
	int numPoints, sizeNeighborhoodList;
	std::tr1::shared_ptr<int> numNeighbors, neighborsPtr, neighbors;
public:

	/*
	 * INPUT to this construct must be identical to
	 * that associated with "PdGridData" and/or that
	 * coming out of "PdQuickGrid"
	 */
	NeighborhoodList
	(
			int numPoints,
			std::tr1::shared_ptr<int> neighborsPtr,
			std::tr1::shared_ptr<int> neighbors
	);

	int getNumPoints() const { return numPoints; }

	std::tr1::shared_ptr<int> getNeighborhoodPtr() const {
		return neighborsPtr;
	}

	std::tr1::shared_ptr<int> getNeighborhood() const {
		return neighbors;
	}

	const int* getNeighborhood (int localId) const {
		int ptr = *(neighborsPtr.get()+localId);
		return neighbors.get()+ptr;
	}

	int getNumNeigh (int localId) const {
		return *(numNeighbors.get()+localId);
	}

	std::tr1::shared_ptr<int> getNumNeigh () const {
		return numNeighbors;
	}

	int getSizeNeighborhoodList () const {
		return sizeNeighborhoodList;
	}



};


Pd_shared_ptr_Array<int> getPointsAxisAlignedMinimum
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon
);


Pd_shared_ptr_Array<int> getPointsAxisAlignedMaximum
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon
);

/*
 * This function finds points that are strictly greater than the
 * input "axisMinimumValue" but are within 1 neighborhood "horizon" of
 * this value
 */

Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon,
		double axisMinimumValue
);

/*
 * This function finds points that are strictly less than the
 * input "axisMaximumValue" but are within 1 neighborhood "horizon" of
 * this value
 */
Pd_shared_ptr_Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CoordinateLabel axis,
		std::tr1::shared_ptr<double> xPtr,
		int numPoints,
		double horizon,
		double axisMaximumValue
);
/*
 * This function produces neighborhoods for each point "i"
 * If point "i" should be included in the list, then add
 * the extra optional argument withSelf=true;  Default
 * behaviour does not include point "i" in the neighborhood
 */
NeighborhoodList getNeighborhoodList
(
		double horizon,
		int numPoints,
		std::tr1::shared_ptr<double> xOwnedPtr,
		vtkSmartPointer<vtkUnstructuredGrid> overlapGrid,
		bool withSelf=false
);

class Coordinates {
public:
	Coordinates(std::tr1::shared_ptr<double> coordinates, int numPoints): numPoints(numPoints), dPtr(coordinates){}
	int getNumPoints() const { return numPoints; }

	class SortComparator;
	friend class SortComparator;
	class SortComparator {
	private:
		std::tr1::shared_ptr<double> dPtr;
		CoordinateLabel label;
		double *d;
	public:
		SortComparator(std::tr1::shared_ptr<double> coordinates, CoordinateLabel label) : dPtr(coordinates), label(label), d(coordinates.get()){}
		bool operator() (int left, int right) { return d[3*left+label]<d[3*right+label]; }
	};
	SortComparator getSortComparator(CoordinateLabel label) const { return SortComparator(dPtr,label); }
	inline std::tr1::shared_ptr<int> getIdentityMap() const;


	class SearchIterator;
	friend class SearchIterator;
	class SearchIterator : public std::iterator<std::forward_iterator_tag,double> {
	public:
		SearchIterator() : numPoints(0), dPtr(), label(PdNeighborhood::X), map(),  P(0), p(0), start(0), end(0) {}
		SearchIterator
		(
				std::tr1::shared_ptr<double> coordinates,
				int numPoints,
				CoordinateLabel label,
				std::tr1::shared_ptr<int> map,
				int P=0
		)
		: numPoints(numPoints),
		  dPtr(coordinates),
		  label(label),
		  map(map),
		  P(P),
		  p(map.get()+P),
		  start(map.get()),
		  end(map.get()+numPoints)
		  {}
		const int* mapIterator() const { return p; }
		const int* mapEnd()      const { return end; }
		const int* mapStart()    const { return start; }
		int numPointsFromStart() const { return P; }
		int numPointsToEnd() const { return numPoints-P; }
	private:
		int numPoints;
		std::tr1::shared_ptr<double> dPtr;
		CoordinateLabel label;
		std::tr1::shared_ptr<int> map;
		int P;
		int *p, *start, *end;
		const double* getD() const { return (dPtr.get()+label); }

	public:

		/*
		 * Need to throw and exception for dereferencing iterator if at "end"
		 */
		double operator*()  const  { const double *d = getD(); return *(d+3*(*p)); }
		double operator->() const  { const double *d = getD(); return *(d+3*(*p)); }
		SearchIterator& operator++()        {if (p!=end) {p++;P++;} return *this;}
		SearchIterator& operator++(int)     { return operator++(); }
		SearchIterator& operator--()        {if (p!=start) {p--;P--;} return *this;}
		SearchIterator& operator--(int)     { return operator--(); }
		const SearchIterator operator+(int n) const {

			if (P+n<=numPoints){
				return SearchIterator(dPtr,numPoints,label,map,P+n);
			}
			else {
				/*
				 * Need to throw an exception here; perhaps just return end
				 */
				return *this;
			}
		}
		SearchIterator operator-(int n) const {

			if (P-n>=0) {
				return SearchIterator(dPtr,numPoints,label,map,P-n);
			}
			else {
				/*
				 * Need to throw an exception here or perhaps just return end
				 */
				return *this;
			}
		}

		SearchIterator& operator+=(int n)    { if (P+n<=numPoints) {p+=n;P+=n;} return *this; }
		SearchIterator& operator-=(int n)    { if (P-n>=0) {p-=n;P-=n;} return *this; }
		bool operator==( const SearchIterator& r) const  { return p==r.p; }
		bool operator!=( const SearchIterator& r) const  { return p!=r.p; }
		bool operator ()  ( const SearchIterator& left,  const SearchIterator& right) const { return *left < *right; }
		bool operator ()  (double left, double right) const { return left < right; }
	};

	SearchIterator begin(CoordinateLabel label,std::tr1::shared_ptr<int> map) const {
		return SearchIterator(dPtr,numPoints,label,map);
	}
	SearchIterator end(CoordinateLabel label,std::tr1::shared_ptr<int> map) const {
		return SearchIterator(dPtr,numPoints,label,map,numPoints);
	}


private:
	int numPoints;
	std::tr1::shared_ptr<double> dPtr;
};

std::tr1::shared_ptr<int> Coordinates::getIdentityMap() const {
	/*
	 * Creates an identity map
	 */
	int N = getNumPoints();
	std::tr1::shared_ptr<int> mapPtr = std::tr1::shared_ptr<int>(new int[N],PdQuickGrid::Deleter<int>());
	int *map = mapPtr.get();
	int *end = map + N;
	for(int j=0;map != end; map++, j++) {*map=j;}
	return mapPtr;
}

}


#endif /* PDNEIGHBORHOOD_H_ */
