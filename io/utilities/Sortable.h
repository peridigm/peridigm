/*
 * SortableCartesianComponents.h
 *
 *  Created on: Mar 11, 2011
 *      Author: jamitch
 */

#ifndef SORTABLE_H_
#define SORTABLE_H_

#include <set>
#include "Array.h"

namespace UTILITIES {

enum CartesianComponent { X=0, Y, Z };

class Sortable : public Array<double> {

public:
	/*
	 * Expectation is that incoming coordinates are an array of tuples; each tuple has length 3; x1,y1,z1, x2,y2,z2, x3,y3,z3, ...
	 */
	Sortable(size_t num_points, shared_ptr<double> coordinates): Array<double>(3*num_points,coordinates), num_points(num_points) {}

	size_t get_num_points() const { return num_points; }

	class Comparator;
	friend class Comparator;
	class Comparator {
	private:
		shared_ptr<double> dPtr;
		CartesianComponent component;
		const double *d;
	public:
		Comparator(shared_ptr<double> coordinates, CartesianComponent xyz) : dPtr(coordinates), component(xyz), d(coordinates.get()){}
		bool operator() (int left, int right) const { return d[3*left+component] < d[3*right+component]; }
	};

	Comparator getComparator(CartesianComponent xyz) const { return Comparator(get_shared_ptr(),xyz); }

	inline Array<int> getIdentityMap() const;


	class SearchIterator;
	friend class SearchIterator;
	class SearchIterator : public std::iterator<std::forward_iterator_tag,double> {
	public:
		SearchIterator() : num_points(0), dPtr(), component(UTILITIES::X), map(),  P(0), p(0), start(0), end(0) {}
		SearchIterator
		(
				shared_ptr<double> coordinates,
				int num_points,
				CartesianComponent component,
				shared_ptr<int> map,
				int P=0
		)
		: num_points(num_points),
		  dPtr(coordinates),
		  component(component),
		  map(map),
		  P(P),
		  p(map.get()+P),
		  start(map.get()),
		  end(map.get()+num_points)
		{}
		const int* mapIterator() const { return p; }
		const int* mapEnd()      const { return end; }
		const int* mapStart()    const { return start; }
		int numPointsFromStart() const { return P; }
		int numPointsToEnd() const { return num_points-P; }
	private:
		int num_points;
		shared_ptr<double> dPtr;
		CartesianComponent component;
		shared_ptr<int> map;
		int P;
		int *p;
		const int *start, *end;
		const double* getD() const { return (dPtr.get()+component); }

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

			if (P+n<=num_points){
				return SearchIterator(dPtr,num_points,component,map,P+n);
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
				return SearchIterator(dPtr,num_points,component,map,P-n);
			}
			else {
				/*
				 * Need to throw an exception here or perhaps just return end
				 */
				return *this;
			}
		}

		SearchIterator& operator+=(int n)    { if (P+n<=num_points) {p+=n;P+=n;} return *this; }
		SearchIterator& operator-=(int n)    { if (P-n>=0) {p-=n;P-=n;} return *this; }
		bool operator==( const SearchIterator& r) const  { return p==r.p; }
		bool operator!=( const SearchIterator& r) const  { return p!=r.p; }
		bool operator ()  ( const SearchIterator& left,  const SearchIterator& right) const { return *left < *right; }
		bool operator ()  (double left, double right) const { return left < right; }
	};

	SearchIterator begin(CartesianComponent component, shared_ptr<int> map) const {
		return SearchIterator(get_shared_ptr(),num_points,component,map);
	}
	SearchIterator end(CartesianComponent component, shared_ptr<int> map) const {
		return SearchIterator(get_shared_ptr(),num_points,component,map,num_points);
	}

	private:
		size_t num_points;
};

Array<int> Sortable::getIdentityMap() const {
	/*
	 * Creates an identity map
	 */
	size_t N = get_num_points();
	Array<int> mapPtr(N);
	int *map = mapPtr.get();
	int *end = map + N;
	for(int j=0;map != end; map++, j++) {*map=j;}
	return mapPtr;
}



Array<int> getSortedMap(const Sortable& c, CartesianComponent axis);

Array<int> getPointsInNeighborhoodOfAxisAlignedMinimumValue
(
		CartesianComponent axis,
		shared_ptr<double> xPtr,
		size_t numPoints,
		double horizon,
		double axisMinimumValue
);

Array<int> getPointsInNeighborhoodOfAxisAlignedMaximumValue
(
		CartesianComponent axis,
		shared_ptr<double> xPtr,
		size_t numPoints,
		double horizon,
		double axisMaximumValue
);

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
);

Array<int> getPointsAxisAlignedMaximum
(
		CartesianComponent axis,
		shared_ptr<double> xPtr,
		size_t numPoints,
		double horizon
);

shared_ptr< std::set<int> > constructFrameSet(size_t num_owned_points, shared_ptr<double> owned_x, double horizon);


}

#endif /* SORTABLE_H_ */
