/*
 * utPdNeighborhoodProcessorFrame.cxx
 *
 *  Created on: Mar 10, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "PdNeighborhood.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include <vector>
#include <set>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <time.h>


using namespace PdQuickGrid;
using namespace PdNeighborhood;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::vector;
using std::set;

const int NUMBER_OF_POINTS_SORTED = 100000;
const int N = NUMBER_OF_POINTS_SORTED;

const int nx = 6;
const int ny = 6;
const int nz = 1;
const double xStart = 0.0;
const double xLength = 1.0;
const double yStart = 0.0;
const double yLength = 1.0;
const double zStart = 0.0;
const double zLength = 1.0;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;



shared_ptr<double> createPoints(){
	/*
	 * Random coordinates between 0 and PI
	 */
	shared_ptr<double> xPtr(new double[3*N],Deleter<double>());
	{
		/*
		 * Initialize random number generator
		 */
		srand ( time(NULL) );
		double *X = xPtr.get();
		double pi = M_PI;
		for(int p=0;p<N;p++,X+=3){
			*(X+0)= pi*(rand()%N)/N;
			*(X+1)= pi*(rand()%N)/N;
			*(X+2)= pi*(rand()%N)/N;
		}
	}

	return xPtr;
}

void leastUpperBound(){


	// Creates a random set of points
	shared_ptr<double> xPtr = createPoints();
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y,PdNeighborhood::Z};

	for(int label=0;label<3;label++){
		/*
		 * create a random value between 0 and pi
		 */
		double value = M_PI*(rand()%N)/N;

		Coordinates c(xPtr,N);
		/*
		 * Sort each coordinate direction
		 */
		Coordinates::SortComparator compare = c.getSortComparator(labels[label]);
		std::tr1::shared_ptr<int> mapX = c.getIdentityMap();
		std::sort(mapX.get(),mapX.get()+N,compare);
		const Coordinates::SearchIterator start=c.begin(labels[label],mapX);
		const Coordinates::SearchIterator end=start+N;
		/*
		 * find least upper bound from coordinates given "value"
		 */
		Coordinates::SearchIterator f = std::upper_bound(start,end,value);

		/*
		 * assert that all points including first are greater than "value"
		 */
		for(;f!=end; f++){
			BOOST_CHECK(value<=*f);
		}

	}

}

void noCellHorizon() {

	int numProcs=1;
	int myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	/*
	 * Assert that there are no neighbors
	 */
	BOOST_CHECK(numCells == decomp.globalNumPoints);
	int size = decomp.sizeNeighborhoodList;
	BOOST_CHECK(size=numCells);
	int *list = decomp.neighborhood.get();
	int *end = list+size;
	for(;list!=end;){
		int numNeigh = *list; list++;
		BOOST_CHECK(0==numNeigh);
		int *e = list+numNeigh;
		for(;list!=e;list++ ){
			BOOST_CHECK(false);
		}
	}

}

void noCellHorizonFrameLeftAndRight() {

	int numProcs=1;
	int myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Coordinates c(xPtr,numCells);


	/*
	 * Find edge points on left and right (xMin,xMax)
	 */

	/*
	 * Point ids on edges
	 */
	int left[]   = {0,6,12,18,24,30};
	int right[]   = {5,11,17,23,29,35};

	int numAxes=2;
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y};
	std::vector<std::tr1::shared_ptr<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Coordinates::SortComparator compare = c.getSortComparator(labels[j]);
		sortedMaps[j] = c.getIdentityMap();
		/*
		 * Sort points
		 */
		std::sort(sortedMaps[j].get(),sortedMaps[j].get()+numCells,compare);

	}

	{
		/*
		 * Now find points along LEFT side;
		 * Find least upper bound of points for xMin+horizon
		 */
		int axis = PdNeighborhood::X;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * First value in map corresponds with minimum value of x
		 */
		int iXMIN = *mapX;
		double xMin = x[3*iXMIN];
		const double tolerance = 1.0e-15;
		BOOST_CHECK_CLOSE(xMin,(xLength/nx/2.0),tolerance);
		double value = xMin + horizon;
		const Coordinates::SearchIterator start=c.begin(labels[axis],sortedMaps[axis]);
		const Coordinates::SearchIterator end=start+numCells;
		Coordinates::SearchIterator upper = std::upper_bound(start,end,value);
		int numPointsInSet = upper.numPointsFromStart();
		BOOST_CHECK(numPointsInSet==6);
		std::set<int> found(upper.mapStart(),upper.mapIterator());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=left;b!=left+numPointsInSet;b++){
			BOOST_CHECK(setEnd!=found.find(*b));
		}
	}

	{
		/*
		 * Now find points along RIGHT side;
		 * Find least upper bound of points for xMin+horizon
		 */
		int axis = PdNeighborhood::X;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * Last point corresponds with minimum value
		 */
		int iXMAX = *(mapX+numCells-1);
		double xMax = x[3*iXMAX];
		const double tolerance = 1.0e-13;
		BOOST_CHECK_CLOSE(xMax,(xLength-xLength/nx/2.0),tolerance);
		double value = xMax - horizon;
		const Coordinates::SearchIterator start=c.begin(labels[axis],sortedMaps[axis]);
		const Coordinates::SearchIterator end=start+numCells;
		Coordinates::SearchIterator lower = std::lower_bound(start,end,value);
		int numPointsInSet = lower.numPointsToEnd();
		BOOST_CHECK(numPointsInSet==6);
		std::set<int> found(lower.mapIterator(),lower.mapEnd());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=right;b!=right+numPointsInSet;b++){
			BOOST_CHECK(setEnd!=found.find(*b));
		}
	}

}

void noCellHorizonFrameTopAndBottom() {

	int numProcs=1;
	int myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Coordinates c(xPtr,numCells);


	/*
	 * Find points on top and bottom edges
	 */

	/*
	 * Point ids on edges
	 */
	int bottom[]  = {0,1,3,4,4,5};
	int top[]     = {30,31,32,33,34,35};

	int numAxes=2;
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y};
	std::vector<std::tr1::shared_ptr<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Coordinates::SortComparator compare = c.getSortComparator(labels[j]);
		sortedMaps[j] = c.getIdentityMap();
		/*
		 * Sort points
		 */
		std::sort(sortedMaps[j].get(),sortedMaps[j].get()+numCells,compare);

	}

	{
		/*
		 * Now find points along BOTTOM side;
		 * Find least upper bound of points for xMin+horizon
		 */
		int axis = PdNeighborhood::Y;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * First value in map corresponds with minimum value of x
		 */
		int iXMIN = *mapX;
		double xMin = x[3*iXMIN+axis];
		const double tolerance = 1.0e-15;
		BOOST_CHECK_CLOSE(xMin,(xLength/nx/2.0),tolerance);
		double value = xMin + horizon;
		const Coordinates::SearchIterator start=c.begin(labels[axis],sortedMaps[axis]);
		const Coordinates::SearchIterator end=start+numCells;
		Coordinates::SearchIterator upper = std::upper_bound(start,end,value);
		int numPointsInSet = upper.numPointsFromStart();
		BOOST_CHECK(numPointsInSet==6);
		std::set<int> found(upper.mapStart(),upper.mapIterator());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=bottom;b!=bottom+numPointsInSet;b++){
			BOOST_CHECK(setEnd!=found.find(*b));
		}
	}

	{
		/*
		 * Now find points along TOP side;
		 * Find least upper bound of points for xMin+horizon
		 */
		int axis = PdNeighborhood::Y;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * Last point corresponds with maximum value
		 */
		int iXMAX = *(mapX+numCells-1);
		double xMax = x[3*iXMAX+axis];
		const double tolerance = 1.0e-13;
		BOOST_CHECK_CLOSE(xMax,(xLength-xLength/nx/2.0),tolerance);
		double value = xMax - horizon;
		const Coordinates::SearchIterator start=c.begin(labels[axis],sortedMaps[axis]);
		const Coordinates::SearchIterator end=start+numCells;
		Coordinates::SearchIterator lower = std::lower_bound(start,end,value);
		int numPointsInSet = lower.numPointsToEnd();
		BOOST_CHECK(numPointsInSet==6);
		std::set<int> found(lower.mapIterator(),lower.mapEnd());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=top;b!=top+numPointsInSet;b++){
			BOOST_CHECK(setEnd!=found.find(*b));
		}
	}

}

void greatestLowerBound(){


	// Creates a random set of points
	shared_ptr<double> xPtr = createPoints();
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y,PdNeighborhood::Z};

	for(int label=0;label<3;label++){
		/*
		 * create a random value between 0 and pi
		 */
		double value = M_PI*(rand()%N)/N;

		Coordinates c(xPtr,N);
		/*
		 * Sort each coordinate direction
		 */
		Coordinates::SortComparator compare = c.getSortComparator(labels[label]);
		std::tr1::shared_ptr<int> mapX = c.getIdentityMap();
		std::sort(mapX.get(),mapX.get()+N,compare);
		const Coordinates::SearchIterator start=c.begin(labels[label],mapX);
		const Coordinates::SearchIterator end=start+N;
		/*
		 * find greatest lower bound from coordinates given "value"
		 * NOTE that "f" is NOT in the set -- this is different from the
		 * upper_bound case which actually includes the point in the set
		 */
		Coordinates::SearchIterator f = std::lower_bound(start,end,value);

		/*
		 * assert that all points less than "value"
		 */
		Coordinates::SearchIterator iter=c.begin(labels[label],mapX);
		for(;iter!=f; iter++){
			BOOST_CHECK(value>=*iter);
		}

	}

}
void constructFrameNoHorizon(){
	int numProcs=1;
	int myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Coordinates c(xPtr,numCells);

	int numAxes=2;
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y};
	std::vector<std::tr1::shared_ptr<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Coordinates::SortComparator compare = c.getSortComparator(labels[j]);
		sortedMaps[j] = c.getIdentityMap();
		/*
		 * Sort points
		 */
		std::sort(sortedMaps[j].get(),sortedMaps[j].get()+numCells,compare);

	}

	unsigned int numFrameCells=20;
	vector<int> frameCells(numFrameCells);
	int frame[] = {
			0,1,2,3,4,5,
			30,31,32,33,34,35,
			6,12,18,24,
			11,17,23,29
	};
	frameCells.assign(frame,frame+numFrameCells);
	BOOST_CHECK(frameCells.capacity()==numFrameCells);

	/*
	 * Loop over axes and collect points at min and max ranges
	 * Add Points to frame set
	 */
	std::set<int> frameSet;
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
			const Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const Coordinates::SearchIterator end=start+numCells;
			Coordinates::SearchIterator lub = std::upper_bound(start,end,value);
			int numPointsInSet = lub.numPointsFromStart();
			BOOST_CHECK(numPointsInSet==6);
			frameSet.insert(lub.mapStart(),lub.mapIterator());

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
			const Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const Coordinates::SearchIterator end=start+numCells;
			Coordinates::SearchIterator glb = std::upper_bound(start,end,value);
			int numPointsInSet = glb.numPointsToEnd();
			BOOST_CHECK(numPointsInSet==6);
			frameSet.insert(glb.mapIterator(),glb.mapEnd());
		}
	}

	/*
	 * Now compare frame cells with computed "frameSet"
	 */
	vector<int>::iterator i = frameCells.begin();
	const vector<int>::iterator end = frameCells.end();
	const set<int>::iterator setEnd = frameSet.end();
	for(;i!=end;i++){
		BOOST_CHECK(setEnd!=frameSet.find(*i));
	}


}

void constructFrameOneCellHorizon(){
	int numProcs=1;
	int myRank=0;
	double SCALE=1.1;
	double horizon = SCALE*xSpec.getCellSize();
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Coordinates c(xPtr,numCells);

	int numAxes=2;
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y};
	std::vector<std::tr1::shared_ptr<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Coordinates::SortComparator compare = c.getSortComparator(labels[j]);
		sortedMaps[j] = c.getIdentityMap();
		/*
		 * Sort points
		 */
		std::sort(sortedMaps[j].get(),sortedMaps[j].get()+numCells,compare);

	}

	unsigned int numFrameCells=32;
	vector<int> frameCells(numFrameCells);
	int frame[] = {
			0,1,2,3,4,5,
			30,31,32,33,34,35,
			6,12,18,24,
			11,17,23,29,
			7,8,9,10,
			25,26,27,28,
			13,16,
			19,22
	};
	frameCells.assign(frame,frame+numFrameCells);
	BOOST_CHECK(frameCells.capacity()==numFrameCells);

	/*
	 * Loop over axes and collect points at min and max ranges
	 * Add Points to frame set
	 */
	std::set<int> frameSet;
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
			const Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const Coordinates::SearchIterator end=start+numCells;
			Coordinates::SearchIterator lub = std::upper_bound(start,end,value);
			int numPointsInSet = lub.numPointsFromStart();
			BOOST_CHECK(numPointsInSet==12);
			frameSet.insert(lub.mapStart(),lub.mapIterator());

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
			const Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const Coordinates::SearchIterator end=start+numCells;
			Coordinates::SearchIterator glb = std::upper_bound(start,end,value);
			int numPointsInSet = glb.numPointsToEnd();
			BOOST_CHECK(numPointsInSet==12);
			frameSet.insert(glb.mapIterator(),glb.mapEnd());
		}
	}

	/*
	 * Now compare frame cells with computed "frameSet"
	 */
	vector<int>::iterator i = frameCells.begin();
	const vector<int>::iterator end = frameCells.end();
	const set<int>::iterator setEnd = frameSet.end();
	for(;i!=end;i++){
		BOOST_CHECK(setEnd!=frameSet.find(*i));
	}


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPdNeighborhoodProcessorFrame" );
	proc->add(BOOST_TEST_CASE( &createPoints ));
	proc->add(BOOST_TEST_CASE( &leastUpperBound ));
	proc->add(BOOST_TEST_CASE( &greatestLowerBound ));
	proc->add(BOOST_TEST_CASE( &noCellHorizon ));
	proc->add(BOOST_TEST_CASE( &noCellHorizonFrameLeftAndRight ));
	proc->add(BOOST_TEST_CASE( &noCellHorizonFrameTopAndBottom ));
	proc->add(BOOST_TEST_CASE( &constructFrameNoHorizon ));
	proc->add(BOOST_TEST_CASE( &constructFrameOneCellHorizon ));
	framework::master_test_suite().add( proc );
	return success;

}

bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(
		int argc,
		char* argv[]
)
{


	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
