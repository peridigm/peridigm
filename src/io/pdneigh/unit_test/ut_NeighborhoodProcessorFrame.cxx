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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_as.hpp"
#include "quick_grid/QuickGrid.h"
#include "Sortable.h"
#include "Array.h"
#include "../NeighborhoodList.h"
#include <vector>
#include <set>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <time.h>


using UTILITIES::CartesianComponent;
using UTILITIES::Sortable;
using namespace PDNEIGH;
using std::shared_ptr;

using std::vector;
using std::set;

const size_t NUMBER_OF_POINTS_SORTED = 100000;
const size_t N = NUMBER_OF_POINTS_SORTED;

const size_t nx = 6;
const size_t ny = 6;
const size_t nz = 1;
const double xStart = 0.0;
const double xLength = 1.0;
const double yStart = 0.0;
const double yLength = 1.0;
const double zStart = 0.0;
const double zLength = 1.0;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;

Array<double> createPoints(){
	/*
	 * Random coordinates between 0 and PI
	 */
	Array<double> xPtr(3*N);
	{
		/*
		 * Initialize random number generator
		 */
		srand ( time(NULL) );
		double *X = xPtr.get();
		double pi = M_PI;
		for(size_t p=0;p<N;p++,X+=3){
			*(X+0)= pi*(rand()%N)/N;
			*(X+1)= pi*(rand()%N)/N;
			*(X+2)= pi*(rand()%N)/N;
		}
	}

	return xPtr;
}


TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame,  CreatePointsTest) {

	/*
	 * Random coordinates between 0 and PI
	 */
	Array<double> xPtr(3*N);
	xPtr = createPoints();
}

TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame, LeastUpperBoundTest) {


	// Creates a random set of points
	Array<double> xPtr = createPoints();
	CartesianComponent labels[] = {UTILITIES::X, UTILITIES::Y, UTILITIES::Z};

	for(int label=0;label<3;label++){
		/*
		 * create a random value between 0 and pi
		 */
		double value = M_PI*(rand()%N)/N;

		Sortable c(N, xPtr.get_shared_ptr());
		/*
		 * Sort each coordinate direction
		 */
		Sortable::Comparator compare = c.getComparator(labels[label]);
		Array<int> mapX = c.getIdentityMap();
		std::sort(mapX.get(),mapX.get()+N,compare);
		const Sortable::SearchIterator start=c.begin(labels[label],mapX.get_shared_ptr());
		const Sortable::SearchIterator end=start+N;
		/*
		 * find least upper bound from coordinates given "value"
		 */
		Sortable::SearchIterator f = std::upper_bound(start,end,value);

		/*
		 * assert that all points including first are greater than "value"
		 */
		for(;f!=end; f++){
			TEST_ASSERT(value<=*f);
		}

	}

}

TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame, NoCellHorizonTest) {


	size_t numProcs=1;
	size_t myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	/*
	 * Assert that there are no neighbors
	 */
	TEST_ASSERT(numCells == decomp.globalNumPoints);
	int size = decomp.sizeNeighborhoodList;
	TEST_ASSERT((unsigned int)size == numCells);
	int *list = decomp.neighborhood.get();
	int *end = list+size;
	for(;list!=end;){
		int numNeigh = *list; list++;
		TEST_ASSERT(0==numNeigh);
		int *e = list+numNeigh;
		for(;list!=e;list++ ){
			TEST_ASSERT(false);
		}
	}

}

TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame, NoCellHorizonFrameLeftAndRightTest) {


	size_t numProcs=1;
	size_t myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Sortable c(numCells,xPtr);


	/*
	 * Find edge points on left and right (xMin,xMax)
	 */

	/*
	 * Point ids on edges
	 */
	int left[]   = {0,6,12,18,24,30};
	int right[]   = {5,11,17,23,29,35};

	int numAxes=2;
	CartesianComponent labels[] = {UTILITIES::X, UTILITIES::Y};
	std::vector<Array<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Sortable::Comparator compare = c.getComparator(labels[j]);
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
		int axis = UTILITIES::X;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * First value in map corresponds with minimum value of x
		 */
		int iXMIN = *mapX;
		double xMin = x[3*iXMIN];
		const double tolerance = 1.0e-15;
		TEST_FLOATING_EQUALITY(xMin,(xLength/nx/2.0),tolerance);
		double value = xMin + horizon;
		const Sortable::SearchIterator start=c.begin(labels[axis],sortedMaps[axis].get_shared_ptr());
		const Sortable::SearchIterator end=start+numCells;
		Sortable::SearchIterator upper = std::upper_bound(start,end,value);
		int numPointsInSet = upper.numPointsFromStart();
		TEST_ASSERT(numPointsInSet==6);
		std::set<int> found(upper.mapStart(),upper.mapIterator());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=left;b!=left+numPointsInSet;b++){
			TEST_ASSERT(setEnd!=found.find(*b));
		}
	}

	{
		/*
		 * Now find points along RIGHT side;
		 * Find least upper bound of points for xMin+horizon
		 */
		int axis = UTILITIES::X;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * Last point corresponds with minimum value
		 */
		int iXMAX = *(mapX+numCells-1);
		double xMax = x[3*iXMAX];
		const double tolerance = 1.0e-13;
		TEST_FLOATING_EQUALITY(xMax,(xLength-xLength/nx/2.0),tolerance);
		double value = xMax - horizon;
		const Sortable::SearchIterator start=c.begin(labels[axis],sortedMaps[axis].get_shared_ptr());
		const Sortable::SearchIterator end=start+numCells;
		Sortable::SearchIterator lower = std::lower_bound(start,end,value);
		int numPointsInSet = lower.numPointsToEnd();
		TEST_ASSERT(numPointsInSet==6);
		std::set<int> found(lower.mapIterator(),lower.mapEnd());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=right;b!=right+numPointsInSet;b++){
			TEST_ASSERT(setEnd!=found.find(*b));
		}
	}

}

TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame, NoCellHorizonFrameTopAndBottomTest) {

	size_t numProcs=1;
	size_t myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Sortable c(numCells,xPtr);


	/*
	 * Find points on top and bottom edges
	 */

	/*
	 * Point ids on edges
	 */
	int bottom[]  = {0,1,3,4,4,5};
	int top[]     = {30,31,32,33,34,35};

	int numAxes=2;
	CartesianComponent labels[] = {UTILITIES::X, UTILITIES::Y};
	std::vector<Array<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Sortable::Comparator compare = c.getComparator(labels[j]);
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
		int axis = UTILITIES::Y;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * First value in map corresponds with minimum value of x
		 */
		int iXMIN = *mapX;
		double xMin = x[3*iXMIN+axis];
		const double tolerance = 1.0e-15;
		TEST_FLOATING_EQUALITY(xMin,(xLength/nx/2.0),tolerance);
		double value = xMin + horizon;
		const Sortable::SearchIterator start=c.begin(labels[axis],sortedMaps[axis].get_shared_ptr());
		const Sortable::SearchIterator end=start+numCells;
		Sortable::SearchIterator upper = std::upper_bound(start,end,value);
		int numPointsInSet = upper.numPointsFromStart();
		TEST_ASSERT(numPointsInSet==6);
		std::set<int> found(upper.mapStart(),upper.mapIterator());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=bottom;b!=bottom+numPointsInSet;b++){
			TEST_ASSERT(setEnd!=found.find(*b));
		}
	}

	{
		/*
		 * Now find points along TOP side;
		 * Find least upper bound of points for xMin+horizon
		 */
		int axis = UTILITIES::Y;
		double *x = xPtr.get();
		// x map
		int *mapX = sortedMaps[axis].get();
		/*
		 * Last point corresponds with maximum value
		 */
		int iXMAX = *(mapX+numCells-1);
		double xMax = x[3*iXMAX+axis];
		const double tolerance = 1.0e-13;
		TEST_FLOATING_EQUALITY(xMax,(xLength-xLength/nx/2.0),tolerance);
		double value = xMax - horizon;
		const Sortable::SearchIterator start=c.begin(labels[axis],sortedMaps[axis].get_shared_ptr());
		const Sortable::SearchIterator end=start+numCells;
		Sortable::SearchIterator lower = std::lower_bound(start,end,value);
		int numPointsInSet = lower.numPointsToEnd();
		TEST_ASSERT(numPointsInSet==6);
		std::set<int> found(lower.mapIterator(),lower.mapEnd());
		std::set<int>::iterator setEnd = found.end();
		for(int *b=top;b!=top+numPointsInSet;b++){
			TEST_ASSERT(setEnd!=found.find(*b));
		}
	}

}

TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame, GreatestLowerBoundTest) {


	// Creates a random set of points
	Array<double> xPtr = createPoints();
	CartesianComponent labels[] = {UTILITIES::X, UTILITIES::Y, UTILITIES::Z};

	for(int label=0;label<3;label++){
		/*
		 * create a random value between 0 and pi
		 */
		double value = M_PI*(rand()%N)/N;

		Sortable c(N, xPtr.get_shared_ptr());
		/*
		 * Sort each coordinate direction
		 */
		Sortable::Comparator compare = c.getComparator(labels[label]);
		Array<int> mapX = c.getIdentityMap();
		std::sort(mapX.get(),mapX.get()+N,compare);
		const Sortable::SearchIterator start=c.begin(labels[label],mapX.get_shared_ptr());
		const Sortable::SearchIterator end=start+N;
		/*
		 * find greatest lower bound from coordinates given "value"
		 * NOTE that "f" is NOT in the set -- this is different from the
		 * upper_bound case which actually includes the point in the set
		 */
		Sortable::SearchIterator f = std::lower_bound(start,end,value);

		/*
		 * assert that all points less than "value"
		 */
		Sortable::SearchIterator iter=c.begin(labels[label],mapX.get_shared_ptr());
		for(;iter!=f; iter++){
			TEST_ASSERT(value>=*iter);
		}

	}

}

TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame, ConstructFrameNoHorizonTest) {

	size_t numProcs=1;
	size_t myRank=0;
	double SCALE=0.4;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Sortable c(numCells,xPtr);

	int numAxes=2;
	CartesianComponent labels[] = {UTILITIES::X, UTILITIES::Y};
	std::vector<Array<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Sortable::Comparator compare = c.getComparator(labels[j]);
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
	TEST_ASSERT(frameCells.capacity()==numFrameCells);

	/*
	 * Loop over axes and collect points at min and max ranges
	 * Add Points to frame set
	 */
	std::set<int> frameSet;
	UTILITIES::CartesianComponent *label = labels;
	UTILITIES::CartesianComponent *endLabel = labels+numAxes;
	for(;label!=endLabel;label++){

		{
			/*
			 * MINIMUM
			 * Find least upper bound of points for Min+horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			Array<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();
			/*
			 * First value in map corresponds with minimum value
			 */
			int iMIN = *map;
			double min = x[3*iMIN+axis];
			double value = min + horizon;
			const Sortable::SearchIterator start=c.begin(*label,mapPtr.get_shared_ptr());
			const Sortable::SearchIterator end=start+numCells;
			Sortable::SearchIterator lub = std::upper_bound(start,end,value);
			int numPointsInSet = lub.numPointsFromStart();
			TEST_ASSERT(numPointsInSet==6);
			frameSet.insert(lub.mapStart(),lub.mapIterator());

		}

		{
			/*
			 * MAXIMUM
			 * Find greatest lower bound glb for Max-horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			Array<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();

			/*
			 * Last value in map corresponds with maximum value
			 */
			int iMAX = *(map+numCells-1);
			double max = x[3*iMAX+axis];
			double value = max - horizon;
			const Sortable::SearchIterator start=c.begin(*label,mapPtr.get_shared_ptr());
			const Sortable::SearchIterator end=start+numCells;
			Sortable::SearchIterator glb = std::upper_bound(start,end,value);
			int numPointsInSet = glb.numPointsToEnd();
			TEST_ASSERT(numPointsInSet==6);
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
		TEST_ASSERT(setEnd!=frameSet.find(*i));
	}
}


TEUCHOS_UNIT_TEST( NeighborhoodProcessorFrame, ConstructFrameOneCellHorizonTest) {

	size_t numProcs=1;
	size_t myRank=0;
	double SCALE=1.1;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	const Sortable c(numCells,xPtr);

	int numAxes=2;
	CartesianComponent labels[] = {UTILITIES::X, UTILITIES::Y};
	std::vector<Array<int> > sortedMaps(2);
	for(int j=0;j<numAxes;j++){

		Sortable::Comparator compare = c.getComparator(labels[j]);
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
	TEST_ASSERT(frameCells.capacity()==numFrameCells);

	/*
	 * Loop over axes and collect points at min and max ranges
	 * Add Points to frame set
	 */
	std::set<int> frameSet;
	UTILITIES::CartesianComponent *label = labels;
	UTILITIES::CartesianComponent *endLabel = labels+numAxes;
	for(;label!=endLabel;label++){

		{
			/*
			 * MINIMUM
			 * Find least upper bound of points for Min+horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			Array<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();
			/*
			 * First value in map corresponds with minimum value
			 */
			int iMIN = *map;
			double min = x[3*iMIN+axis];
			double value = min + horizon;
			const Sortable::SearchIterator start=c.begin(*label,mapPtr.get_shared_ptr());
			const Sortable::SearchIterator end=start+numCells;
			Sortable::SearchIterator lub = std::upper_bound(start,end,value);
			int numPointsInSet = lub.numPointsFromStart();
			TEST_ASSERT(numPointsInSet==12);
			frameSet.insert(lub.mapStart(),lub.mapIterator());

		}

		{
			/*
			 * MAXIMUM
			 * Find greatest lower bound glb for Max-horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			Array<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();

			/*
			 * Last value in map corresponds with maximum value
			 */
			int iMAX = *(map+numCells-1);
			double max = x[3*iMAX+axis];
			double value = max - horizon;
			const Sortable::SearchIterator start=c.begin(*label,mapPtr.get_shared_ptr());
			const Sortable::SearchIterator end=start+numCells;
			Sortable::SearchIterator glb = std::upper_bound(start,end,value);
			int numPointsInSet = glb.numPointsToEnd();
			TEST_ASSERT(numPointsInSet==12);
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
		TEST_ASSERT(setEnd!=frameSet.find(*i));
	}


}



int main
(
		int argc,
		char* argv[]
)
{


	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
