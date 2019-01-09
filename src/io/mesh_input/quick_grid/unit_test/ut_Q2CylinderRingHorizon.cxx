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

#include "../QuickGrid.h"
#include "Vector3D.h"
#include <iostream>
#include <cmath>


using namespace QUICKGRID;
using UTILITIES::Vector3D;
using std::shared_ptr;
using std::size_t;


/*
 * For the Q2 cylinder, and mesh discretization used in this file, this
 * function returns the number of neighbors for a point in slab kZ along the axis.
 */
size_t numNeighbors(size_t kZ){

	if(kZ<3)
		return (size_t)(83+kZ*21);
	else if(kZ>=3 && kZ <=55)
		return  (size_t)(83+3*21);
	else
		return  (size_t)(83+(59-kZ)*21);
}


TEUCHOS_UNIT_TEST( Q2CylinderRingHorizon, RingHorizonTest) {


	/*
	 * Construct ring spec
	 */
	// set ring at point = (0,0,0)
	// Note that ring is only translated to the location
	Vector3D center;
	/*
	 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
	 * 	 double innerRadius = 0.020;
	 *   double outerRadius = 0.025;
	 *   int numRings = 3;
	 */
	double innerRadius = 0.020;
	double outerRadius = 0.025;
	double a = .050;
	double cylinderLength = 2.0*a;
	int numRings = 3;
	SpecRing2D ring2dSpec(center,innerRadius,outerRadius,numRings);
	/*
	 * This produces 85 cells in a ring for the Q2 cylinder geometry
	 */
	Spec1D thetaSpec(ring2dSpec.getNumRays(), 0.0, 2.0*M_PI);
	/*
	 * Compute number of cells along length of cylinder so that aspect ratio
	 * is cells is approximately 1
	 *
	 */
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	int numCellsAxis = (int)(cylinderLength/cellSize) + 1;
	Spec1D axisSpec(numCellsAxis,0.0,cylinderLength);

	double SCALE=2.51;
	double horizon = SCALE*cellSize;
	RingHorizon ringHorizon = thetaSpec.getRingCellHorizon(horizon,innerRadius);

	/*
	 * Testing
	 */
	/*
	 * Check cells in ring horizon
	 */
	{
		int cell=80;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(77 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(78 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(79 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(80 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(81 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(82 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(83 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=81;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(78 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(79 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(80 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(81 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(82 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(83 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(84 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=82;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(79 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(80 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(81 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(82 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(83 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(84 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(0 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=83;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(80 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(81 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(82 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(83 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(84 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(0 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(1 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=84;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(81 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(82 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(83 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(84 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(0 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(1 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(2 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=0;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(82 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(83 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(84 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(0 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(1 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(2 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(3 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}
	{
		int cell=1;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(83 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(84 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(0 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(1 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(2 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(3 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(4 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=2;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(84 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(0 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(1 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(2 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(3 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(4 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(5 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=3;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(0 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(1 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(2 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(3 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(4 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(5 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(6 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}

	{
		int cell=4;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		TEST_ASSERT(7 == hIter.numCells());
		TEST_ASSERT(hIter.hasNextCell());
		// loop over cells in horizon
		TEST_ASSERT(1 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(2 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(3 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(4 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(5 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(6 == hIter.nextCell());
		TEST_ASSERT(hIter.hasNextCell());
		TEST_ASSERT(7 == hIter.nextCell());
		TEST_ASSERT(!hIter.hasNextCell());
	}


}


TEUCHOS_UNIT_TEST( Q2CylinderRingHorizon, NeighborhoodSizesTest) {

	/*
	 * Construct ring spec
	 */
	// set ring at point = (0,0,0)
	// Note that ring is only translated to the location
	Vector3D center;
	/*
	 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
	 * 	 double innerRadius = 0.020;
	 *   double outerRadius = 0.025;
	 *   int numRings = 3;
	 */
	double innerRadius = 0.020;
	double outerRadius = 0.025;
	double a = .050;
	double cylinderLength = 2.0*a;
	size_t numRings = 3;
	SpecRing2D ring2dSpec(center,innerRadius,outerRadius,numRings);
	size_t numRays = ring2dSpec.getNumRays();
	/*
	 * This produces 85 cells in a ring for the Q2 cylinder geometry
	 */
	Spec1D thetaSpec(ring2dSpec.getNumRays(), 0.0, 2.0*M_PI);
	/*
	 * Compute number of cells along length of cylinder so that aspect ratio
	 * is cells is approximately 1
	 */
	double SCALE=2.51;
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	int numCellsAxis = (int)(cylinderLength/cellSize)+1;
	Spec1D axisSpec(numCellsAxis,0.0,cylinderLength);
	double horizon = SCALE*cellSize;

	// Create decomposition iterator
	size_t numProcs = 1;
	TensorProductCylinderMeshGenerator cellPerProcIter(numProcs, horizon,ring2dSpec, axisSpec);

	/*
	 * Testing
	 */
	{
		TEST_ASSERT(15300==cellPerProcIter.getNumGlobalCells());
		TEST_ASSERT(85==numRays);
		TEST_ASSERT(60==numCellsAxis);
		size_t i=0,j=0,k=0;
		size_t proc=0;
		Cell3D cellLocator(i,j,k);
		/*
		 * The number below can be hand calculated
		 * k=0 -->83*255
		 * k=1 -->(83+21)*255
		 * k=2 -->(83+21+21)*255
		 * k>=3 -->(83+21+21+21)*255
		 * ...
		 * k=57 -->(83+21+21)*255
		 * k=58 -->(83+21)*255
		 * k=59 -->83*255
		 *
		 * sum all of the above, then add numPoints for storing length of list at each point
		 */
		TEST_ASSERT((255*8508+15300)==cellPerProcIter.getSizeNeighborList(proc,cellLocator));
	}
	/*
	 * On one processor, its possible to calculate the length of the neighborhood list by hand
	 * size = 83
	 */
	/*
	 * Check neighborhoods
	 */
	{
		/*
		 * Every point for k=0 should have 83 neighbors
		 */
		size_t k=0;
		for(size_t j=0;j<numRays;j++){
			for(size_t i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				TEST_ASSERT(83==cellPerProcIter.computeNumNeighbors(i,j,k));
                                TEST_ASSERT(k<60);
				TEST_ASSERT(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k=1 should have 83+21 neighbors
		 */
		size_t k=1;
		for(size_t j=0;j<numRays;j++){
			for(size_t i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				TEST_ASSERT((83+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
                                TEST_ASSERT(k<60);
				TEST_ASSERT(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k=2 should have 83+21+21 neighbors
		 */
		size_t k=2;
		for(size_t j=0;j<numRays;j++){
			for(size_t i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				TEST_ASSERT((83+21+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
                                TEST_ASSERT(k<60);
				TEST_ASSERT(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k=3 should have 83+21+21+21 neighbors
		 */
		size_t k=3;
		for(size_t j=0;j<numRays;j++){
			for(size_t i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				TEST_ASSERT((83+21+21+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
                                TEST_ASSERT(k<60);
				TEST_ASSERT(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k>=3 and k<?? should have 83+21+21+21 neighbors
		 */
		size_t k=4;
		for(size_t j=0;j<numRays;j++){
			for(size_t i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				TEST_ASSERT((83+21+21+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
                                TEST_ASSERT(k<60);
				TEST_ASSERT(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k<=60 up to three planes from end
		 */
		for(size_t k=59;k>55;k--){
			for(size_t j=0;j<numRays;j++){
				for(size_t i=0;i<numRings;i++){
					Cell3D cellLocator(i,j,k);
					TEST_ASSERT((83+(59-k)*21)==cellPerProcIter.computeNumNeighbors(i,j,k));
                                        TEST_ASSERT(k<60);
					TEST_ASSERT(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
				}
			}
		}
	}

}


TEUCHOS_UNIT_TEST( Q2CylinderRingHorizon, NeighborhoodsTest) {


	/*
	 * Construct ring spec
	 */
	// set ring at point = (0,0,0)
	// Note that ring is only translated to the location
	Vector3D center;
	/*
	 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
	 * 	 double innerRadius = 0.020;
	 *   double outerRadius = 0.025;
	 *   int numRings = 3;
	 */
	double innerRadius = 0.020;
	double outerRadius = 0.025;
	double a = .050;
	double cylinderLength = 2.0*a;
	size_t numRings = 3;
	SpecRing2D ring2dSpec(center,innerRadius,outerRadius,numRings);
	/*
	 * This produces 85 cells in a ring for the Q2 cylinder geometry
	 */
	Spec1D thetaSpec(ring2dSpec.getNumRays(), 0.0, 2.0*M_PI);
	/*
	 * Compute number of cells along length of cylinder so that aspect ratio
	 * is cells is approximately 1
	 */
	double SCALE=2.51;
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	size_t numCellsAxis = (int)(cylinderLength/cellSize)+1;
	Spec1D axisSpec(numCellsAxis,0.0,cylinderLength);
	double dz = axisSpec.getCellSize();
	double cellRads = thetaSpec.getCellSize();
	double dr = cellSize;
	double horizon = SCALE*cellSize;

	Array<double> rPtr = getDiscretization(ring2dSpec.getRaySpec());
	Array<double> thetaPtr = getDiscretization(ring2dSpec.getRingSpec());
	Array<double> zPtr = getDiscretization(axisSpec);

	// Create decomposition iterator
	size_t numProcs = 1;
	TensorProductCylinderMeshGenerator cellPerProcIter(numProcs, horizon,ring2dSpec, axisSpec);
	QuickGridData gridData;
	/*
	 * Testing
	 */
	double cylinderVolume = M_PI*(outerRadius*outerRadius-innerRadius*innerRadius)*cylinderLength;
	double sumCellVol=0.0;
	{
		size_t i=0,j=0,k=0,proc=0;
		Cell3D cellLocator(i,j,k);
		QuickGridData pdGridDataProc0 = cellPerProcIter.allocatePdGridData();
		std::pair<Cell3D,QuickGridData> p0Data = cellPerProcIter.computePdGridData(proc,cellLocator,pdGridDataProc0);
		gridData = p0Data.second;
//		Cell3D nextCellLocator = p0Data.first;

		TEST_ASSERT(3==gridData.dimension);
		TEST_ASSERT(15300==gridData.globalNumPoints);
		TEST_ASSERT(15300==gridData.numPoints);
		TEST_ASSERT((255*8508+15300)==gridData.sizeNeighborhoodList);
		TEST_ASSERT(0==gridData.numExport);
		int *gIds = gridData.myGlobalIDs.get();
		for(size_t p=0;p<gridData.numPoints;p++,gIds++)
			TEST_ASSERT((int)p==*gIds);


		/*
		 * Generate a random number of points and assert their coordinates and neighborhoods
		 */
		double *X = gridData.myX.get();
		double *r = rPtr.get();
		double *theta = thetaPtr.get();
		double *z = zPtr.get();
		size_t nx = ring2dSpec.getNumRings();
		size_t ny = thetaSpec.getNumCells();
		int *neighborhoodPtr = gridData.neighborhoodPtr.get();
		int *neighborhood = gridData.neighborhood.get();
		double *vol = gridData.cellVolume.get();

		for(size_t k=0;k<axisSpec.getNumCells();k++){
			size_t ranZid = k; // rand() % axisSpec.getNumCells();
			double ranZ = z[ranZid];
			for(size_t j=0;j<thetaSpec.getNumCells();j++){
				size_t ranThetaid = j; // rand() % thetaSpec.getNumCells();
				double ranTheta = theta[ranThetaid];
				for(size_t i=0;i<ring2dSpec.getNumRings();i++){
					size_t ranRid = i; // rand() % ring2dSpec.getNumRings();
					double ranR = r[ranRid];
					double x = ranR*cos(ranTheta);
					double y = ranR*sin(ranTheta);
					double z = ranZ;
					int gId =  ranRid + ranThetaid * nx + ranZid * nx * ny;
					TEST_ASSERT(x==X[3*gId]);
					TEST_ASSERT(y==X[3*gId+1]);
					TEST_ASSERT(z==X[3*gId+2]);
					int ptr = neighborhoodPtr[gId];
                                        TEST_ASSERT( ranZid < 60);
					TEST_ASSERT((int)numNeighbors(ranZid)==neighborhood[ptr]);

					/*
					 * Volume
					 */
					double v = ranR*dr*cellRads*dz;
					TEST_ASSERT(v==vol[gId]);
					sumCellVol+= vol[gId];
				}
			}
		}
		const double tolerance = 1.0e-10;
		TEST_FLOATING_EQUALITY(sumCellVol,cylinderVolume,tolerance);
	//	std::cout << "cylinderVolume = " << cylinderVolume << "; sumCellVol = " << sumCellVol << std::endl;
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

