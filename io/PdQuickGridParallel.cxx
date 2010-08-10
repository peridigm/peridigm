/*
 * PdQuickGrid3D.cxx
 *
 *  Created on: Nov 9, 2009
 *      Author: jamitch
 */

#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdNeighborhood.h"
#include "mpi.h"
#include<iostream>

namespace PdQuickGrid {

using namespace std;
using std::tr1::shared_ptr;
using namespace PdNeighborhood;

PointBB::PointBB(const double *x, double h)
: xMin(*(x+0)-h), xMax(*(x+0)+h),
  yMin(*(x+1)-h), yMax(*(x+1)+h),
  zMin(*(x+2)-h), zMax(*(x+2)+h)
{}

/**
 * This function produces an unbalanced discretization (although for some geometries it
 * may not be too bad
 */
PdGridData getDiscretization(int rank, PdQuickGridMeshGenerationIterator &cellIter)
{
	MPI_Status status;
	int ack = 0;
	int ackTag = 0;
	int numPointsTag=1;
	int idsTag = 2;
	int coordinatesTag = 3;
	int globalNumPointsTag=4;
	int sizeNeighborhoodListTag=5;
	int neighborhoodTag=6;
	int volumeTag=7;
	int neighborhoodPtrTag=8;
	PdGridData gridData;
	int dimension = cellIter.getDimension();

	if(0 == rank){

		PdGridData pdGridDataProc0 = cellIter.allocatePdGridData();
		PdGridData  pdGridDataProcN = cellIter.allocatePdGridData();
		std::pair<Cell3D,PdGridData> p0Data = cellIter.beginIterateProcs(pdGridDataProc0);
		gridData = p0Data.second;
		Cell3D nextCellLocator = p0Data.first;


		while(cellIter.hasNextProc()){
			int proc = cellIter.proc();
			std::pair<Cell3D,PdGridData> data = cellIter.nextProc(nextCellLocator,pdGridDataProcN);
			PdGridData gridData = data.second;
			nextCellLocator = data.first;


			// Need to send this data to proc
			int globalNumPoints = gridData.globalNumPoints;
			int numPoints = gridData.numPoints;
			int sizeNeighborhoodList = gridData.sizeNeighborhoodList;
			shared_ptr<int> gIds = gridData.myGlobalIDs;
			shared_ptr<double> X = gridData.myX;
			shared_ptr<double> V = gridData.cellVolume;
			shared_ptr<int> neighborhood = gridData.neighborhood;
			shared_ptr<int> neighborhoodPtr = gridData.neighborhoodPtr;

			MPI_Send(&numPoints, 1, MPI_INT, proc, numPointsTag, MPI_COMM_WORLD);
			MPI_Recv(&ack, 1, MPI_INT, proc, ackTag, MPI_COMM_WORLD, &status);
			MPI_Send(&globalNumPoints, 1, MPI_INT, proc, globalNumPointsTag, MPI_COMM_WORLD);
			MPI_Send(&sizeNeighborhoodList, 1, MPI_INT, proc, sizeNeighborhoodListTag, MPI_COMM_WORLD);
			MPI_Send(gIds.get(), numPoints, MPI_INT, proc, idsTag, MPI_COMM_WORLD);
			MPI_Send(X.get(), dimension*numPoints, MPI_DOUBLE, proc, coordinatesTag, MPI_COMM_WORLD);
			MPI_Send(V.get(),numPoints, MPI_DOUBLE, proc, volumeTag, MPI_COMM_WORLD);
			MPI_Send(neighborhood.get(), sizeNeighborhoodList, MPI_INT, proc, neighborhoodTag, MPI_COMM_WORLD);
			MPI_Send(neighborhoodPtr.get(), numPoints, MPI_INT, proc, neighborhoodPtrTag, MPI_COMM_WORLD);

		}
	    /* signal all procs it is OK to go on */
	    ack = 0;
	    for(int proc=1;proc<cellIter.getNumProcs();proc++){
	      MPI_Send(&ack, 1, MPI_INT, proc, ackTag, MPI_COMM_WORLD);
	    }

	}
	else {
		// Receive data from processor 0
		// Create this procs 'GridData'
		int numPoints=0;
		int globalNumPoints = 0;
		int sizeNeighborhoodList = 0;
		MPI_Recv(&numPoints, 1, MPI_INT, 0, numPointsTag, MPI_COMM_WORLD, &status);
		ack = 0;
		if (numPoints > 0){
			PdGridData 	gData = PdQuickGrid::allocatePdGridData(numPoints,dimension);
			std::tr1::shared_ptr<double> g=gData.myX;
			std::tr1::shared_ptr<double> cellVolume=gData.cellVolume;
			std::tr1::shared_ptr<int> gIds=gData.myGlobalIDs;
			std::tr1::shared_ptr<int> neighborhoodPtr=gData.neighborhoodPtr;
			MPI_Send(&ack, 1, MPI_INT, 0, ackTag, MPI_COMM_WORLD);
			MPI_Recv(&globalNumPoints, 1, MPI_INT, 0, globalNumPointsTag, MPI_COMM_WORLD, &status);
			MPI_Recv(&sizeNeighborhoodList, 1, MPI_INT, 0, sizeNeighborhoodListTag, MPI_COMM_WORLD, &status);
			std::tr1::shared_ptr<int> neighborhood(new int[sizeNeighborhoodList],Deleter<int>());
			MPI_Recv(gIds.get(), numPoints, MPI_INT, 0, idsTag, MPI_COMM_WORLD, &status);
			MPI_Recv(g.get(), dimension*numPoints, MPI_DOUBLE, 0, coordinatesTag, MPI_COMM_WORLD, &status);
			MPI_Recv(cellVolume.get(), numPoints, MPI_DOUBLE, 0,volumeTag, MPI_COMM_WORLD, &status);
			MPI_Recv(neighborhood.get(), sizeNeighborhoodList, MPI_INT, 0, neighborhoodTag, MPI_COMM_WORLD, &status);
			MPI_Recv(neighborhoodPtr.get(), numPoints, MPI_INT, 0, neighborhoodPtrTag, MPI_COMM_WORLD, &status);

			gData.dimension = dimension;
			gData.globalNumPoints = globalNumPoints;
			gData.numPoints = numPoints;
			gData.sizeNeighborhoodList = sizeNeighborhoodList;
			gData.neighborhood=neighborhood;

			gridData = gData;

		}
		else if (numPoints == 0){
			MPI_Send(&ack, 1, MPI_INT, 0, ackTag, MPI_COMM_WORLD);
		}
		else{
			MPI_Finalize();
			exit(1);
		}

	    MPI_Recv(&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	    if (ack < 0){
	      MPI_Finalize();
	      exit(1);
	    }

	}
	return gridData;
}


shared_ptr< std::set<int> > constructParallelDecompositionFrameSet(PdGridData& decomp, double horizon) {
	shared_ptr<double> xPtr = decomp.myX;
	shared_ptr<int> gIdsPtr = decomp.myGlobalIDs;
	int numCells = decomp.numPoints;
	const Coordinates c(xPtr,numCells);

	int numAxes=3;
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y, PdNeighborhood::Z};
	std::vector<std::tr1::shared_ptr<int> > sortedMaps(numAxes);
	for(int j=0;j<numAxes;j++){

		Coordinates::SortComparator compare = c.getSortComparator(labels[j]);
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
			const Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const Coordinates::SearchIterator end=start+numCells;
			Coordinates::SearchIterator lub = std::upper_bound(start,end,value);
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
			const Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const Coordinates::SearchIterator end=start+numCells;
			Coordinates::SearchIterator glb = std::upper_bound(start,end,value);
			frameSetPtr->insert(glb.mapIterator(),glb.mapEnd());
		}
	}

	return frameSetPtr;

}

} // namespace PdQuickGrid3D

