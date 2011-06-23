/*
 * PdGridData.h
 *
 *  Created on: Nov 20, 2009
 *      Author: jamitch
 */

#ifndef QUICKGRIDDATA_H_
#define QUICKGRIDDATA_H_
#include <tr1/memory>
using std::tr1::shared_ptr;

struct Zoltan_Struct;

namespace QUICKGRID {

typedef struct Data {
	int dimension;
	int globalNumPoints;
	int numPoints;
	int sizeNeighborhoodList;
	int numExport;
	bool unPack;
	shared_ptr<int> myGlobalIDs;
	shared_ptr<double> myX;
	shared_ptr<double> cellVolume;
	shared_ptr<int> neighborhood;
	shared_ptr<int> neighborhoodPtr;
	shared_ptr<char> exportFlag;
	shared_ptr<struct Zoltan_Struct> zoltanPtr;
	Data() : dimension(-1), globalNumPoints(-1), numPoints(-1), sizeNeighborhoodList(-1), numExport(0) {}
	Data(int d, int numPoints, int myNumPts) : dimension(d), globalNumPoints(numPoints), numPoints(myNumPts) {}
} QuickGridData;

}

#endif /* QUICKGRIDDATA_H_ */
