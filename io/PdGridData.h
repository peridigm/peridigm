/*
 * PdGridData.h
 *
 *  Created on: Nov 20, 2009
 *      Author: jamitch
 */

#ifndef PDGRIDDATA_H_
#define PDGRIDDATA_H_
#include <tr1/memory>
using std::tr1::shared_ptr;

struct Zoltan_Struct;

typedef struct Data {
	int dimension;
	int globalNumPoints;
	int numPoints;
	int sizeNeighborhoodList;
	int numExport;
	bool unPack;
	std::tr1::shared_ptr<int> myGlobalIDs;
	std::tr1::shared_ptr<double> myX;
	std::tr1::shared_ptr<double> cellVolume;
	std::tr1::shared_ptr<int> neighborhood;
	std::tr1::shared_ptr<int> neighborhoodPtr;
	std::tr1::shared_ptr<char> exportFlag;
	std::tr1::shared_ptr<struct Zoltan_Struct> zoltanPtr;
	Data() : dimension(-1), globalNumPoints(-1), numPoints(-1), sizeNeighborhoodList(-1), numExport(0) {}
	Data(int d, int numPoints, int myNumPts) : dimension(d), globalNumPoints(numPoints), numPoints(myNumPts) {}
} PdGridData;

#endif /* PDGRIDDATA_H_ */
