/*
 * PdQuickGrid3D.h
 *
 *  Created on: Nov 9, 2009
 *      Author: jamitch
 */

#ifndef PDQUICKGRID3D_H_
#define PDQUICKGRID3D_H_

#include "PdQuickGrid.h"
#include "PdGridData.h"
#include <tr1/memory>
#include <set>
using std::tr1::shared_ptr;
using std::set;


namespace PdQuickGrid {

enum Pd_MPI_TAGS {commCreate=9,commDo=10};

PdGridData getDiscretization(int rank, PdQuickGridMeshGenerationIterator &meshGenerationIterator);

shared_ptr< std::set<int> > constructParallelDecompositionFrameSet(PdGridData& decomp, double horizon);

class PointBB {
public:
	PointBB (const double* center, double horizon);
	inline double get_xMin() const { return xMin; }
	inline double get_yMin() const { return yMin; }
	inline double get_zMin() const { return zMin; }
	inline double get_xMax() const { return xMax; }
	inline double get_yMax() const { return yMax; }
	inline double get_zMax() const { return zMax; }

private:
	double xMin, xMax, yMin, yMax, zMin, zMax;
};

} // namespace PdQuickGrid3D

#endif /* PDQUICKGRID3D_H_ */
