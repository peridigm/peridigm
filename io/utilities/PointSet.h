/*
 * PointSet.h
 *
 *  Created on: Mar 11, 2011
 *      Author: wow
 */

#ifndef POINTSET_H_
#define POINTSET_H_

#include "Array.h"

namespace UTILITIES {

class PointSet : public Array<double> {

public:

	PointSet() : Array<double>(), num_points(0) {}

	PointSet(size_t num_points) : Array<double>(3*num_points), num_points(num_points) {}

	size_t get_num_points() const { return num_points; }

private:

	size_t num_points;

};

}
#endif /* POINTSET_H_ */
