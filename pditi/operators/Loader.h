/*
 * Loader.h
 *
 *  Created on: Jun 18, 2010
 *      Author: jamitch
 */

#ifndef LOADER_H_
#define LOADER_H_

#include "vtk/Field.h"

namespace PdImp {

class Loader {
public:
	virtual ~Loader() {}
	virtual void computeOwnedExternalForce(double lambda, Field_NS::Field<double> force) const=0;

};

}


#endif /* LOADER_H_ */
