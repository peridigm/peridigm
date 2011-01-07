/*
 * DirichletBc.h
 *
 *  Created on: Jun 21, 2010
 *      Author: jamitch
 */

#ifndef STAGE_DIRICHLETBC_H_
#define STAGE_DIRICHLETBC_H_
#include "Field.h"

namespace PdImp {

class DirichletBcSpec;

class StageDirichletBc {
public:
	virtual ~StageDirichletBc() {}
	virtual const DirichletBcSpec& getSpec() const = 0;
	virtual void applyHomogeneousForm(Field_NS::Field<double>& residual) const = 0;
	virtual void applyKinematics(double lambda, Field_NS::Field<double>& displacement) const = 0;
};

}


#endif /* STAGE_DIRICHLETBC_H_ */
