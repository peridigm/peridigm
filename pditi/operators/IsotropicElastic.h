/*
 * PimpIsotropicElastic.h
 *
 *  Created on: Apr 3, 2010
 *      Author: awesome
 */

#ifndef ISOTROPICELASTIC_H_
#define ISOTROPICELASTIC_H_
#include "ConstitutiveModel.h"

namespace PdImp {

class IsotropicElastic : public ConstitutiveModel {

public:
	virtual ~IsotropicElastic() {}
	virtual void computeInternalForce(Field<double>& displacement, Field<double>& velocity, double delta_T, Field<double>& fInt) = 0;

};

}
#endif /* ISOTROPICELASTIC_H_ */
