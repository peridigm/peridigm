/*
 * Q2Cylinder.hpp
 *
 *  Created on: Jan 4, 2011
 *      Author: jamitch
 */

#ifndef Q2CYLINDER_HPP_
#define Q2CYLINDER_HPP_

#include "InitialCondition.hpp"
#include "VectorUtils.h"

namespace PeridigmNS {

namespace InitialConditionsNS {

class Q2Cylinder : public InitialCondition {
public:
	virtual ~Q2Cylinder() {}
	Q2Cylinder(const Teuchos::RCP<Teuchos::ParameterList>& peridigmParams);
	virtual void apply() {}

private:
	double vr0, vr1, vz0, z0, a;
	VectorUtilsNS::Vector3D center;
};

}  // namespace InitialConditionsNS

} // namespace PeridigmNS

#endif /* Q2CYLINDER_HPP_ */
