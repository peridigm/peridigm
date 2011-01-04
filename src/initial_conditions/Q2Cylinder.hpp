/*
 * Q2Cylinder.hpp
 *
 *  Created on: Jan 4, 2011
 *      Author: jamitch
 */

#ifndef Q2CYLINDER_HPP_
#define Q2CYLINDER_HPP_

#include "InitialCondition.hpp"

namespace PeridigmNS {

namespace InitialConditionsNS {

class Q2Cylinder : public InitialCondition {
public:
	virtual ~Q2Cylinder() {}
	Q2Cylinder(const Teuchos::RCP<Teuchos::ParameterList>& peridigmParams) {}
	virtual void apply() {}

private:
	double z0, a, vr0, vr1, vz0;
};

}  // namespace InitialConditionsNS

} // namespace PeridigmNS

#endif /* Q2CYLINDER_HPP_ */
