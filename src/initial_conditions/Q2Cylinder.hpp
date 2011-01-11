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


namespace InitialConditionsNS {

class Q2Cylinder : public InitialCondition {
public:
	virtual ~Q2Cylinder() {}
	Q2Cylinder(double _vr0, double _vr1, double _vz0, double _z0, double _a, const VectorUtilsNS::Vector3D& _center);
	virtual void apply(const Epetra_Vector& x, Epetra_Vector& u, Epetra_Vector& v);

private:
	double vr0, vr1, vz0, z0, a;
	VectorUtilsNS::Vector3D center;
};

}  // namespace InitialConditionsNS



#endif /* Q2CYLINDER_HPP_ */
