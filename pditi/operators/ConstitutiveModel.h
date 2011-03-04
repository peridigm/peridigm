/*
 * PimpConstitutiveModel.h
 *
 *  Created on: Apr 3, 2010
 *      Author: awesome
 */

#ifndef CONSTITUTIVEMODEL_H_
#define CONSTITUTIVEMODEL_H_
#include "Field.h"
#include <vector>

using std::vector;
using Field_NS::FieldSpec;
using Field_NS::Field;
using Field_NS::TemporalField;

namespace PdITI {

class ConstitutiveModel {

public:

	virtual ~ConstitutiveModel() {}
	/*
	 * Bond variables integrated with time; ie Temporal fields associated with each bond
	 */
	virtual vector<FieldSpec> registerTemporalBondVariables() const = 0;

	/*
	 * COMPUTE INTERNAL FORCE
	 */
	virtual void computeInternalForce
	(
			const double* xOverlapPtr,
			const double* yOverlapPtr,
			const double* mOwned,
			const double* volumeOverlapPtr,
			const double* dilatationOwned,
			const double* bondDamage,
			const double* dsf,
			double* fInternalOverlapPtr,
			const int*  localNeighborList,
			int numOwnedPoints,
			vector< TemporalField<double> > temporalFields = vector< TemporalField<double> >()
	) = 0;

	/*
	 * Compute 3x3 tangent stiffness associated with I, P, Q
	 */
	virtual double* kIPQ3x3(
			int I,
			int P,
			int Q,
			const double* x,
			const double* u,
			const double* volume,
			const double* m,
			const double *dilatation,
			double *k,
			double horizon,
			double dsf_I=1.0
	) = 0;

};

}
#endif /* CONSTITUTIVEMODEL_H_ */
