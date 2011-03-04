/*
 * IsotropicElastic_No_DSF.h
 *
 *  Created on: Mar 3, 2011
 *      Author: jamitch
 */

#ifndef ISOTROPICELASTIC_NO_DSF_H_
#define ISOTROPICELASTIC_NO_DSF_H_

#include <vector>
#include "Field.h"
#include "PdImpMaterials.h"
#include "ConstitutiveModel.h"

using std::vector;
using Field_NS::FieldSpec;
using Field_NS::Field;
using Field_NS::TemporalField;
using PdImp::IsotropicHookeSpec;

namespace PdITI {


class IsotropicElastic_No_DSF : public ConstitutiveModel {

public:
	~IsotropicElastic_No_DSF() {}

	IsotropicElastic_No_DSF(IsotropicHookeSpec matSpec);
	/*
	 * Bond variables integrated with time; ie Temporal fields associated with each bond
	 */
	vector<FieldSpec> registerTemporalBondVariables() const;

	/*
	 * COMPUTE INTERNAL FORCE
	 */
	void computeInternalForce
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
	);

	/*
	 * Calculation includes DSF factor
	 */
	void computeInternalForceLinearElastic
	(
			const double* xOverlap,
			const double* yOverlap,
			const double* mOwned,
			const double* volumeOverlap,
			const double* dilatationOwned,
			const double* bondDamage,
			const double* dsf,
			double* fInternalOverlap,
			const int*  localNeighborList,
			int numOwnedPoints,
			double BULK_MODULUS,
			double SHEAR_MODULUS
	);

	/*
	 * Compute 3x3 tangent stiffness associated with I, P, Q
	 */
	double* kIPQ3x3(
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
	);

private:
	PdImp::IsotropicHookeSpec matSpec;

};

}

#endif /* ISOTROPICELASTIC_NO_DSF_H_ */
