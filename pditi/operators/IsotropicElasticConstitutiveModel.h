/*
 * IsotropicElasticConstitutiveModel.h
 *
 *  Created on: Jun 23, 2010
 *      Author: jamitch
 */

#ifndef ISOTROPICELASTICCONSTITUTIVEMODEL_H_
#define ISOTROPICELASTICCONSTITUTIVEMODEL_H_
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


class IsotropicElasticConstitutiveModel : public ConstitutiveModel {

public:
	~IsotropicElasticConstitutiveModel() {}

	IsotropicElasticConstitutiveModel(IsotropicHookeSpec matSpec);
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
			double* fInternalOverlapPtr,
			const int*  localNeighborList,
			int numOwnedPoints,
			vector< TemporalField<double> > temporalFields = vector< TemporalField<double> >()
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
			double horizon
	);

private:
	PdImp::IsotropicHookeSpec matSpec;

};

}
#endif /* ISOTROPICELASTICCONSTITUTIVEMODEL_H_ */
