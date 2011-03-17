/*
 * BodyLoadSpec.cxx
 *
 *  Created on: Jun 18, 2010
 *      Author: jamitch
 */
#include "BodyLoadSpec.h"
#include "PdITI_Utilities.h"

namespace PdImp {

using Field_NS::Field;

BodyLoadSpec::BodyLoadSpec(double unitDirection[3], const UTILITIES::Array<int> pointIds)
:
		localIds(pointIds)
{
	u[0]=unitDirection[0];
	u[1]=unitDirection[1];
	u[2]=unitDirection[2];
}

shared_ptr<Loader> BodyLoadSpec::getStageLoader(const StageFunction& stageFunction) const {
	return shared_ptr<Loader>(new BodyLoader(*this,stageFunction));
}

BodyLoadSpec::BodyLoader::BodyLoader(const BodyLoadSpec& spec, const StageFunction& stageFunction)
:
		stageFunction(stageFunction), localIds(spec.localIds)
{
	u[0] = spec.u[0];
	u[1] = spec.u[1];
	u[2] = spec.u[2];
}

void BodyLoadSpec::BodyLoader::computeOwnedExternalForce(double lambda, Field<double> force) const {
	double bF[3];
	double *fHead = force.get();
	double *f;
	double magnitude = stageFunction.value(lambda);
	for(const int *p = localIds.get();p!=localIds.end();p++){

		int id = *p;
		/*
		 * bF = u
		 */
		PdITI::COPY(u,u+3,bF);

		/*
		 * bF = b * bF = b * u
		 */
		PdITI::SCALE_BY_VALUE(bF,bF+3,magnitude);

		f = fHead + 3 * id;

		/*
		 * f = f + bF
		 */
		PdITI::SUMINTO(bF,bF+3,f);

	}

}
}
