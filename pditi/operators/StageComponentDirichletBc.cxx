/*
 * StageComponentDirichletBc.cxx
 *
 *  Created on: Jun 21, 2010
 *      Author: jamitch
 */
#include "StageComponentDirichletBc.h"
#include "ComponentDirichletBcSpec.h"
#include "StageFunction.h"
#include "PdImpOperatorUtilities.h"

namespace PdImp {


StageComponentDirichletBc::StageComponentDirichletBc(const ComponentDirichletBcSpec& spec, const StageFunction& stageFunction)
: spec(spec), stageFunction(stageFunction)
{}

void StageComponentDirichletBc::applyHomogeneousForm(Field_NS::Field<double>& residualField) const {
	vector< vector<double> > dirs = spec.getUnitDirections();
	const Pd_shared_ptr_Array<int>& ids = spec.getPointIds();
	int numDirs = dirs.size();
	double vec[3];
	double *residual = residualField.getArray().get();
	for(const int *idsPtr=ids.get();idsPtr!=ids.end();idsPtr++){
		/*
		 * NOTE:
		 * This implementation assumes all coordinate directions are orthogonal to each other
		 */
		double *r = residual+3*(*idsPtr);
		for(int i=0;i<numDirs;i++){
			vector<double> u = dirs[i];
			double rn = PdITI::DOT(&u[0],r);
			// Set vec=u
			PdITI::COPY(&u[0],&u[0]+3,vec);

			// Set vec = rn*vec = rn*u
			PdITI::SCALE_BY_VALUE(vec,vec+3,rn);

			// Set R = R - rn*vec
			PdITI::SUBTRACTINTO(vec,vec+3,r);

		}
	}
}


void StageComponentDirichletBc::applyKinematics(double lambda, Field_NS::Field<double>& displacement) const {
	vector< vector<double> > dirs = spec.getUnitDirections();
	const Pd_shared_ptr_Array<int>& ids = spec.getPointIds();
	int numDirs = dirs.size();
	double vec[3];
	double *uHead = displacement.getArray().get();
	double val = stageFunction.value(lambda);
	for(const int *idsPtr=ids.get();idsPtr!=ids.end();idsPtr++){
		/*
		 * NOTE:
		 * This implementation assumes all coordinate directions are orthogonal to each other
		 */
		double *u = uHead+3*(*idsPtr);
		for(int i=0;i<numDirs;i++){

			vector<double> n = dirs[i];
			double nx=n[0], ny = n[1], nz = n[2];
			double ux=*(u), uy = *(u+1), uz = *(u+2);
			double un = ux*nx + uy*ny + uz*nz;
			*(u+0) = *(u+0) - un * nx + val * nx;
			*(u+1) = *(u+1) - un * ny + val * ny;
			*(u+2) = *(u+2) - un * nz + val * nz;

		}
	}
}

}
