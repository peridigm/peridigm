/*
 * StageComponentDirichletBc.cxx
 *
 *  Created on: Jun 21, 2010
 *      Author: jamitch
 */
#include "StageComponentDirichletBc.h"
#include "ComponentDirichletBcSpec.h"
#include "StageFunction.h"
#include "PdITI_Utilities.h"

namespace PdImp {

using Field_NS::Field;

StageComponentDirichletBc::StageComponentDirichletBc(const ComponentDirichletBcSpec& spec, const StageFunction& stageFunction)
: spec(spec), stageFunction(stageFunction)
{}

void StageComponentDirichletBc::imprint_bc(Field_NS::Field<char>& maskField) const {

	char* b = maskField.get();
	char mask = spec.get_mask();
	const UTILITIES::Array<int>& ids = spec.getPointIds();
	/*
	 * NOTE bitwise 'or' since we are essentially concatenating boundary conditions
	 */
	for(const int *idsPtr=ids.get();idsPtr!=ids.end();idsPtr++, b++){
		*b |= mask;
	}
}


void StageComponentDirichletBc::applyHomogeneousForm(Field<double>& residualField) const {
	vector< vector<double> > dirs = spec.getUnitDirections();
	const UTILITIES::Array<int>& ids = spec.getPointIds();
	int numDirs = dirs.size();
	double vec[3];
	double *residual = residualField.get();
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


void StageComponentDirichletBc::applyKinematics(double lambda, Field<double>& displacement) const {
	vector< vector<double> > dirs = spec.getUnitDirections();
	const UTILITIES::Array<int>& ids = spec.getPointIds();
	int numDirs = dirs.size();
	double vec[3];
	double *uHead = displacement.get();
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
