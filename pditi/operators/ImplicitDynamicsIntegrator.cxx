/*
 * ImplicitDynamicsResidual.cxx
 *
 *  Created on: Apr 3, 2010
 *      Author: awesome
 */

#include "ImplicitDynamicsIntegrator.h"


namespace PdImp {

using Field_NS::Field;
using Field_NS::FieldSpec;
using Field_NS::TemporalField;
using PdITI::PdITI_Operator;

void ImplicitDynamicsIntegrator::computeResidual
(
		TemporalField<double>& displacement,
		TemporalField<double>& velocity,
		TemporalField<double>& acceleration,
		double delta_T,
		TemporalField<double>& residual
) const
{


	/**
	 * Initialize residual to zero
	 * Compute internal force for step "N+1"
	 */
	double zero(0.0);
	FieldSpec::FieldStep NP1 = FieldSpec::STEP_NP1;
	FieldSpec::FieldStep N = FieldSpec::STEP_N;
	Field<double> rNp1Field = residual.getField(NP1);
	rNp1Field.set(zero);
	fIntOperator->computeInternalForce(displacement.getField(NP1),residual.getField(NP1));

	double* fN    = residual.getField(N).get();
	double* resid = residual.getField(NP1).get();
	double* uN    = displacement.getField(N).get();
	double* uNp1  = displacement.getField(NP1).get();
	double* vN    = velocity.getField(N).get();

	double beta = timeIntegrator.getBeta();
	size_t n = rNp1Field.get_size();
	computeResidual(rho.getValue(),n,beta,fN,resid,uN,uNp1,vN,delta_T);
}

void ImplicitDynamicsIntegrator::computeResidual
(
		double rho,
		std::size_t n,
		double beta,
		const double *fN,
		double *fNp1,
		const double *uN,
		const double *uNp1,
		const double* vN,
		double dt
){
	double dt2 = dt*dt;
	double bdt2 = beta*dt2;
	double c = (.5-beta)*dt2;
	for(std::size_t j=0;j<n;j++,fN++,fNp1++,uN++,uNp1++,vN++)
		*fNp1 = -( rho * (*uNp1) - bdt2*(*fNp1) - rho * (*uN + *vN*dt ) - c * (*fN) );
}

}
