/*
 * ImplicitDynamicsResidual.cxx
 *
 *  Created on: Apr 3, 2010
 *      Author: awesome
 */

#include "ImplicitDynamicsIntegrator.h"
#include "PdImpOperatorUtilities.h"


namespace PdImp {

using Field_NS::Field;
using Field_NS::FieldSpec;
using Field_NS::TemporalField;

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
	PdImp::SET(residual.getField(NP1).getArray().get(),residual.getField(NP1).getArray().end(),zero);
	fIntOperator->computeInternalForce(displacement.getField(NP1),residual.getField(NP1));

	Pd_shared_ptr_Array<double> fN    = residual.getField(N).getArray();
	Pd_shared_ptr_Array<double> resid = residual.getField(NP1).getArray();
	Pd_shared_ptr_Array<double> uN    = displacement.getField(N).getArray();
	Pd_shared_ptr_Array<double> uNp1  = displacement.getField(NP1).getArray();
	Pd_shared_ptr_Array<double> vN    = velocity.getField(N).getArray();

	double beta = timeIntegrator.getBeta();
	computeResidual(rho.getValue(),fN.getSize(),beta,fN.get(),resid.get(),uN.get(),uNp1.get(),vN.get(),delta_T);
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
