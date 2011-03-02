/*
 * ImplicitLinearDynamicsResidual.cxx
 *
 *  Created on: Jun 17, 2010
 *      Author: jamitch
 */

#include "Field.h"
#include "ImplicitLinearDynamicsIntegrator.h"
#include "PdImpOperatorUtilities.h"


namespace PdImp {

using Field_NS::Field;
using Field_NS::FieldSpec;
using Field_NS::TemporalField;

void ImplicitLinearDynamicsIntegrator::computeResidual
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
	 * Compute right hand side vector due to inertial terms
	 */
	double zero(0.0);
	FieldSpec::FieldStep NP1 = FieldSpec::STEP_NP1;
	FieldSpec::FieldStep N = FieldSpec::STEP_N;
	PdITI::SET(residual.getField(NP1).getArray().get(),residual.getField(NP1).getArray().end(),zero);

	Pd_shared_ptr_Array<double> uN    = displacement.getField(N).getArray();
	Pd_shared_ptr_Array<double> vN    = velocity.getField(N).getArray();
	Pd_shared_ptr_Array<double> aN    = acceleration.getField(N).getArray();
	Pd_shared_ptr_Array<double> rNp1  = residual.getField(NP1).getArray();

	double beta = timeIntegrator.getBeta();
	computeResidual(rho.getValue(),rNp1.getSize(),beta,rNp1.get(),uN.get(),vN.get(),aN.get(),delta_T);

}

void ImplicitLinearDynamicsIntegrator::computeResidual
(
		double rho,
		std::size_t n,
		double beta,
		double* rNp1,
		const double *uN,
		const double* vN,
		const double* aN,
		double dt
){
	double c2 = 1.0/2.0/beta - 1.0;
	double c0 = beta * dt;
	double c1 = c0 * dt;
	for(std::size_t j=0;j<n;j++,rNp1++,uN++,vN++,aN++)
		*rNp1 = rho * (*uN/c1 + *vN/c0 + *aN*c2);
}

}
