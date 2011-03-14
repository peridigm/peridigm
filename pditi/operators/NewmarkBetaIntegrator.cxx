/*
 * NewmarkBetaIntegrator.cxx
 *
 *  Created on: Jun 16, 2010
 *      Author: jamitch
 */
#include "NewmarkBetaIntegrator.h"

namespace PdImp {

using Field_NS::Field;
using Field_NS::FieldSpec;
using Field_NS::TemporalField;

void NewmarkBetaIntegrator::integrateStep
(
		Field_NS::TemporalField<double>& displacement,
		Field_NS::TemporalField<double>& velocity,
		Field_NS::TemporalField<double>& acceleration,
		double delta_T
) const
{

	FieldSpec::FieldStep NP1 = FieldSpec::STEP_NP1;
	FieldSpec::FieldStep N = FieldSpec::STEP_N;

	Pd_shared_ptr_Array<double> uN    = displacement.getField(N).getArray();
	Pd_shared_ptr_Array<double> uNp1  = displacement.getField(NP1).getArray();
	Pd_shared_ptr_Array<double> vN    = velocity.getField(N).getArray();
	Pd_shared_ptr_Array<double> vNp1  = velocity.getField(NP1).getArray();
	Pd_shared_ptr_Array<double> aN    = acceleration.getField(N).getArray();
	Pd_shared_ptr_Array<double> aNp1  = acceleration.getField(NP1).getArray();

	std::size_t n = uN.getSize();
	integrateStep(n,alpha,beta,uN.get(),uNp1.get(),vN.get(),vNp1.get(),aN.get(),aNp1.get(),delta_T);


}

void NewmarkBetaIntegrator::advanceState
(
		Field_NS::TemporalField<double>& displacement,
		Field_NS::TemporalField<double>& velocity,
		Field_NS::TemporalField<double>& acceleration
)
{
	displacement.advanceStep();
	velocity.advanceStep();
	acceleration.advanceStep();
}

void NewmarkBetaIntegrator::integrateStep
(
		std::size_t n,
		double alpha,
		double beta,
		const double *uN,
		const double *uNp1,
		const double* vN,
		double* vNp1,
		const double* aN,
		double* aNp1,
		double dt
)
{
	double c2 = 1.0/2.0/beta - 1.0;
	double c0 = beta * dt;
	double c1 = c0 * dt;
	for(std::size_t j=0;j<n;j++,uN++,uNp1++,vN++,vNp1++,aN++,aNp1++){
		*aNp1 = (*uNp1-*uN)/c1 -*vN/c0-*aN*c2;
		*vNp1 = *vN + *aN*(1-alpha)*dt + *aNp1*alpha*dt;
	}

}

}
