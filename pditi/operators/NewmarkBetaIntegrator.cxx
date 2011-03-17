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
		TemporalField<double>& displacement,
		TemporalField<double>& velocity,
		TemporalField<double>& acceleration,
		double delta_T
) const
{

	FieldSpec::FieldStep NP1 = FieldSpec::STEP_NP1;
	FieldSpec::FieldStep N = FieldSpec::STEP_N;

	const double* uN    = displacement.getField(N).get();
	const double* uNp1  = displacement.getField(NP1).get();
	const double* vN    = velocity.getField(N).get();
	double* vNp1  = velocity.getField(NP1).get();
	const double* aN    = acceleration.getField(N).get();
	double* aNp1  = acceleration.getField(NP1).get();

	/*
	 * Note that 'get_size' is not the same as 'get_num_points'
	 */
	std::size_t n = displacement.getField(N).get_size();
	integrateStep(n,alpha,beta,uN,uNp1,vN,vNp1,aN,aNp1,delta_T);

}

void NewmarkBetaIntegrator::advanceState
(
		TemporalField<double>& displacement,
		TemporalField<double>& velocity,
		TemporalField<double>& acceleration
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
