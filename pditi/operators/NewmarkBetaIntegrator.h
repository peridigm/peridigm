/*
 * NewmarkBetaIntegrator.h
 *
 *  Created on: Jun 16, 2010
 *      Author: jamitch
 */

#ifndef NEWMARKBETAINTEGRATOR_H_
#define NEWMARKBETAINTEGRATOR_H_
#include "vtk/Field.h"

namespace PdImp {


class NewmarkBetaIntegrator {

public:
	NewmarkBetaIntegrator(double alpha, double beta) : alpha(alpha), beta(beta) {}
	void integrateStep
	(
			Field_NS::TemporalField<double>& displacement,
			Field_NS::TemporalField<double>& velocity,
			Field_NS::TemporalField<double>& acceleration,
			double delta_T
	) const;

	static void advanceState
	(
			Field_NS::TemporalField<double>& displacement,
			Field_NS::TemporalField<double>& velocity,
			Field_NS::TemporalField<double>& acceleration
	);

	double getAlpha() const { return alpha; }
	double getBeta() const { return beta; }

	static void integrateStep
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
	);


private:
	/**
	 * Integration parameters for Newmark-Beta method
	 */
	double alpha, beta;
};

}
#endif /* NEWMARKBETAINTEGRATOR_H_ */
