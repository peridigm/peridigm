/*
 * PimpImplicitResidual.h
 *
 *  Created on: Apr 3, 2010
 *      Author: awesome
 */

#ifndef IMPLICIT_DYNAMICS_RESIDUAL_H_
#define IMPLICIT_DYNAMICS_RESIDUAL_H_
#include "PdImpOperator.h"
#include "PdImpMaterials.h"
#include "NewmarkBetaIntegrator.h"
#include <tr1/memory>
#include "Field.h"


namespace PdImp {

using Field_NS::TemporalField;

class ImplicitDynamicsIntegrator {
public:
	/**
	 * Temporary Hack Until I can get the complete abstraction finished;
	 * Future ImplicitResidual will take std::vector< std::tr1::shared_ptr<ContitutiveModel> >
	 * @param alpha=1/2
	 * @param beta=1/4
	 * Above parameters create the "average acceleration" method
	 */
	ImplicitDynamicsIntegrator(PdImp::MassDensity val, std::tr1::shared_ptr<PdImpOperator> fIntOp) : timeIntegrator(0.5, 0.25), rho(val), fIntOperator(fIntOp) {}
	virtual ~ImplicitDynamicsIntegrator() {}
	void computeResidual
	(
			TemporalField<double>& displacement,
			TemporalField<double>& velocity,
			TemporalField<double>& acceleration,
			double delta_T,
			TemporalField<double>& residual
	) const;

	static void computeResidual
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
	);

private:
	//	std::tr1::shared_ptr<ConstitutiveModel> fInt;
	//	std::tr1::shared_ptr<ExternalForce> fExt;

	NewmarkBetaIntegrator timeIntegrator;
	MassDensity rho;
//	std::tr1::shared_ptr<PdImpOperator> fIntOperator;
};

}
#endif /* IMPLICIT_DYNAMICS_RESIDUAL_H_ */
