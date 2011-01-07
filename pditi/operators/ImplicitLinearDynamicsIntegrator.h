/*
 * ImplicitLinearDynamicsIntegrator.h
 *
 *  Created on: Jun 17, 2010
 *      Author: jamitch
 */
#include "PdImpMaterials.h"
#include "NewmarkBetaIntegrator.h"
#include "Field.h"
#include <tr1/memory>

namespace PdImp {

using Field_NS::TemporalField;

class ImplicitLinearDynamicsIntegrator {
public:
	/**
	 * Temporary Hack Until I can get the complete abstraction finished;
	 * Future ImplicitResidual will take std::vector< std::tr1::shared_ptr<ContitutiveModel> >
	 * @param alpha=1/2
	 * @param beta=1/4
	 * Above parameters create the "average acceleration" method
	 */
	ImplicitLinearDynamicsIntegrator(PdImp::MassDensity val) : timeIntegrator(0.5, 0.25), rho(val) {}
	virtual ~ImplicitLinearDynamicsIntegrator() {}

	const NewmarkBetaIntegrator& getNewmarkIntegrator() const { return timeIntegrator; }
	const MassDensity& getDensity() const { return rho; }

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
			double* rNp1,
			const double *uN,
			const double* vN,
			const double* aN,
			double dt
	);

private:

	const NewmarkBetaIntegrator timeIntegrator;
	const MassDensity rho;
};

}
