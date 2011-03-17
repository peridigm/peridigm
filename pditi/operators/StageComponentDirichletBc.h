/*
 * StageComponentDirichletBc.h
 *
 *  Created on: Jun 21, 2010
 *      Author: jamitch
 */

#ifndef STAGECOMPONENTDIRICHLETBC_H_
#define STAGECOMPONENTDIRICHLETBC_H_

#include "vtk/Field.h"
#include "StageDirichletBc.h"
#include "StageFunction.h"
#include "ComponentDirichletBcSpec.h"

namespace PdImp {


class StageComponentDirichletBc : public StageDirichletBc {

public:
	StageComponentDirichletBc(const ComponentDirichletBcSpec& spec, const StageFunction& stageFunction);
	const DirichletBcSpec& getSpec() const { return spec; }
	void applyHomogeneousForm(Field_NS::Field<double>& residual) const;
	void applyKinematics(double lambda, Field_NS::Field<double>& displacement) const;

private:
	const ComponentDirichletBcSpec spec;
	const StageFunction stageFunction;
};

}
#endif /* STAGECOMPONENTDIRICHLETBC_H_ */
