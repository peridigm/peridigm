/*
 * DirichletBcSpec.h
 *
 *  Created on: Apr 26, 2010
 *      Author: jamitch
 */

#ifndef DIRICHLETBCSPEC_H_
#define DIRICHLETBCSPEC_H_
#include "Array.h"
#include "StageDirichletBc.h"
#include <vector>

using std::vector;

namespace PdImp {

class StageFunction;

class DirichletBcSpec {

public:
	enum ComponentLabel {X=0,Y,Z};
	virtual ~DirichletBcSpec() {}
	virtual vector<ComponentLabel> getComponents() const = 0;
	virtual vector< vector<double> > getUnitDirections() const = 0;
	virtual const UTILITIES::Array<int>& getPointIds() const = 0;
	virtual const StageDirichletBc& getStageDirichletBc(const StageFunction& stageFunction) const = 0;
};

}

#endif /* DIRICHLETBCSPEC_H_ */
