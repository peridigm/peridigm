/*
 * ComponentDirichletBcSpec.h
 *
 *  Created on: Jun 21, 2010
 *      Author: jamitch
 */

#ifndef COMPONENTDIRICHLETBCSPEC_H_
#define COMPONENTDIRICHLETBCSPEC_H_
#include "DirichletBcSpec.h"
#include "StageFunction.h"
#include <vector>

using std::vector;

namespace PdImp {

class ComponentDirichletBcSpec : public DirichletBcSpec {
public:
	ComponentDirichletBcSpec(const vector<DirichletBcSpec::ComponentLabel>& components, const Pd_shared_ptr_Array<int>& pointIds);
	~ComponentDirichletBcSpec();
	vector<DirichletBcSpec::ComponentLabel> getComponents() const;
	vector< vector<double> > getUnitDirections() const;
	const Pd_shared_ptr_Array<int>& getPointIds() const;
	const StageDirichletBc& getStageDirichletBc(const StageFunction& stageFunction) const;
	static ComponentDirichletBcSpec getAllComponents(const Pd_shared_ptr_Array<int>& pointIds);

private:
	vector<DirichletBcSpec::ComponentLabel> components;
	vector< vector<double> > directions;
	/*
	 * Set of point Ids -- these Ids should be relevant to the fields passed in to the constructor
	 */
	const Pd_shared_ptr_Array<int> localIds;

};



}
#endif /* COMPONENTDIRICHLETBCSPEC_H_ */
