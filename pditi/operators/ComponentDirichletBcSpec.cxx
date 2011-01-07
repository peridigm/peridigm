/*
 * ComponentDirichletBcSpec.cxx
 *
 *  Created on: Jun 21, 2010
 *      Author: jamitch
 */

#include "ComponentDirichletBcSpec.h"
#include "StageComponentDirichletBc.h"

namespace PdImp {

using Field_NS::Field;

ComponentDirichletBcSpec::
ComponentDirichletBcSpec(const vector<DirichletBcSpec::ComponentLabel>& components,const Pd_shared_ptr_Array<int>& pointIds)
:
components(components),
directions(components.size()),
localIds(pointIds)
{
	vector<ComponentLabel>::iterator compIter = this->components.begin();
	vector< vector<double> >::iterator dirIter = directions.begin();

	for(;compIter!=this->components.end();compIter++,dirIter++){
		*dirIter=vector<double>(3);
		std::fill(dirIter->begin(),dirIter->end(),0.0);

		switch(*compIter){

		case X:
			(*dirIter)[0]=1.0;
			break;
		case Y:
			(*dirIter)[1]=1.0;
			break;

		case Z:
			(*dirIter)[2]=1.0;
			break;
		}

	}

}


ComponentDirichletBcSpec::~ComponentDirichletBcSpec() {}

vector<DirichletBcSpec::ComponentLabel> ComponentDirichletBcSpec::getComponents() const { return components; }


vector< vector<double> >ComponentDirichletBcSpec::getUnitDirections() const { return directions; }

const Pd_shared_ptr_Array<int>& ComponentDirichletBcSpec::getPointIds() const { return localIds; }

const StageDirichletBc& ComponentDirichletBcSpec::getStageDirichletBc(const StageFunction& stageFunction) const {
	return StageComponentDirichletBc(*this,stageFunction);
}

ComponentDirichletBcSpec ComponentDirichletBcSpec::getAllComponents(const Pd_shared_ptr_Array<int>& pointIds) {
	std::vector< DirichletBcSpec::ComponentLabel > c(3);
	c[0] = DirichletBcSpec::X;
	c[1] = DirichletBcSpec::Y;
	c[2] = DirichletBcSpec::Z;
	return ComponentDirichletBcSpec(c,pointIds);
}





}
