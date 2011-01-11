/*
 * InitialCondition.cpp
 *
 *  Created on: Jan 10, 2011
 *      Author: jamitch
 */

#include "InitialCondition.hpp"
#include "Q2Cylinder.hpp"
#include <map>
#include <string>



namespace PeridigmNS {

namespace InitialConditionsNS {

using std::map;
using std::string;

/*
 * This is a function pointer to a generic builder
 */
typedef RCP<InitialCondition> (*BuilderFunctionPointer)(const Teuchos::RCP<Teuchos::ParameterList>&);



/*
 * Forward declarations of function pointers to builders
 */
RCP<InitialCondition> getQ2CylinderIC(const Teuchos::RCP<Teuchos::ParameterList>& params);

/*
 * Creates map of support functions
 */
map<string,BuilderFunctionPointer> getSupport() {
	map<string,BuilderFunctionPointer> icSupport;
	icSupport["Q2Cylinder"] = getQ2CylinderIC;
	return icSupport;
}

RCP<InitialCondition> getInstance(const Teuchos::RCP<Teuchos::ParameterList>& peridigmParams) {
	map<string,BuilderFunctionPointer> icSupport = getSupport();

	Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(&(peridigmParams->sublist("Problem").sublist("Initial Conditions")), false);
	string icType = params->get<string>("Type");
	TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "Q2Cylinder", "PeridigmNS::InitialConditionsNS::Q2Cylinder -- invalid \'Type\'");

	if (!params->isSublist("IC Params")){
		TEST_FOR_EXCEPT_MSG(true, "PeridigmNS::InitialConditionsNS::getInstance-->missing \'IC Params\' sublist.");
	}
	Teuchos::RCP<Teuchos::ParameterList> ICParams = Teuchos::rcp(&(params->sublist("IC Params")), false);
	BuilderFunctionPointer builder = icSupport[icType];
	return builder(ICParams);
}


RCP<InitialCondition> getQ2CylinderIC(const Teuchos::RCP<Teuchos::ParameterList>& params) {
	return RCP<InitialCondition>(0);
}

}  // namespace InitialConditionsNS

} // namespace PeridigmNS
