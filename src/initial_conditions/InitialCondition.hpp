/*
 * InitialCondition.hpp
 *
 *  Created on: Jan 4, 2011
 *      Author: jamitch
 */

#ifndef INITIALCONDITION_HPP_
#define INITIALCONDITION_HPP_

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace PeridigmNS {

namespace InitialConditionsNS {

using Teuchos::RCP;





class InitialCondition {
public:
	virtual ~InitialCondition() {}
	virtual void apply() = 0;
};


RCP<InitialCondition> getInstance(const Teuchos::RCP<Teuchos::ParameterList>& params);


}  // namespace InitialConditionsNS

} // namespace PeridigmNS

#endif /* INITIALCONDITION_HPP_ */
