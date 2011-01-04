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

class InitialCondition {
	virtual ~InitialCondition() {}
public:
	virtual void apply() = 0;
};

}  // namespace InitialConditionsNS

} // namespace PeridigmNS

#endif /* INITIALCONDITION_HPP_ */
