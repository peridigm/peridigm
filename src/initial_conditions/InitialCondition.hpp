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

class Epetra_Vector;

namespace InitialConditionsNS {

using Teuchos::RCP;

class InitialCondition {
public:
	InitialCondition() {}
	virtual ~InitialCondition() {}
	virtual void apply(const Epetra_Vector& x, Epetra_Vector& u, Epetra_Vector& v) = 0;
};

RCP<InitialCondition> getInstance(const Teuchos::ParameterList& peridigmParams);


}  // namespace InitialConditionsNS


#endif /* INITIALCONDITION_HPP_ */
