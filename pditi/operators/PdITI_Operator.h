/*
 * PdITI_Operator.h
 *
 *  Created on: Mar 2, 2011
 *      Author: jamitch
 */

#ifndef PDITI_OPERATOR_H_
#define PDITI_OPERATOR_H_

#include "../pdneigh/NeighborhoodList.h"
class Epetra_Comm;

namespace PdITI {

class PdITI_Operator {
public:
	PdITI_Operator(const Epetra_Comm& comm, const PDNEIGH::NeighborhoodList neighborhoodList);

private:
	const Epetra_Comm& epetraComm;
	PDNEIGH::NeighborhoodList list;

};

}


#endif /* PDITI_OPERATOR_H_ */
