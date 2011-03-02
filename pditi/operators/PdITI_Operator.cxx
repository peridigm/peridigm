/*
 * PdITI_Operator.cxx
 *
 *  Created on: Mar 2, 2011
 *      Author: jamitch
 */
#include "PdITI_Operator.h"
#include "Epetra_Comm.h"

namespace PdITI {


PdITI_Operator::PdITI_Operator(const Epetra_Comm& comm, const PDNEIGH::NeighborhoodList neighborhoodList)
:
epetraComm(comm),
list(neighborhoodList)
{


}

}
