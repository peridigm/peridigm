/*
 * utPdITI.h
 *
 *  Created on: Feb 22, 2011
 *      Author: jamitch
 */

#ifndef UTPDITI_H_
#define UTPDITI_H_

#include "../PdImpOperator.h"
#include "../RowStiffnessOperator.h"
#include "../PdImpOperator.h"
#include "../PdImpMaterials.h"
#include "../StageComponentDirichletBc.h"
#include <tr1/memory>
using std::tr1::shared_ptr;

class Epetra_MpiComm;
class Epetra_CrsGraph;
class Epetra_RowMatrix;


namespace utPdITI {
	shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian);
	shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm);
	PdImp::IsotropicHookeSpec getMaterialSpec(double e, double nu);
	shared_ptr<Epetra_RowMatrix>getOperator
	(
			const vector<shared_ptr<PdImp::StageComponentDirichletBc> >& bcArray,
			shared_ptr<Epetra_CrsGraph>& graphPtr,
			shared_ptr<RowStiffnessOperator>& jacobian
	);

}


#endif /* UTPDITI_H_ */
