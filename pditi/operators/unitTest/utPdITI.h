/*
 * utPdITI.h
 *
 *  Created on: Feb 22, 2011
 *      Author: jamitch
 */

#ifndef UTPDITI_H_
#define UTPDITI_H_

#include "../RowStiffnessOperator.h"
#include "../PdImpMaterials.h"
#include "../StageComponentDirichletBc.h"
#include "../DirichletBcSpec.h"
#include "vtk/Field.h"
#include <tr1/memory>
using std::tr1::shared_ptr;

class Epetra_MpiComm;
class Epetra_CrsGraph;
class Epetra_RowMatrix;


namespace utPdITI {
	shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian);
	vector<PdImp::DirichletBcSpec::ComponentLabel> getComponents(char mask);
	PdImp::IsotropicHookeSpec getMaterialSpec(double e, double nu);
	shared_ptr<Epetra_RowMatrix>getOperator
	(
	        const Field_NS::Field<char> bcMaskFieldOverlap,
	        shared_ptr<RowStiffnessOperator>& jacobian
	);

}


#endif /* UTPDITI_H_ */
