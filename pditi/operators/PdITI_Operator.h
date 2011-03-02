/*
 * PdITI_Operator.h
 *
 *  Created on: Mar 2, 2011
 *      Author: jamitch
 */

#ifndef PDITI_OPERATOR_H_
#define PDITI_OPERATOR_H_

#include "../pdneigh/NeighborhoodList.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Field.h"
#include "ConstitutiveModel.h"
#include <tr1/memory>

using std::tr1::shared_ptr;
using Field_NS::TemporalField;
class Epetra_Comm;

namespace PdITI {

class PdITI_Operator {
public:
	PdITI_Operator(const Epetra_Comm& comm, const PDNEIGH::NeighborhoodList neighborhoodList);

private:
	const Epetra_Comm& epetraComm;
	PDNEIGH::NeighborhoodList list;
	const Epetra_BlockMap ownedMapScalar,ownedMapNDF;
	const Epetra_BlockMap overlapMapScalar,overlapMapNDF;
	Epetra_Import importScalar,importNDF;
	Epetra_Export exportAssembly;
	Field_NS::Field<double> mOwnedField, dilatationOwnedField;
	shared_ptr<double> bondDamage;
	shared_ptr<double> xOverlapPtr, uOverlapPtr, yOverlapPtr, fInternalOverlapPtr, volumeOverlapPtr;
//	shared_ptr<ConstitutiveModel> fIntPtr;
	vector< TemporalField<double> > temporalBondFields;

};

}


#endif /* PDITI_OPERATOR_H_ */
