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
	PdITI_Operator(const Epetra_Comm& comm, const PDNEIGH::NeighborhoodList neighborhoodList, shared_ptr<double> ownedCellVolume);
	void addConstitutiveModel(shared_ptr<ConstitutiveModel>& model);
	void advanceStateVariables();
	Field_NS::Field<double> computeOwnedDilatation(Field_NS::Field<double> uOwnedField);
	void computeInternalForce(Field_NS::Field<double> uOwned, Field_NS::Field<double> fIntOwned, bool withDilatation=true);
	Field_NS::Field<double> getWeightedVolume() const { return mOwnedField; }
	const Epetra_Comm& getEpetra_Comm() const { return epetraComm; }
	const Epetra_BlockMap& getOwnedMapScalar() { return  ownedMapScalar; }
	const Epetra_BlockMap& getOverlapMapScalar() { return  overlapMapScalar; }
	const Epetra_BlockMap& getOwnedMapNDF()   const { return  ownedMapNDF; }
	const Epetra_BlockMap& getOverlapMapNDF() const { return  overlapMapNDF; }

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
	shared_ptr<ConstitutiveModel> fIntPtr;
	vector< TemporalField<double> > temporalBondFields;

};

}


#endif /* PDITI_OPERATOR_H_ */
