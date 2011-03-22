/*
 * PdITI_Operator.h
 *
 *  Created on: Mar 2, 2011
 *      Author: jamitch
 */

#ifndef PDITI_OPERATOR_H_
#define PDITI_OPERATOR_H_

#include "NeighborhoodList.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Array.h"
#include "vtk/Field.h"
#include "ConstitutiveModel.h"
#include "RowStiffnessOperator.h"
#include <tr1/memory>

using std::tr1::shared_ptr;
using Field_NS::TemporalField;
using Field_NS::Field;
using UTILITIES::Array;
class Epetra_Comm;

namespace PdITI {

class PdITI_Operator {
public:
	PdITI_Operator(shared_ptr<Epetra_Comm> comm, const PDNEIGH::NeighborhoodList neighborhoodList, shared_ptr<double> ownedCellVolume);
	const PDNEIGH::NeighborhoodList get_row_matrix_neighborhood() const { return row_matrix_list; }
	void addConstitutiveModel(shared_ptr<ConstitutiveModel>& model);
	void advanceStateVariables();
	Field<double> computeOwnedDilatation(Field<double> uOwnedField);
	void computeInternalForce(Field<double> uOwned, Field<double> fIntOwned, bool withDilatation=true);
	shared_ptr<RowStiffnessOperator> getJacobian(Field<double> uOwned);
	Field<double> getWeightedVolume() const { return mOwnedField; }
	shared_ptr<Epetra_Comm> getEpetra_Comm() const { return epetraComm; }
	const Epetra_BlockMap& getOwnedMapScalar() { return  ownedMapScalar; }
	const Epetra_BlockMap& getOverlapMapScalar() { return  overlapMapScalar; }
	const Epetra_BlockMap& getOwnedMapNDF()   const { return  ownedMapNDF; }
	const Epetra_BlockMap& getOverlapMapNDF() const { return  overlapMapNDF; }

	class RowOperator;
	friend class RowOperator;
	class RowOperator : public RowStiffnessOperator {
	public:
		virtual ~RowOperator() {}
		int getNumRows() const { return row_matrix_list.get_num_owned_points(); }
		Array<int> getColumnLIDs(int localRowID) const;
		const Array<int>& getNumColumnsPerRow() const;
		const Array<double>& computeRowStiffness(int localRowID, Array<int> rowGIDs);
		const Epetra_BlockMap& getRowMap() const { return rowMapNDF; }
		const Epetra_BlockMap& getColMap() const { return colMapNDF; }
		RowOperator
		(
				shared_ptr<Epetra_Comm> comm,
				const PDNEIGH::NeighborhoodList row_matrix_list_2_horizon,
				shared_ptr<double> ownedCellVolume,
				shared_ptr<double> mOwnedPtr,
				shared_ptr<double> dsfOwnedPtr
		);

		void initialize(Field<double> uOwnedField, Field<double> ownedDilatationField, shared_ptr<ConstitutiveModel> fIntPtr);

	private:
		PDNEIGH::NeighborhoodList row_matrix_list;
		Epetra_BlockMap rowMapNDF, rowMapScalar, colMapNDF, colMapScalar;
		shared_ptr<double> xOverlapPtr, uOverlapPtr;
		shared_ptr<double> dsfOverlapPtr, mOverlapPtr, volOverlapPtr, dilatationOverlapPtr;
		Array<double>  rowStiffnessPtr;
		Array<int> numColumnsPerRow;
		double k3x3[27];
		shared_ptr<ConstitutiveModel> matPtr;
		Array<int> computeNumColumnsPerRow() const;
	};


private:
	shared_ptr<Epetra_Comm> epetraComm;
	PDNEIGH::NeighborhoodList list,row_matrix_list;
	const Epetra_BlockMap ownedMapScalar,ownedMapNDF;
	const Epetra_BlockMap overlapMapScalar,overlapMapNDF;
	Epetra_Import importScalar,importNDF;
	Epetra_Export exportAssembly;
	Field_NS::Field<double> mOwnedField, dilatationOwnedField;
	shared_ptr<double> ownedVolPtr;
	shared_ptr<double> bondDamagePtr, ownedDSF_Ptr;
	shared_ptr<double> xOverlapPtr, uOverlapPtr, yOverlapPtr, fInternalOverlapPtr, volumeOverlapPtr;
	shared_ptr<ConstitutiveModel> fIntPtr;
	vector< TemporalField<double> > temporalBondFields;
	shared_ptr<RowOperator> rowStiffnessOperatorPtr;
};

}


#endif /* PDITI_OPERATOR_H_ */
