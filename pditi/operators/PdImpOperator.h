/*
 * PdOperator.cxx
 *
 *  Created on: Feb 10, 2010
 *      Author: jamitch
 */

#ifndef PDIMP_OPERATOR_H_
#define PDIMP_OPERATOR_H_

#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_BlockMap.h"
#include "PdGridData.h"
#include "PdImpMaterials.h"
#include "Pd_shared_ptr_Array.h"
#include "RowStiffnessOperator.h"
#include "PdNeighborhood.h"
#include "PdITI_Utilities.h"
#include "ConstitutiveModel.h"
#include "Field.h"
#include <tr1/memory>
#include <vector>

using std::tr1:: shared_ptr;
using std::vector;
using Field_NS::Field;
using Field_NS::TemporalField;
class Epetra_BlockMap;
class Epetra_Comm;

namespace PdImp {

class PdImpOperator {
public:
	PdImpOperator(const Epetra_Comm& comm, PdGridData& gridData);
	void addConstitutiveModel(shared_ptr<PdITI::ConstitutiveModel>& model);
	void advanceStateVariables();
	Field_NS::Field<double> getWeightedVolume() const { return mOwnedField; }
	const Epetra_Comm& getEpetra_Comm() const { return epetraComm; }
	const Epetra_BlockMap& getOwnedMapScalar() { return  ownedMapScalar; }
	const Epetra_BlockMap& getOverlapMapScalar() { return  overlapMapScalar; }
	const Epetra_BlockMap& getOwnedMapNDF()   const { return  ownedMapNDF; }
	const Epetra_BlockMap& getOverlapMapNDF() const { return  overlapMapNDF; }
	Field_NS::Field<double> computeOwnedDilatation(Field_NS::Field<double> uOwned);
	void computeInternalForce(Field_NS::Field<double> uOwned, Field_NS::Field<double> fIntOwned, bool withDilatation=true);
	std::size_t getNumPointsOverlap() const;
	std::size_t getNumOwnedPoints() const;
	std::tr1::shared_ptr<double> allocOwned(int ndf) const;
	std::tr1::shared_ptr<double> allocOverlap(int ndf) const;
	std::tr1::shared_ptr<int> allocOwned() const;
	std::tr1::shared_ptr<int> allocOverlap() const;
	std::tr1::shared_ptr<RowStiffnessOperator> getRowStiffnessOperator(Field_NS::Field<double> uOwned, double horizon);

	class RowOperator;
	friend class RowOperator;
	class RowOperator : public RowStiffnessOperator {
	public:
		~RowOperator() {}
		int getNumRows() const { return op->getNumOwnedPoints(); }
		Pd_shared_ptr_Array<int> getColumnLIDs(int localRowID) const;
		const Pd_shared_ptr_Array<int>& getNumColumnsPerRow() const;
		const Pd_shared_ptr_Array<double>& computeRowStiffness(int localRowID, Pd_shared_ptr_Array<int> rowGIDs);
		const Epetra_BlockMap& getRowMap() const { return op->getOwnedMapNDF(); }
		const Epetra_BlockMap& getColMap() const { return colMapNDF; }
		RowOperator(PdImpOperator* parentOperator, double horizon);
		void initialize(Field_NS::Field<double> uOwnedField);

	private:
		PdImpOperator *op;
		double horizon;
		Field_NS::Field<double>  xOverlapField, uOverlapField;
		Field_NS::Field<double> mOverlapField, volOverlapField, dilatationOverlapField;
		Pd_shared_ptr_Array<double>  rowStiffnessPtr;
		Pd_shared_ptr_Array<int> numColumnsPerRow;
		PdNeighborhood::NeighborhoodList rowMatrixNeighborhoodList;
		Epetra_BlockMap colMapNDF, colMapScalar;
		double k3x3[27];
		Pd_shared_ptr_Array<int> computeNumColumnsPerRow() const;
	};

private:
	const Epetra_Comm& epetraComm;
	int numOwnedPoints;
	const Epetra_BlockMap ownedMapScalar,ownedMapNDF;
	const Epetra_BlockMap overlapMapScalar,overlapMapNDF;
	Epetra_Import importScalar,importNDF;
	Epetra_Export exportAssembly;
	std::tr1::shared_ptr<double> xOwnedPtr;
	std::tr1::shared_ptr<double> ownedDSF_Ptr;
	Field_NS::Field<double> mOwnedField, dilatationOwnedField;
	Pd_shared_ptr_Array<double> bondDamage;
	std::tr1::shared_ptr<double> xOverlapPtr, uOverlapPtr, volumeOverlapPtr;
	std::tr1::shared_ptr<double> yOverlapPtr, fInternalOverlapPtr;
	PdNeighborhood::NeighborhoodList localList;
	std::tr1::shared_ptr<RowOperator> rowStiffnessOperatorPtr;
	PdGridData rowMatrixPdGridData;
	shared_ptr<PdITI::ConstitutiveModel> fIntPtr;
	vector< TemporalField<double> > temporalBondFields;

};

}
#endif //PDIMP_OPERATOR_H_
