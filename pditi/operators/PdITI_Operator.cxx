/*
 * PdITI_Operator.cxx
 *
 *  Created on: Mar 2, 2011
 *      Author: jamitch
 */
#include "PdITI_Operator.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "PdMaterialUtilities.h"
#include <iostream>
using std::cout;
using std::endl;
using namespace PdMaterialUtilities;

namespace PdITI {

const int scalarNDF = 1;
const int vectorNDF=3;

template<class T> struct ArrayDeleter{
	void operator()(T* d) {
		delete [] d;
	}
};

//PdITI_Operator::
//RowOperator::RowOperator
//(
//		const Epetra_Comm& comm,
//		const PDNEIGH::NeighborhoodList row_matrix_list_2_horizon,
//		shared_ptr<double> ownedCellVolume
//):
//row_matrix_list(row_matrix_list_2_horizon),
//rowMapNDF(row_matrix_list.getOwnedMap(comm,vectorNDF)),
//rowMapScalar(row_matrix_list.getOwnedMap(comm,scalarNDF)),
//colMapNDF(row_matrix_list.getOverlapMap(comm,vectorNDF)),
//colMapScalar(row_matrix_list.getOverlapMap(comm,scalarNDF)),
//xOverlapPtr(new double[colMapNDF.NumMyElements()*vectorNDF],ArrayDeleter<double>()),
//uOverlapPtr(new double[colMapNDF.NumMyElements()*vectorNDF],ArrayDeleter<double>()),
//mOverlapPtr(new double[colMapNDF.NumMyElements()*scalarNDF],ArrayDeleter<double>()),
//volOverlapPtr(new double[colMapNDF.NumMyElements()*scalarNDF],ArrayDeleter<double>()),
//dilatationOverlapPtr(new double[colMapNDF.NumMyElements()*scalarNDF],ArrayDeleter<double>())
//{
//
//	/*
//	 * Communicate coordinates, volume, and weighted volume overlap vector
//	 */
//	Epetra_Import importNDF(colMapNDF,rowMap);
//	Epetra_Vector xOverlap(View,colMapNDF,xOverlapPtr.get());
//	Epetra_Vector xOwned(View,rowMap,row_matrix_list.get_owned_x().get());
//	xOverlap.Import(xOwned,importNDF,Insert);
////	/*
////	 *Volume
////	 */
////	Epetra_Import importScalar(colMapScalar,rowMapScalar);
////	Epetra_Vector vOverlap(View,colMapScalar,volOverlapPtr.get());
////	Epetra_Vector vOwned(View,rowMapScalar,ownedCellVolume.get());
////	vOverlap.Import(vOwned,importScalar,Insert);
////	/*
////	 * Weighted volume
////	 */
////	Epetra_Vector mOverlap(View,colMapScalar,mOverlapPtr.get());
////	Epetra_Vector mOwned(View,rowMapScalar,op->mOwnedPtr.getArray().get());
////	mOverlap.Import(mOwned,importScalar,Insert);
//
//
//}
//
///**
// * This function initializes the jacobian operator for the input displacement field
// * 1) computes overlap vector for "displacement"
// * 2) computes dilatation for input displacement and then computes the overlap vector for dilatation
// */
//void PdITI_Operator::RowOperator::initialize(Field_NS::Field<double> uOwnedField) {
//
//	/*
//	 * Broadcast displacements to "overlap" vector
//	 */
//	Epetra_Import importNDF(colMapNDF,op->ownedMapNDF);
//	Epetra_Vector uOverlap(View,colMapNDF,uOverlapField.getArray().get());
//	Epetra_Vector uOwned(View,op->ownedMapNDF,uOwnedField.getArray().get());
//	uOverlap.Import(uOwned,importNDF,Insert);
//
//	/*
//	 * Compute dilatation and broadcast
//	 */
//	Epetra_Import importScalar(colMapScalar,op->ownedMapScalar);
//	Field<double> ownedDilatationField = op->computeOwnedDilatation(uOwnedField);
//	Epetra_Vector dilatationOverlap(View,colMapScalar,dilatationOverlapField.getArray().get());
//	Epetra_Vector ownedDilatation(View,op->ownedMapScalar,ownedDilatationField.getArray().get());
//	dilatationOverlap.Import(ownedDilatation,importScalar,Insert);
//}
//
///**
// * This function returns the local column ids in a particular row of the matrix
// * @param localRowID -- this should be an "owned point"
// * @return Pd_shared_ptr_Array<int>-- array of global ids in neighborhood
// */
//Pd_shared_ptr_Array<int> PdITI_Operator::RowOperator::getColumnLIDs(int localRowID) const {
//	const int *neigh = row_matrix_list.getgetNeighborhood(localRowID);
//	int numCols = *neigh; neigh++;
//	Pd_shared_ptr_Array<int> lIds(numCols);
//	int *idsPtr = lIds.get();
//	for(;idsPtr != lIds.end();idsPtr++,neigh++){
//		*idsPtr=*neigh;
//	}
//	return lIds;
//}
//
//
//const Pd_shared_ptr_Array<int>& PdITI_Operator::RowOperator::getNumColumnsPerRow() const {
//	return numColumnsPerRow;
//}
//
///**
// * This function returns the number of points in a neighborhood -- which is also equal to number of columns in a row
// * @return Pd_shared_ptr_Array<int> -- array of number of columns associated
// * with each owned point (row)
// */
//Pd_shared_ptr_Array<int> PdITI_Operator::RowOperator::computeNumColumnsPerRow() const {
//	int numPoints = op->getNumOwnedPoints();
//	Pd_shared_ptr_Array<int> numCols(numPoints);
//	int *numColsPtr = numCols.get();
//	for(int I=0;I<numPoints;I++,numColsPtr++){
//		*numColsPtr = rowMatrixNeighborhoodList.getNumNeigh(I);
//	}
//	return numCols;
//}
//
//
//const Pd_shared_ptr_Array<double>& PdITI_Operator::RowOperator::computeRowStiffness(int localRowID, Pd_shared_ptr_Array<int> rowLIDs){
//
////	std::cout << "PimpOperator::RowOperator::computeRowStiffness; localRowID = " << localRowID << std::endl;
//
//	int numColsRowI = rowLIDs.getSize();
//
//	std::vector< std::pair<int,int> > pairs(numColsRowI);
//	{
//		int *cols = rowLIDs.get();
//		for(std::size_t i=0;i<pairs.size();i++,cols++){
//			pairs[i] = std::make_pair(*cols,i);
//		}
//	}
//	std::map<int,int> map(pairs.begin(),pairs.end());
//
//	const double zero = 0.0;
//	double *kRowI = rowStiffnessPtr.get();
//	double *k3x3Ptr = k3x3;
//	double *x = xOverlapField.getArray().get();
//	double *u = uOverlapField.getArray().get();
//	double *volume = volOverlapField.getArray().get();
//	double *m = mOverlapField.getArray().get();
//	double *dilatationOverlap = dilatationOverlapField.getArray().get();
//	/*
//	 * This neighborhood includes 'I'
//	 */
//	int I = localRowID;
//
//	/*
//	 * Initialize row stiffness to zero
//	 */
//	PdITI::SET(kRowI,kRowI+9*numColsRowI,zero);
//	/*
//	 * Loop over neighbors of cell I and assemble stiffness
//	 */
//	for(int *colsI = rowLIDs.get();colsI!=rowLIDs.end();colsI++){
//		int Q = *colsI;
//
//		if(Q==I) continue;
//
//		for(int *colsP = rowLIDs.get();colsP!=rowLIDs.end();colsP++){
//			int P = *colsP;
//
//			/*
//			 * Initialize local 3x3 matrices to zero
//			 * Assemble into 3x3
//			 * Copy into row
//			 */
//			PdITI::SET(k3x3Ptr,k3x3Ptr+27,zero);
//			double *k1 = k3x3Ptr;
//			double *k2 = k3x3Ptr+9;
//			double *k3 = k3x3Ptr+18;
//			/*
//			 * Check for existence of bond IQ
//			 */
//			(I==Q || I==P) ? NULL :              op->fIntPtr->kIPQ3x3(I,P,Q,x,u,volume,m,dilatationOverlap,k1,horizon);
//			(P==I || P==Q) ? NULL : PdITI::SUBTRACTINTO(op->fIntPtr->kIPQ3x3(P,I,Q,x,u,volume,m,dilatationOverlap,k2,horizon),k2+9,k1);
//			(Q==I || Q==P) ? NULL :      PdITI::SUMINTO(op->fIntPtr->kIPQ3x3(Q,I,P,x,u,volume,m,dilatationOverlap,k3,horizon),k3+9,k1);
//
//			/*
//			 * Assembly into Row Matrix
//			 * Skip Diagonal : Will sum rows at bottom to place diagonal
//			 */
//			double *kIQ = kRowI + 9 * map[Q];
//			PdITI::SUMINTO(k3x3Ptr,k3x3Ptr+9,kIQ);
//
//		}
//
//	}
//
//	/*
//	 * Final Step: Set diagonal entry by summing row and placing on diagonal
//	 */
//	double *kII = kRowI + 9 * map[I];
//	for(int q = 0;q<numColsRowI;q++){
//		if(map[I]==q) continue;
//		double *kIQ = kRowI + 9 * q;
//		PdITI::SUBTRACTINTO(kIQ,kIQ+9,kII);
//	}
//
//	return rowStiffnessPtr;
//}


PdITI_Operator::PdITI_Operator(const Epetra_Comm& comm, const PDNEIGH::NeighborhoodList neighborhoodList, shared_ptr<double> ownedCellVolume)
:
epetraComm(comm),
list(neighborhoodList),
row_matrix_list(list.cloneAndShare(2.0*neighborhoodList.get_horizon())),
ownedMapScalar(neighborhoodList.getOwnedMap(comm,scalarNDF)),
ownedMapNDF(neighborhoodList.getOwnedMap(comm,vectorNDF)),
overlapMapScalar(neighborhoodList.getOverlapMap(comm,scalarNDF)),
overlapMapNDF(neighborhoodList.getOverlapMap(comm,vectorNDF)),
importScalar(overlapMapScalar,ownedMapScalar),
importNDF(overlapMapNDF,ownedMapNDF),
exportAssembly(overlapMapNDF,ownedMapNDF),
mOwnedField(Field_NS::getWEIGHTED_VOLUME(neighborhoodList.get_num_owned_points())),
dilatationOwnedField(Field_NS::Field<double>(Field_NS::DILATATION,neighborhoodList.get_num_owned_points())),
bondDamage(new double[neighborhoodList.get_size_neighborhood_list()-neighborhoodList.get_num_owned_points()],ArrayDeleter<double>()),
xOverlapPtr(new double[overlapMapScalar.NumMyElements()*vectorNDF],ArrayDeleter<double>()),
uOverlapPtr(new double[overlapMapScalar.NumMyElements()*vectorNDF],ArrayDeleter<double>()),
yOverlapPtr(new double[overlapMapScalar.NumMyElements()*vectorNDF],ArrayDeleter<double>()),
fInternalOverlapPtr(new double[overlapMapScalar.NumMyElements()*vectorNDF],ArrayDeleter<double>()),
volumeOverlapPtr(new double[overlapMapScalar.NumMyElements()*scalarNDF],ArrayDeleter<double>()),
fIntPtr(),
temporalBondFields(),
rowStiffnessOperatorPtr()
{


	/*
	 * Pull in shared coordinates
	 * NOTE View versions of Epetra_Vectors
	 */
	Epetra_Vector xOverlap(View,overlapMapNDF,xOverlapPtr.get());
	Epetra_Vector xOwned(View,ownedMapNDF,list.get_owned_x().get());
	xOverlap.Import(xOwned,importNDF,Insert);

	/*
	 * Pull in shared volume
	 * NOTE View versions of Epetra_Vectors
	 */
	Epetra_Vector volumeOverlap(View,overlapMapScalar,volumeOverlapPtr.get());
	Epetra_Vector volOwned(View,ownedMapScalar,ownedCellVolume.get());
	volumeOverlap.Import(volOwned,importScalar,Insert);

	/*
	 * Compute weighted volume
	 */
	computeWeightedVolume(xOverlapPtr.get(),volumeOverlapPtr.get(),mOwnedField.getArray().get(),list.get_num_owned_points(),list.get_local_neighborhood().get());

}

Field_NS::Field<double> PdITI_Operator::computeOwnedDilatation(Field_NS::Field<double> uOwnedField) {

	/*
	 * UPDATE GEOMETRY FOR GIVEN DISPLACMENT FIELD
	 *
	 *
	 * Broadcast  displacement to sharing processors
	 * Note VIEW Variety of Epetra_Vectors
	 */
	Epetra_Vector uOverlap(View,overlapMapNDF,uOverlapPtr.get());
	Epetra_Vector uOwned(View,ownedMapNDF,uOwnedField.getArray().get());
	uOverlap.Import(uOwned,importNDF,Insert);

	double *x = xOverlapPtr.get();
	double *u = uOverlapPtr.get();
	double *y = yOverlapPtr.get();
	double *end = x+overlapMapNDF.NumMyPoints();
	for(;x!=end;x++,u++,y++){
		*y = *x + *u;
	}

	int numOwnedPoints = list.get_num_owned_points();
	/*
	 *
	 * COMPUTE DILATATION
	 */
	double *theta = dilatationOwnedField.getArray().get();

	computeDilatation(xOverlapPtr.get(),yOverlapPtr.get(),mOwnedField.getArray().get(),volumeOverlapPtr.get(),bondDamage.get(),theta,list.get_local_neighborhood().get(),numOwnedPoints);

	return dilatationOwnedField;
}

/**
 * Assembles into the fIntOwnedField
 */
void PdITI_Operator::computeInternalForce(Field_NS::Field<double> uOwnedField, Field_NS::Field<double> fIntOwnedField, bool withDilatation) {
	/*
	 * UPDATE GEOMETRY FOR GIVEN DISPLACMENT FIELD
	 *
	 *
	 * Broadcast  displacement to sharing processors
	 * Note VIEW Variety of Epetra_Vectors
	 */
	Epetra_Vector uOverlap(View,overlapMapNDF,uOverlapPtr.get());
	Epetra_Vector uOwned(View,ownedMapNDF,uOwnedField.getArray().get());
	uOverlap.Import(uOwned,importNDF,Insert);

	double *x = xOverlapPtr.get();
	double *u = uOverlapPtr.get();
	double *y = yOverlapPtr.get();
	double *f = fInternalOverlapPtr.get();
	double *end = x+overlapMapNDF.NumMyPoints();
	for(;x!=end;x++,u++,y++,f++){
		*y = *x + *u;
		*f = 0;
	}

	int numOwnedPoints = list.get_num_owned_points();
	/*
	 *
	 * COMPUTE DILATATION
	 */
	/*
	 * THIS IS KIND OF A HACK for developing the plasticity model
	 */
	double *theta = dilatationOwnedField.getArray().get();
	if(withDilatation){
		computeDilatation(xOverlapPtr.get(),yOverlapPtr.get(),mOwnedField.getArray().get(),volumeOverlapPtr.get(),bondDamage.get(),theta,list.get_local_neighborhood().get(),numOwnedPoints);
	} else {
		dilatationOwnedField.setValue(0.0);
	}

	/*
	 * Create array of pointers to temporal bond variables
	 */
//	vector<double*> temporalBondVarsN   = vector<double*>(temporalBondFields.size());
//	vector<double*> temporalBondVarsNp1 = vector<double*>(temporalBondFields.size());
//	for(std::size_t v=0;v<temporalBondVars.size();v++){
//		Field<double> fieldN = temporalBondFields[v].getField(Field_NS::FieldSpec::STEP_N);
//		temporalBondVarsN[v] = fieldN.getArray().get();
//		Field<double> fieldNP1 = temporalBondFields[v].getField(Field_NS::FieldSpec::STEP_NP1);
//		temporalBondVarsNp1[v] = fieldNP1.getArray().get();
//	}


	/*
	 * COMPUTE INTERNAL FORCE
	 */

	fIntPtr->computeInternalForce
	(
			xOverlapPtr.get(),
			yOverlapPtr.get(),
			mOwnedField.getArray().get(),
			volumeOverlapPtr.get(),
			theta,
			bondDamage.get(),
			fInternalOverlapPtr.get(),
			list.get_local_neighborhood().get(),
			numOwnedPoints,
			temporalBondFields
	);

	/*
	 * ASSEMBLE COMPUTED FORCES across all processors
	 */
	Epetra_Vector fOverlap(View,overlapMapNDF,fInternalOverlapPtr.get());
	Epetra_Vector fOwned(View,ownedMapNDF,fIntOwnedField.getArray().get());

	/*
	 * Export "fOverlap" to fOwned
	 */
	fOwned.Export(fOverlap,exportAssembly,Add);

}

void PdITI_Operator::addConstitutiveModel(shared_ptr<PdITI::ConstitutiveModel>& model) {

	fIntPtr = model;
	vector<FieldSpec> temporalBondSpecs = fIntPtr->registerTemporalBondVariables();
	std::size_t numPoints = list.get_num_owned_points();
	std::size_t allocateBondVarSize = list.get_size_neighborhood_list() - numPoints;

	/*
	 * Create and allocate temporal bond fields
	 */
	temporalBondFields = vector<TemporalField<double> >(temporalBondSpecs.size());
	for(std::size_t f=0;f<temporalBondSpecs.size();f++){
		temporalBondFields[f] = TemporalField<double>(temporalBondSpecs[f],allocateBondVarSize);
		Field<double> fN = temporalBondFields[f].getField(FieldSpec::STEP_N);
		Field<double> fNp1 = temporalBondFields[f].getField(FieldSpec::STEP_NP1);
		fN.setValue(0.0);
		fNp1.setValue(0.0);

	}

}

void PdITI_Operator::advanceStateVariables() {

	/*
	 * Loop over temporal fields and advance state
	 */
	for(std::size_t f=0;f<temporalBondFields.size();f++){
		temporalBondFields[f].advanceStep();
	}

}


/**
 * @param u Linearization about this displacement field
 */
shared_ptr<RowStiffnessOperator> PdITI_Operator::getJacobian(Field_NS::Field<double> uOwned) {
	rowStiffnessOperatorPtr = shared_ptr<RowStiffnessOperator>();
//	if (rowStiffnessOperatorPtr==shared_ptr<RowStiffnessOperator>()){
//		rowStiffnessOperatorPtr = shared_ptr<RowOperator>(new RowOperator(this,horizon));
//	}

	/*
	 * Reset row counter on operator
	 */
//	rowStiffnessOperatorPtr.get()->initialize(uOwned);

	return rowStiffnessOperatorPtr;
}


}
