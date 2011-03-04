/*
 * PimpOperator.cxx
 *
 *  Created on: Feb 10, 2010
 *      Author: jamitch
 */
#include "PdImpOperator.h"
#include "PdQuickGrid.h"
#include <iostream>
#include "Pd_shared_ptr_Array.h"
#include "PdBondFilter.h"
#include "PdMaterialUtilities.h"
#include "PdVTK.h"
#include "PdZoltan.h"
#include "Epetra_Comm.h"
#include <Teuchos_RCP.hpp>

#include <map>
#include <set>

using std::cout;
using std::endl;
using Teuchos::RCP;
using namespace PdBondFilter;

namespace  PdImp{

using namespace PdMaterialUtilities;
using std::tr1::shared_ptr;
using Field_NS::Field;
using Field_NS::TemporalField;
using Field_NS::COORD3D;
using Field_NS::DISPL3D;
using Field_NS::WEIGHTED_VOLUME;
using Field_NS::VOLUME;
using Field_NS::DILATATION;

const int scalarNDF = 1;
const int vectorNDF=3;

PdImpOperator::RowOperator::RowOperator(PdImpOperator *parentOp, double h)
:
op(parentOp),
horizon(h),
/*
 * These are just initialized here with length 1; they require colMap lengths; correctly constructed below
 * after colMaps have been created
 */
xOverlapField(COORD3D,1),
uOverlapField(DISPL3D,1),
mOverlapField(WEIGHTED_VOLUME,1),
volOverlapField(VOLUME,1),
dilatationOverlapField(DILATATION,1),
rowStiffnessPtr(),
numColumnsPerRow(),
rowMatrixNeighborhoodList(0,shared_ptr<int>(),shared_ptr<int>()),
/*
 * NOTE: THESE MAPS are just place holders and initialized here;  they will be overwritten below
 */
colMapNDF(PdQuickGrid::getOverlapMap(parentOp->getEpetra_Comm(),op->rowMatrixPdGridData,vectorNDF)),
colMapScalar(PdQuickGrid::getOverlapMap(parentOp->getEpetra_Comm(),op->rowMatrixPdGridData,scalarNDF))
{

	/*
	 * Re-set neighborhood list
	 * 1) Set sizeNeighborhoodList=1
	 * 2) Create a new and empty neighborhood
	 * 3) Set neighborhood pointer for each point to 0
	 */
	op->rowMatrixPdGridData.sizeNeighborhoodList=1;
	shared_ptr<int> neighborhoodList(new int[op->rowMatrixPdGridData.sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
	op->rowMatrixPdGridData.neighborhood = neighborhoodList;
	int *neighborhood = neighborhoodList.get();
	/*
	 * number of neighbors for every point is zero
	 */
	*neighborhood = 0;

	/*
	 * 1) Create a new neighborhood pointer
	 * 2) Re-set neighborhood pointer to point to first entry in list above
	 * 3) Create new column Map
	 */
	shared_ptr<int> neighborhoodPtr(new int[op->getNumOwnedPoints()],PdQuickGrid::Deleter<int>());
	op->rowMatrixPdGridData.neighborhoodPtr = neighborhoodPtr;
	int *neighPtr = neighborhoodPtr.get();
	for(int p=0;p<op->numOwnedPoints;p++,neighPtr++)
		*neighPtr=0;

	/*
	 * Note multiplication of horizon by "2"
	 */
	RCP<BondFilter> filterPtr=RCP<BondFilter>(new BondFilterWithSelf());
	op->rowMatrixPdGridData = createAndAddNeighborhood(op->rowMatrixPdGridData,2*h,filterPtr);
	/*
	 * Need to create map before computing "local" neighborhood list since map is created from "global" ids
	 * NOTE: this is really the "overlap" map
	 */
	colMapNDF = PdQuickGrid::getOverlapMap(parentOp->getEpetra_Comm(),op->rowMatrixPdGridData,vectorNDF);
	colMapScalar = PdQuickGrid::getOverlapMap(parentOp->getEpetra_Comm(),op->rowMatrixPdGridData,scalarNDF);

	/*
	 * Initialize fields with newly created maps for RowMatrix
	 */
	xOverlapField = Field<double>(COORD3D,colMapScalar.NumMyElements());
	uOverlapField = Field<double>(DISPL3D,colMapScalar.NumMyElements());
	mOverlapField = Field<double>(WEIGHTED_VOLUME,colMapScalar.NumMyElements());
	volOverlapField = Field<double>(VOLUME,colMapScalar.NumMyElements());
	dilatationOverlapField = Field<double>(DILATATION,colMapScalar.NumMyElements());
	/*
	 * Communicate coordinates, volume, and weighted volume overlap vector
	 */
	Epetra_Import importNDF(colMapNDF,op->ownedMapNDF);
	Epetra_Vector xOverlap(View,colMapNDF,xOverlapField.getArray().get());
	Epetra_Vector xOwned(View,op->ownedMapNDF,op->xOwnedPtr.get());
	xOverlap.Import(xOwned,importNDF,Insert);
	/*
	 *Volume
	 */
	Epetra_Import importScalar(colMapScalar,op->ownedMapScalar);
	Epetra_Vector vOverlap(View,colMapScalar,volOverlapField.getArray().get());
	Epetra_Vector vOwned(View,op->ownedMapScalar,op->rowMatrixPdGridData.cellVolume.get());
	vOverlap.Import(vOwned,importScalar,Insert);
	/*
	 * Weighted volume
	 */
	Epetra_Vector mOverlap(View,colMapScalar,mOverlapField.getArray().get());
	Epetra_Vector mOwned(View,op->ownedMapScalar,op->mOwnedField.getArray().get());
	mOverlap.Import(mOwned,importScalar,Insert);

	/*
	 * Convert neighborhood list to a "local" list
	 */
	rowMatrixNeighborhoodList = PdNeighborhood::NeighborhoodList(op->getNumOwnedPoints(),op->rowMatrixPdGridData.neighborhoodPtr,PdQuickGrid::getLocalNeighborList(op->rowMatrixPdGridData,colMapScalar));
	/*
	 * Computes number of columns per row
	 */
	numColumnsPerRow = computeNumColumnsPerRow();
	/*
	 * Loop over rows and determine the maximum row size;
	 * Use maximum row size to allocate enough room for all row stiffnesses
	 */
	int max = 0;
	int *numCols = numColumnsPerRow.get();
	for(std::size_t r=0;r<numColumnsPerRow.getSize();r++, numCols++){
		if (max < *numCols) max = *numCols;
	}

	rowStiffnessPtr = Pd_shared_ptr_Array<double>(9*max);

}

/**
 * This function initializes the jacobian operator for the input displacement field
 * 1) computes overlap vector for "displacement"
 * 2) computes dilatation for input displacement and then computes the overlap vector for dilatation
 */
void PdImpOperator::RowOperator::initialize(Field_NS::Field<double> uOwnedField) {

	/*
	 * Broadcast displacements to "overlap" vector
	 */
	Epetra_Import importNDF(colMapNDF,op->ownedMapNDF);
	Epetra_Vector uOverlap(View,colMapNDF,uOverlapField.getArray().get());
	Epetra_Vector uOwned(View,op->ownedMapNDF,uOwnedField.getArray().get());
	uOverlap.Import(uOwned,importNDF,Insert);

	/*
	 * Compute dilatation and broadcast
	 */
	Epetra_Import importScalar(colMapScalar,op->ownedMapScalar);
	Field<double> ownedDilatationField = op->computeOwnedDilatation(uOwnedField);
	Epetra_Vector dilatationOverlap(View,colMapScalar,dilatationOverlapField.getArray().get());
	Epetra_Vector ownedDilatation(View,op->ownedMapScalar,ownedDilatationField.getArray().get());
	dilatationOverlap.Import(ownedDilatation,importScalar,Insert);
}

/**
 * This function returns the local column ids in a particular row of the matrix
 * @param localRowID -- this should be an "owned point"
 * @return Pd_shared_ptr_Array<int>-- array of global ids in neighborhood
 */
Pd_shared_ptr_Array<int> PdImpOperator::RowOperator::getColumnLIDs(int localRowID) const {
	const int *neigh = rowMatrixNeighborhoodList.getNeighborhood(localRowID);
	int numCols = *neigh; neigh++;
	Pd_shared_ptr_Array<int> lIds(numCols);
	int *idsPtr = lIds.get();
	for(;idsPtr != lIds.end();idsPtr++,neigh++){
		*idsPtr=*neigh;
	}
	return lIds;
}


const Pd_shared_ptr_Array<int>& PdImpOperator::RowOperator::getNumColumnsPerRow() const {
	return numColumnsPerRow;
}

/**
 * This function returns the number of points in a neighborhood -- which is also equal to number of columns in a row
 * @return Pd_shared_ptr_Array<int> -- array of number of columns associated
 * with each owned point (row)
 */
Pd_shared_ptr_Array<int> PdImpOperator::RowOperator::computeNumColumnsPerRow() const {
	int numPoints = op->getNumOwnedPoints();
	Pd_shared_ptr_Array<int> numCols(numPoints);
	int *numColsPtr = numCols.get();
	for(int I=0;I<numPoints;I++,numColsPtr++){
		*numColsPtr = rowMatrixNeighborhoodList.getNumNeigh(I);
	}
	return numCols;
}


const Pd_shared_ptr_Array<double>& PdImpOperator::RowOperator::computeRowStiffness(int localRowID, Pd_shared_ptr_Array<int> rowLIDs){

//	std::cout << "PimpOperator::RowOperator::computeRowStiffness; localRowID = " << localRowID << std::endl;

	int numColsRowI = rowLIDs.getSize();

	std::vector< std::pair<int,int> > pairs(numColsRowI);
	{
		int *cols = rowLIDs.get();
		for(std::size_t i=0;i<pairs.size();i++,cols++){
			pairs[i] = std::make_pair(*cols,i);
		}
	}
	std::map<int,int> map(pairs.begin(),pairs.end());

	const double zero = 0.0;
	double *kRowI = rowStiffnessPtr.get();
	double *k3x3Ptr = k3x3;
	double *x = xOverlapField.getArray().get();
	double *u = uOverlapField.getArray().get();
	double *volume = volOverlapField.getArray().get();
	double *m = mOverlapField.getArray().get();
	double *dilatationOverlap = dilatationOverlapField.getArray().get();
	/*
	 * This neighborhood includes 'I'
	 */
	int I = localRowID;

	/*
	 * Initialize row stiffness to zero
	 */
	PdITI::SET(kRowI,kRowI+9*numColsRowI,zero);
	/*
	 * Loop over neighbors of cell I and assemble stiffness
	 */
	for(int *colsI = rowLIDs.get();colsI!=rowLIDs.end();colsI++){
		int Q = *colsI;

		if(Q==I) continue;

		for(int *colsP = rowLIDs.get();colsP!=rowLIDs.end();colsP++){
			int P = *colsP;

			/*
			 * Initialize local 3x3 matrices to zero
			 * Assemble into 3x3
			 * Copy into row
			 */
			PdITI::SET(k3x3Ptr,k3x3Ptr+27,zero);
			double *k1 = k3x3Ptr;
			double *k2 = k3x3Ptr+9;
			double *k3 = k3x3Ptr+18;
			/*
			 * Check for existence of bond IQ
			 */
			(I==Q || I==P) ? NULL :              op->fIntPtr->kIPQ3x3(I,P,Q,x,u,volume,m,dilatationOverlap,k1,horizon);
			(P==I || P==Q) ? NULL : PdITI::SUBTRACTINTO(op->fIntPtr->kIPQ3x3(P,I,Q,x,u,volume,m,dilatationOverlap,k2,horizon),k2+9,k1);
			(Q==I || Q==P) ? NULL :      PdITI::SUMINTO(op->fIntPtr->kIPQ3x3(Q,I,P,x,u,volume,m,dilatationOverlap,k3,horizon),k3+9,k1);

			/*
			 * Assembly into Row Matrix
			 * Skip Diagonal : Will sum rows at bottom to place diagonal
			 */
			double *kIQ = kRowI + 9 * map[Q];
			PdITI::SUMINTO(k3x3Ptr,k3x3Ptr+9,kIQ);

		}

	}

	/*
	 * Final Step: Set diagonal entry by summing row and placing on diagonal
	 */
	double *kII = kRowI + 9 * map[I];
	for(int q = 0;q<numColsRowI;q++){
		if(map[I]==q) continue;
		double *kIQ = kRowI + 9 * q;
		PdITI::SUBTRACTINTO(kIQ,kIQ+9,kII);
	}

	return rowStiffnessPtr;
}



/**
 * @param u Linearization about this displacement field
 * @param horizon -- Horizon used in the usual Pd Operator;  For stiffness
 * calculations, the horizon is doubled internally
 */
shared_ptr<RowStiffnessOperator> PdImpOperator::getRowStiffnessOperator(Field_NS::Field<double> uOwned, double horizon) {

	if (rowStiffnessOperatorPtr==shared_ptr<RowStiffnessOperator>()){
		rowStiffnessOperatorPtr = shared_ptr<RowOperator>(new RowOperator(this,horizon));
	}

	/*
	 * Reset row counter on operator
	 */
	rowStiffnessOperatorPtr.get()->initialize(uOwned);

	return rowStiffnessOperatorPtr;
}

PdImpOperator::PdImpOperator(const Epetra_Comm& comm, PdGridData& gridData)
:
epetraComm(comm),
numOwnedPoints(gridData.numPoints),
ownedMapScalar(PdQuickGrid::getOwnedMap(comm,gridData,scalarNDF)),
ownedMapNDF(PdQuickGrid::getOwnedMap(comm,gridData,vectorNDF)),
overlapMapScalar(PdQuickGrid::getOverlapMap(comm,gridData,scalarNDF)),
overlapMapNDF(PdQuickGrid::getOverlapMap(comm,gridData,vectorNDF)),
importScalar(overlapMapScalar,ownedMapScalar),
importNDF(overlapMapNDF,ownedMapNDF),
exportAssembly(overlapMapNDF,ownedMapNDF),
xOwnedPtr(gridData.myX),
ownedDSF_Ptr(new double[ownedMapScalar.NumMyElements()*scalarNDF],PdQuickGrid::Deleter<double>()),
mOwnedField(Field_NS::getWEIGHTED_VOLUME(gridData.numPoints)),
dilatationOwnedField(Field_NS::Field<double>(Field_NS::DILATATION,gridData.numPoints)),
bondDamage(gridData.sizeNeighborhoodList-gridData.numPoints),
xOverlapPtr(),
uOverlapPtr(),
volumeOverlapPtr(),
yOverlapPtr(),
fInternalOverlapPtr(),
localList(gridData.numPoints,gridData.neighborhoodPtr,PdQuickGrid::getLocalNeighborList(gridData,overlapMapScalar)),
rowStiffnessOperatorPtr(),
rowMatrixPdGridData(gridData),
fIntPtr(),
temporalBondFields()
{

	/*
	 * Pull in shared coordinates
	 * NOTE View versions of Epetra_Vectors
	 */
	xOverlapPtr = allocOverlap(vectorNDF);
	Epetra_Vector xOverlap(View,overlapMapNDF,xOverlapPtr.get());
	Epetra_Vector xOwned(View,ownedMapNDF,xOwnedPtr.get());
	xOverlap.Import(xOwned,importNDF,Insert);

	/*
	 * Allocated overlap vector for displacements
	 */
	uOverlapPtr = allocOverlap(vectorNDF);

	/*
	 * Pull in shared volume
	 * NOTE View versions of Epetra_Vectors
	 */
	volumeOverlapPtr = allocOverlap(scalarNDF);
	Epetra_Vector volumeOverlap(View,overlapMapScalar,volumeOverlapPtr.get());
	Epetra_Vector volOwned(View,ownedMapScalar,gridData.cellVolume.get());
	volumeOverlap.Import(volOwned,importScalar,Insert);

	/*
	 * Compute weighted volume
	 */
	computeWeightedVolume(xOverlapPtr.get(),volumeOverlapPtr.get(),mOwnedField.getArray().get(),numOwnedPoints,localList.getNeighborhood().get());
	/*
	 * Allocate overlap pointer for "y"
	 */
	yOverlapPtr = allocOverlap(vectorNDF);

	/*
	 * Allocate overlap pointer for the internal force vector
	 */
	fInternalOverlapPtr = allocOverlap(vectorNDF);

	/*
	 * Set DSF=1
	 */
	double* dsf = ownedDSF_Ptr.get();
	double* end = dsf+gridData.numPoints;
	for(;dsf!=end;dsf++)
		*dsf=1.0;


}


Field_NS::Field<double> PdImpOperator::computeOwnedDilatation(Field_NS::Field<double> uOwnedField) {

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

	shared_ptr<int> localNeighborList = localList.getNeighborhood();
	/*
	 *
	 * COMPUTE DILATATION
	 */
	double *theta = dilatationOwnedField.getArray().get();

	computeDilatation(xOverlapPtr.get(),yOverlapPtr.get(),mOwnedField.getArray().get(),volumeOverlapPtr.get(),bondDamage.get(),theta,localNeighborList.get(),numOwnedPoints);

	return dilatationOwnedField;
}


/**
 * Assembles into the fIntOwnedField
 */
void PdImpOperator::computeInternalForce(Field_NS::Field<double> uOwnedField, Field_NS::Field<double> fIntOwnedField, bool withDilatation) {

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

	shared_ptr<int> localNeighborList = localList.getNeighborhood();
	/*
	 *
	 * COMPUTE DILATATION
	 */
	/*
	 * THIS IS KIND OF A HACK for developing the plasticity model
	 */
	double *theta = dilatationOwnedField.getArray().get();
	if(withDilatation){
		computeDilatation(xOverlapPtr.get(),yOverlapPtr.get(),mOwnedField.getArray().get(),volumeOverlapPtr.get(),bondDamage.get(),theta,localNeighborList.get(),numOwnedPoints);
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
			ownedDSF_Ptr.get(),
			fInternalOverlapPtr.get(),
			localNeighborList.get(),
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


void PdImpOperator::addConstitutiveModel(shared_ptr<PdITI::ConstitutiveModel>& model) {

	fIntPtr = model;
	vector<FieldSpec> temporalBondSpecs = fIntPtr->registerTemporalBondVariables();
	std::size_t numPoints = getNumOwnedPoints();
	std::size_t allocateBondVarSize = localList.getSizeNeighborhoodList() - numPoints;

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

void PdImpOperator::advanceStateVariables() {

	/*
	 * Loop over temporal fields and advance state
	 */
	for(std::size_t f=0;f<temporalBondFields.size();f++){
		temporalBondFields[f].advanceStep();
	}

}

std::size_t PdImpOperator::getNumPointsOverlap() const { return overlapMapScalar.NumMyElements(); }
std::size_t PdImpOperator::getNumOwnedPoints() const { return numOwnedPoints; }

shared_ptr<double> PdImpOperator::allocOwned(int ndf) const {
	return shared_ptr<double>(new double[ownedMapScalar.NumMyElements()*ndf],PdQuickGrid::Deleter<double>());
}

shared_ptr<double> PdImpOperator::allocOverlap(int ndf) const {
	return shared_ptr<double>(new double[overlapMapScalar.NumMyElements()*ndf],PdQuickGrid::Deleter<double>());
}

shared_ptr<int> PdImpOperator::allocOwned() const {
	return shared_ptr<int>(new int[ownedMapScalar.NumMyElements()],PdQuickGrid::Deleter<int>());
}

shared_ptr<int> PdImpOperator::allocOverlap() const {
	return shared_ptr<int>(new int[overlapMapScalar.NumMyElements()],PdQuickGrid::Deleter<int>());
}


}
