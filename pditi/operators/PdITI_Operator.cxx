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

PdITI_Operator::PdITI_Operator(const Epetra_Comm& comm, const PDNEIGH::NeighborhoodList neighborhoodList, shared_ptr<double> ownedCellVolume)
:
epetraComm(comm),
list(neighborhoodList),
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
temporalBondFields()
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
	computeWeightedVolume(xOverlapPtr.get(),volumeOverlapPtr.get(),mOwnedField.getArray().get(),neighborhoodList.get_num_owned_points(),list.get_neighborhood().get());

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

	shared_ptr<int> localNeighborList = list.get_neighborhood();
	int numOwnedPoints = list.get_num_owned_points();
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

	shared_ptr<int> localNeighborList = list.get_neighborhood();
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


}
