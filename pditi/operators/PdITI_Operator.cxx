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

}
