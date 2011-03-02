/*
 * PdITI_Operator.cxx
 *
 *  Created on: Mar 2, 2011
 *      Author: jamitch
 */
#include "PdITI_Operator.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"


namespace PdITI {

const int scalarNDF = 1;
const int vectorNDF=3;

template<class T> struct ArrayDeleter{
	void operator()(T* d) {
		delete [] d;
	}
};

PdITI_Operator::PdITI_Operator(const Epetra_Comm& comm, const PDNEIGH::NeighborhoodList neighborhoodList)
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
temporalBondFields()
{


}

}
