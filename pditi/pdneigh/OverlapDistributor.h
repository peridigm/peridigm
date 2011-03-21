/*
 * Overlap_Distributor.h
 *
 *  Created on: Mar 18, 2011
 *      Author: wow
 */

#ifndef OVERLAP_DISTRIBUTOR_H_
#define OVERLAP_DISTRIBUTOR_H_

#include "vtk/Field.h"
#include "Array.h"
#include "NeighborhoodList.h"

#include "Epetra_Distributor.h"


#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


#include <tr1/memory>
#include <iostream>

namespace PDNEIGH {

using Field_NS::Field;
using Field_NS::FieldSpec;
using UTILITIES::Array;
using std::tr1::shared_ptr;
using std::cout;
using std::endl;

template<typename T>
inline Field<T> createOverlapField(NeighborhoodList& list, Field<T> ownedField) {

	/*
	 * Initialize distributor object
	 */
	shared_ptr<Epetra_Comm> comm = list.get_Epetra_Comm();
	int myRank = comm->MyPID();
	int numProcs = comm->NumProc();
	shared_ptr<Epetra_Distributor> distributor(comm->CreateDistributor());

	shared_ptr<int> sharedGIDs = list.get_shared_gids();
	size_t num_import = list.get_num_shared_points();
	size_t num_owned_points = ownedField.get_num_points();

	shared_ptr<Epetra_BlockMap> ownedMap = list.getOwnedMap(1);
	Array<int> PIDs(num_import), LIDs(num_import);
	ownedMap->RemoteIDList(num_import,sharedGIDs.get(),PIDs.get(),LIDs.get());

	int num_export, *exportGIDs, *exportPIDs;
	distributor->CreateFromRecvs(num_import, sharedGIDs.get(), PIDs.get(),true, num_export, exportGIDs, exportPIDs);

	/*
	 * Gather data for export
	 */
	int *eGIDs=exportGIDs;
	Array<T> exportData(num_export);
	for(size_t i=0;i<num_export;i++,eGIDs++){
		int lid = ownedMap->LID(*eGIDs);
		exportData[i] = ownedField[lid];
	}

	/*
	 * Apply communication plan
	 */
	int size_import=0;
	char *data = 0;
	distributor->Do((char*)exportData.get(),sizeof(T),size_import,data);
	T* importData = reinterpret_cast<T*>(data);

	/*
	 * Gather data received and stuff into 'overlapField'
	 */
	const FieldSpec spec = ownedField.get_spec();
	Field<T> overlapField (spec,num_import+num_owned_points);
	shared_ptr<Epetra_BlockMap> overlapMap = list.getOverlapMap(spec.getLength());

	/*
	 * Copy owned data over
	 */
	for(size_t n=0;n<num_owned_points;n++)
		overlapField[n] = ownedField[n];

	/*
	 * Now import shared id data
	 */
	int *sGIDs = sharedGIDs.get();
	T* iData = importData;
	for(size_t i=0;i<num_import;i++,sGIDs++){
		int lid = overlapMap->LID(*sGIDs);
		overlapField[lid]=iData[i];
	}

	delete [] exportGIDs;
	delete [] exportPIDs;
	delete [] data;
	return overlapField;
}

}

#endif /* OVERLAP_DISTRIBUTOR_H_ */
