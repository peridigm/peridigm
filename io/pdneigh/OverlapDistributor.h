/*
 * Overlap_Distributor.h
 *
 *  Created on: Mar 18, 2011
 *      Author: jamitch
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

namespace PDNEIGH {

using Field_NS::Field;
using Field_NS::FieldSpec;
using UTILITIES::Array;
using std::tr1::shared_ptr;

template<typename T>
inline Field<T> createOverlapField(NeighborhoodList& list, Field<T> ownedField) {

    /*
     * Initialize distributor object
     */
    shared_ptr<const Epetra_Comm> comm = list.get_Epetra_Comm();
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

    /*
     * Data associated with each GID sent will consist of:
     * 1) GID
     * 2) data
     * This allows the receiving processor to associate the data with the proper GID
     */
    int object_size = sizeof(int) + ownedField.getLength() * sizeof(T);
    int *eGIDs=exportGIDs;
    Array<char> exportData(num_export*object_size);
    char *exportDataBuf = exportData.get();
    T* ownedData = ownedField.get();
    for(size_t i=0;i<num_export;i++,eGIDs++){

        /*
         * Use this temporary pointer for copying data of single point
         */
        char *dest = exportDataBuf + i * object_size;

        /*
         * Copy GID into buffer
         */
        int numBytes = sizeof(int);
        memcpy((void*)dest,(void*)eGIDs,numBytes);

        /*
         * Advance pointer
         */
        dest += numBytes;

        /*
         * Copy from owned field into buffer
         */
        int lid = ownedMap->LID(*eGIDs);
        numBytes = sizeof(T) * ownedField.getLength();
        T* o = ownedData + lid * ownedField.getLength();
        memcpy((void*)dest,(void*)o,numBytes);
    }

    /*
     * Apply communication plan
     */
    int size_import=0;
    char *importData = 0;
    distributor->Do(exportData.get(),object_size,size_import,importData);

    /*
     * Gather data received and stuff into 'overlapField'
     */
    Field<T> overlapField (ownedField.get_overlap_spec(),num_import+num_owned_points);
    shared_ptr<Epetra_BlockMap> overlapMap = list.getOverlapMap(1);

    /*
     * Copy owned data over
     * NOTE: for each point, there is an 'element' of length 'spec.getLength()'
     */
    for(size_t n=0;n<num_owned_points*ownedField.getLength();n++)
        overlapField[n] = ownedField[n];

    /*
     * Now import shared id data
     */
    char* iDataBuff = importData;
    T* overlapData = overlapField.get();
    for(size_t i=0;i<num_import;i++){
        /*
         * position pointer at start of point packet
         */
        char *tmp = iDataBuff + i * object_size;
        /*
         * extract GID
         */
        int numBytes = sizeof(int);
        int GID = *((int*)tmp);
        tmp += numBytes;
        int lid = overlapMap->LID(GID);

        /*
         * Extract incoming data
         */
        numBytes = sizeof(T) * ownedField.getLength();
        T* dest = overlapData + lid * ownedField.getLength();
        memcpy((void*)dest,(void*)tmp,numBytes);
    }

    delete [] exportGIDs;
    delete [] exportPIDs;
    delete [] importData;
    return overlapField;
}

}

#endif /* OVERLAP_DISTRIBUTOR_H_ */

