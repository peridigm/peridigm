//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#ifndef OVERLAP_DISTRIBUTOR_H_
#define OVERLAP_DISTRIBUTOR_H_

#include "Field.h"
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


namespace PDNEIGH {

using Field_NS::Field;
using Field_NS::FieldSpec;
using UTILITIES::Array;
using std::shared_ptr;

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

