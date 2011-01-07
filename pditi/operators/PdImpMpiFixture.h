/*
 * PimpMpiFixture.h
 *
 *  Created on: Jan 28, 2010
 *      Author: jamitch
 */

#ifndef PIMPMPIFIXTURE_H_
#define PIMPMPIFIXTURE_H_
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace PdImpRunTime {
struct PimpMpiFixture {
private:
	static bool finalized;
#ifdef HAVE_MPI
	Epetra_MpiComm Comm;
#else
	Epetra_SerialComm Comm;
#endif

#ifdef HAVE_MPI
	PimpMpiFixture(int argc, char* argv[] ) : Comm(MPI_COMM_WORLD) {
		Comm = Epetra_MpiComm(MPI_COMM_WORLD);
	}
#else
	PimpMpiFixture() : Comm()  {
		Comm = Epetra_MpiComm(MPI_COMM_WORLD);
	}
#endif
public:

	static PimpMpiFixture getPimpMPI(int argc, char* argv[]){
#ifdef HAVE_MPI
		MPI_Init(&argc, &argv);
		return PimpMpiFixture(argc,argv);
#else
		return PimpMpiFixture();
#endif
	}
	~PimpMpiFixture(){
#ifdef HAVE_MPI
		if(!finalized){
			MPI_Finalize();
			finalized=true;
		}
#endif
	}

	const Epetra_Comm& getEpetra_Comm() { return Comm; }
};

bool PimpMpiFixture::finalized=false;

} // namespace PimpRunTime

#endif // PIMPMPIFIXTURE_H_
