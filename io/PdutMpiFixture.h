/*
 * utMpiFixture.h
 *
 *  Created on: Oct 28, 2009
 *      Author: jamitch
 */

#ifndef PDUTMPIFIXTURE_H_
#define PDUTMPIFIXTURE_H_
#include "mpi.h"
#include <iostream>

namespace Pdut {
struct PdutMpiFixture {
	int rank;
	int numProcs;
	bool init;
	PdutMpiFixture(): rank(-1), numProcs(0), init(false) {}
	explicit PdutMpiFixture(int argc, char* argv[] ) {
		MPI_Init(&argc,&argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		if(0 == rank)
		{
			std::cout << "PdutMpiFixture setup" << std::endl;
			std::cout << "MPI::Number of processors = " <<  numProcs << std::endl;
		}
		init = true;
	}

	~PdutMpiFixture(){
		if(0 == rank)
		{
			std::cout << "PdutMpiFixture teardown" << std::endl;;

		}
		if(init) MPI_Finalize();
	}

};
}

#endif // PDUTMPIFIXTURE_H_
