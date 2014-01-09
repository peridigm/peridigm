/*! \file Peridigm_Main.cpp
 *
 * File containing main routine for Peridigm: A parallel, multi-physics,
 * peridynamics simulation code.
 */

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

#include <iostream>

#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif
#include <Teuchos_RCP.hpp>

#include "Peridigm_Factory.hpp"
#include "Peridigm_Timer.hpp"

using namespace std;

/*!
 * \brief The main routine for Peridigm: A parallel, multi-physics,
 * peridynamics simulation code.
 */
int main(int argc, char *argv[]) {

  // Initialize MPI and timer
  int mpi_id = 0;
  int mpi_size = 1;
  #ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  #endif

  // Set up communicators
  MPI_Comm peridigmComm = MPI_COMM_WORLD;

  PeridigmNS::Timer::self().startTimer("Total");

  // Banner
  if(mpi_id == 0){
    cout << "\n--Peridigm--\n" << endl ;
    if(mpi_size > 1)
      cout << "MPI initialized on " << mpi_size << " processors.\n" << endl;
  }

  int status = 0;
  try {
    // input file
    if(argc != 2){
      if(mpi_id == 0)
      cout << "Usage:  Peridigm <input.xml>\n" << endl;
      #ifdef HAVE_MPI
        MPI_Finalize();
      #endif
      return 1;
    }

    string xml_file_name(argv[1]);

    // Create factory object to produce main Peridigm object
    PeridigmNS::PeridigmFactory peridigmFactory;
    // Create peridigm object
    Teuchos::RCP<PeridigmNS::Peridigm> peridigm = peridigmFactory.create(xml_file_name, peridigmComm);

    // Solve the problem
    peridigm->executeSolvers();

    peridigm->printMemoryStats();

/****************************
	EpetraExt::ModelEvaluator::InArgs params_in = App->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs responses_out = App->createOutArgs();

	// Give OutArgs someplace to put the response
//  	Teuchos::RCP<Epetra_Vector> g = Teuchos::rcp(new Epetra_Vector(*App->get_g_map(0)));
//  	responses_out.set_g(0, g);

    // Evaluate model
	App->evalModel(params_in, responses_out);
****************************/
  }

  // handle any exceptions thrown
  catch (std::exception& e) {
    cout << e.what() << endl;
    status = 10;
  }
  catch (string& s) {
    cout << s << endl;
    status = 20;
  }
  catch (char *s) {
    cout << s << endl;
    status = 30;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    status = 40;
  }

  PeridigmNS::Timer::self().stopTimer("Total");
  PeridigmNS::Timer::self().printTimingData(cout);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return status;
}
