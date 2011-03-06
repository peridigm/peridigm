/*! \file Peridigm.cpp
 *
 * File containing main routine for Peridigm: A parallel, multi-physics,
 * peridynamics simulation code.
 */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? 
// David J. Littlewood   djlittl@sandia.gov 
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

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
    peridigm->executeExplicit();
//     peridigm->executeImplicit();

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
