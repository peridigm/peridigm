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

#include "Peridigm_SolverFactory.hpp"
#include <iostream>

//! Create silent run option.
//#define SILENT

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
  double total_time = -MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Set up communicators
  MPI_Comm appComm = MPI_COMM_WORLD;
#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> appEpetraComm = Teuchos::rcp(new Epetra_MpiComm(appComm));
#else
  Teuchos::RCP<Epetra_MpiComm> appEpetraComm = Teuchos::rcp(new Epetra_SerialComm);
#endif

  // Banner
#ifndef SILENT
  if(mpi_id == 0){
    cout << "\n--Peridigm--\n" << endl ;
    if(mpi_size > 1)
      cout << "MPI initialized on " << mpi_size << " processors.\n" << endl;
  }
#endif

  int status = 0;
  try {

	// input file
	if(argc != 2){
#ifndef SILENT
	  if(mpi_id == 0)
		cout << "Usage:  Peridigm <input.xml>\n" << endl;
#endif
#ifdef HAVE_MPI
	  MPI_Finalize();
#endif
	  return 1;
	}
	string xml_file_name(argv[1]);

    // Create an instance of SolverFactory.  SolverFactory will be used to 
    // generate an application ModelEvaluator.
	Peridigm::SolverFactory slvrfctry(xml_file_name, appComm);

    // Create the application ModelEvaluator.  This ModelEvaluator provides the interface
    // that drives the application itself and is distinct from other instances of 
    // ModelEvaluator in the code.
 	Teuchos::RCP<EpetraExt::ModelEvaluator> App = slvrfctry.create();

	EpetraExt::ModelEvaluator::InArgs params_in = App->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs responses_out = App->createOutArgs();

	// Give OutArgs someplace to put the response
//  	Teuchos::RCP<Epetra_Vector> g = Teuchos::rcp(new Epetra_Vector(*App->get_g_map(0)));
//  	responses_out.set_g(0, g);

    // Evaluate model
	App->evalModel(params_in, responses_out);
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

#ifndef SILENT
  if(mpi_id == 0)
	cout << "\nComplete." << endl;
#endif
#ifdef HAVE_MPI
  total_time += MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
#ifndef SILENT
  if(mpi_id == 0)
    cout << "\nTotal time: " << total_time << endl;
#endif
  MPI_Finalize() ;
#endif
  if(mpi_id == 0)
	cout << endl;

  return status;
}
