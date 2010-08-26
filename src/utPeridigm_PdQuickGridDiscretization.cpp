/*! \file utPeridigm_PdQuickGridDiscretization.cpp */

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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif
#include "Peridigm_PdQuickGridDiscretization.hpp"

using namespace boost::unit_test;
using namespace Teuchos;
using namespace PeridigmNS;

void initialize()
{
  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif
  RCP<ParameterList> discParams = rcp(new ParameterList);

  discParams->set("Type", "PdQuickGrid");
  discParams->set("Horizon", 0.31);
  ParameterList& quickGridParams = discParams->sublist("TensorProduct3DMeshGenerator");
  quickGridParams.set("Type", "PdQuickGrid");
  quickGridParams.set("X Origin", 0.0);
  quickGridParams.set("Y Origin", 0.0);
  quickGridParams.set("Z Origin", 0.0);
  quickGridParams.set("X Length", 1.0);
  quickGridParams.set("Y Length", 1.0);
  quickGridParams.set("Z Length", 1.0);
  quickGridParams.set("Number Points X", 10);
  quickGridParams.set("Number Points Y", 10);
  quickGridParams.set("Number Points Z", 10);

return;
  RCP<PdQuickGridDiscretization> discretization =
    rcp(new PdQuickGridDiscretization(comm, discParams));
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success = true;

	test_suite* proc = BOOST_TEST_SUITE("utPeridigm_PdQuickGridDiscretization");
	proc->add(BOOST_TEST_CASE(&initialize));
	framework::master_test_suite().add(proc);

	return success;
}

bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(int argc, char* argv[])
{
  #ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
  #endif

  // Initialize UTF
  return unit_test_main(init_unit_test, argc, argv);

  #ifdef HAVE_MPI
    MPI_Finalize();
  #endif
}
