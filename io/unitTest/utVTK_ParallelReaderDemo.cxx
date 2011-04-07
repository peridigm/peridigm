/*! \file utVTK_ParallelReaderDemo.cxx */

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

#include "vtkExecutive.h"
#include "vtkInformation.h"
#include "vtkMPIController.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPUnstructuredGridReader.h"

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

/*
 * THIS function demonstrates input of a vtk file  -- note that
 * this can read a mesh that is distributed across multiple proces
 */
int main (int argc, char *argv[])
{
  VTK_CREATE (vtkMPIController, mpi);
  mpi->Initialize (&argc, &argv);
  vtkMultiProcessController::SetGlobalController (mpi);

  cerr << "proc " << mpi->GetLocalProcessId () << " nprocs " << mpi->GetNumberOfProcesses () << endl;
  VTK_CREATE (vtkXMLPUnstructuredGridReader, reader);
  reader->SetFileName (argv[1]);
  vtkStreamingDemandDrivenPipeline* exec =
          vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive ());
  exec->UpdateInformation ();
  exec->SetUpdateExtent (0, mpi->GetLocalProcessId (), mpi->GetNumberOfProcesses (), 1);
  exec->Update ();
  vtkUnstructuredGrid *grid =
          vtkUnstructuredGrid::SafeDownCast (reader->GetOutput ());
  cerr << "grid cells " << grid->GetNumberOfCells () << " points " << grid->GetNumberOfPoints () << endl;

  mpi->Finalize ();
  return 0;
}
