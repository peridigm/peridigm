/*
 * utVTK_ParallelReaderDemo.cxx
 *
 *  Created on: Feb 12, 2011
 *      Author: jamitch
 */

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
