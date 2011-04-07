/*! \file vtkExodusReaderDemo.cxx */

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
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationDoubleVectorKey.h"
#include "vtkMultiBlockDataSet.h"
#include <vtkPExodusIIReader.h>
#include <vtkExodusIIReader.h>
#include "vtkMPIController.h"
#include "vtkType.h"
#include <stdio.h>
#include "vtkSmartPointer.h"
#include "vtkCompositeDataIterator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkCommunicator.h"
#include "../PdVTK.h"
#include <iostream>
#include <fstream>


#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

double compute_hex8_volume(vtkPoints* points);


int main (int argc, char *argv[])
{

	VTK_CREATE (vtkMPIController, mpi);
	mpi->Initialize (&argc, &argv);
	vtkMultiProcessController::SetGlobalController (mpi);
	VTK_CREATE (vtkExodusIIReader, reader);
	char fileName[50];
	const char *prefix="case_4.e.16";
	const char *pattern="%s.%02d";
	sprintf (fileName, pattern, prefix, mpi->GetLocalProcessId ());
	reader->SetFileName(fileName);
//	printf ("[%s]\n",fileName);
//	const char* fileNames[] =
//	{
//			"case_4.e.16.00","case_4.e.16.01","case_4.e.16.02","case_4.e.16.03",
//			"case_4.e.16.04","case_4.e.16.05","case_4.e.16.06","case_4.e.16.07",
//			"case_4.e.16.08","case_4.e.16.09","case_4.e.16.10","case_4.e.16.11",
//			"case_4.e.16.12","case_4.e.16.13","case_4.e.16.14","case_4.e.16.15"
//	};
//	reader->SetFileNames(16,fileNames);



//	reader->SetFilePrefix (prefix);
//	reader->SetFilePattern(pattern);
//	reader->SetFileRange(0,15);
//

	vtkStreamingDemandDrivenPipeline* exec =
	          vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive ());
	vtkInformation* info = exec->GetOutputInformation(0);
	vtkInformationDoubleVectorKey * timeKey = exec->TIME_STEPS();
//	reader->SetAllArrayStatus(vtkExodusIIReader::ELEM_BLOCK,1);
//	reader->SetAllArrayStatus(vtkExodusIIReader::NODAL_COORDS,1);
//	reader->ExodusModelMetadataOn();
//	reader->SetPackExodusModelOntoOutput(1);
//	exec->SetUpdateExtent (0, mpi->GetLocalProcessId (), mpi->GetNumberOfProcesses (), 0);
	exec->UpdateInformation ();
	exec->Update();
//	reader->Broadcast(mpi);
//	exec->Update();
	int nem_total = reader->GetTotalNumberOfElements();
	/*
	 * This is the number of blocks
	 */
	int num_blocks = reader->GetNumberOfObjects(vtkExodusIIReader::ELEM_BLOCK);
	vtkMultiBlockDataSet *out = reader->GetOutput();
	if(0==mpi->GetLocalProcessId()){
		cout << "proc 0: num_blocks = " << num_blocks << endl;

		for(int b=0;b<num_blocks;b++){
			/*
			 * number of elements in block on this processor --
			 * NOTE: this is different than what the reader output is ie reader->GetOutput
			 * may have a different number of elements in the block.
			 * The reader returns number of elements in the block associated with
			 * the decomposition that was read in.
			 */
			int nem_block = reader->GetNumberOfEntriesInObject(vtkExodusIIReader::ELEM_BLOCK,b);
			cout << "proc 0: num element in block = " << nem_block << endl;
		}

	}

	/*
	 * Loop over blocks and get 'unstructured grids'
	 * This iterator loops only over 'UnstructuredGrids'
	 * Compute total volume as a sanity check on element volumes
	 */
	double localVolume=0.0, volume=0.0;
	vtkCompositeDataIterator *i = out->NewIterator();
	i->InitTraversal();
	while(!i->IsDoneWithTraversal()){

		vtkUnstructuredGrid *g = vtkUnstructuredGrid::SafeDownCast(i->GetCurrentDataObject());
//		cout << "Number of cells = " << g->GetNumberOfCells() << endl;
		vtkCellArray* cells = g->GetCells();
//		cout << "cells->GetNumberOfCells()  = " << cells->GetNumberOfCells() << endl;
//		cout << "cells->GetNumberOfConnectivityEntries() = " << cells->GetNumberOfConnectivityEntries() << endl;
		for(int c=0;c<cells->GetNumberOfCells();c++){
			vtkPoints* points = g->GetCell(c)->GetPoints();
			localVolume += PdVTK::compute_hex8_volume(points);
		}
		i->GoToNextItem();
	}
	i->Delete();
	mpi->AllReduce(&localVolume,&volume,1,vtkCommunicator::SUM_OP);
	if(0==mpi->GetLocalProcessId()){
		cout << "total mesh volume = " << volume << endl;
	}

//	for(int b=0;b<numBlocks;b++){
//		if(0==mpi->GetLocalProcessId()){
//			const char *blockName=reader->GetElementBlockArrayName(b);
//			cout << "block, block name = " << b << ", " << blockName << endl;
//		}
//		reader->SetObjectStatus(vtkExodusIIReader::ELEM_BLOCK,1);
//	}
//	exec->SetUpdateExtent (0, mpi->GetLocalProcessId (), mpi->GetNumberOfProcesses (), 1);
//	exec->Update ();
//	int numSteps = reader->GetNumberOfTimeSteps();
//	for(int i=0;i<numSteps;i++){
//		reader->SetTimeStep(i);
//		reader->Update();
//		vtkUnstructuredGrid* grid = reader->GetOutput();
//		double t = timeKey->Get(info,i);
//		if(0==mpi->GetLocalProcessId())
//			std::cout << "time, " << t << std::endl;
//	}

	mpi->Finalize ();
	return 0;
}



