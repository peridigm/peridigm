//*! \file Peridigm_OutputManager.cpp */
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

#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "Peridigm_OutputManager_VTK_XML.hpp"
#include <Epetra_Comm.h>
#include <Teuchos_TestForException.hpp>

PeridigmNS::OutputManager_VTK_XML::OutputManager_VTK_XML(const Teuchos::RCP<Teuchos::ParameterList>& params) {
  // MLP: Insert call to parameterlist validator here

  // Throws exception if parameters not present
  try {
    numProc = params->INVALID_TEMPLATE_QUALIFIER get<int>("NumProc");
  }
  catch ( const std::exception& e) {
    TEST_FOR_EXCEPTION(1, std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- numProc not present.");
  }

  try {
    myPID = params->INVALID_TEMPLATE_QUALIFIER get<int>("MyPID");
  }
  catch ( const std::exception& e) {
    TEST_FOR_EXCEPTION(1,  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- myPID not present.");
  }

  // Default to no output
  frequency = params->get<int>("Output Frequency",-1); 
  // Default to parallel i/o
  parallelWrite = params->get<bool>("Parallel Write",true); 
  TEST_FOR_EXCEPTION( (numProc != 1) && (!parallelWrite),  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- Multiprocessor serial I/O not currently supported.");
  // Default to ASCII output
  outputFormat = params->get<string>("Output Format","ASCII"); 
  TEST_FOR_EXCEPTION( outputFormat != "ASCII",  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- Binary I/O not currently supported.");
  // Default to not write full neighborlist
  writeNeighborlist = params->get<bool>("Bond Family",false); 
  TEST_FOR_EXCEPTION( (numProc != 1) && (writeNeighborlist),  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- Parallel write of bond families not currently supported.");
  // Output filename base
  filenameBase = params->get<string>("Output Filename","dump"); 
  // User-requested fields for output 
  materialOutputFields = sublist(params, "Material Output Fields");

  // Initialize count (number of times write() has been called)
  // Initialize to -1 because first call to write() corresponds to timstep 0
  count = -1;

  // Set flags for serial/parallel I/O
  if (parallelWrite)
    iWrite = true;
  else {
    if (myPID == 0) 
      iWrite = true;
    else 
      iWrite = false;
  }

  // Always assume no .pvd file open
  pvd_open = false;

}

PeridigmNS::OutputManager_VTK_XML::~OutputManager_VTK_XML() {
  // Cleanup writing of time series container file (.pvd file)
  std::ofstream timeContainerOutfile;
  char timeContainerFilename[50];
  sprintf(timeContainerFilename,"%s.pvd",filenameBase.c_str());
  if (pvd_open && myPID == 0) {
    timeContainerOutfile.open(timeContainerFilename,std::ios_base::app);
    timeContainerOutfile << "  </Collection>" << endl;
    timeContainerOutfile << "</VTKFile>" << endl;
    timeContainerOutfile.close();
  }
}

void PeridigmNS::OutputManager_VTK_XML::write(Teuchos::RCP<const Epetra_Vector> x,
					      Teuchos::RCP<const Epetra_Vector> u,
					      Teuchos::RCP<const Epetra_Vector> v,
						Teuchos::RCP<const Epetra_MultiVector> scalarConstitutiveData,
						Teuchos::RCP<const NeighborhoodData> neighborhoodData,
						Teuchos::RCP<Teuchos::ParameterList>& forceStateDesc) {
  if (!iWrite) return;

  // increment file index count
  count = count + 1;

  // Only write if frequency count match
  if (frequency<=0 || count%frequency!=0) return;

  // Flag to see if I write container file (.pvtu file)
  bool iWriteContainer = parallelWrite && (numProc>1) && (myPID==0);

  std::ofstream outfile;
  std::ofstream containerOutfile;
  std::ofstream timeContainerOutfile;
  char filename[50];
  char containerFilename[50];
  char timeContainerFilename[50];
  sprintf(timeContainerFilename,"%s.pvd",filenameBase.c_str());
  if ( parallelWrite && (numProc>1) ) {
    sprintf(filename,"%s.p%i.t%i.vtu",filenameBase.c_str(),myPID,count);
  }
  else {
    sprintf(filename,"%s.t%i.vtu",filenameBase.c_str(),count);
  }
  outfile.open(filename);
  if ( iWriteContainer ){
    sprintf(containerFilename,"%s.t%i.pvtu",filenameBase.c_str(),count);
    containerOutfile.open(containerFilename);
  }

  // get current system time
  time_t rawtime;
  struct tm *timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  // Initialize .pvd file (contains names of timeseries data files)
  if (!pvd_open && (myPID==0) )  {
    timeContainerOutfile.open(timeContainerFilename);
    timeContainerOutfile << "<?xml version=\"1.0\"?>" << endl;
    timeContainerOutfile << "<!--" << endl;
    timeContainerOutfile << "Peridigm Version XXX: " << asctime(timeinfo);
    timeContainerOutfile << "-->" << endl;
    timeContainerOutfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    timeContainerOutfile << "  <Collection>" << endl;
    timeContainerOutfile.close();
    pvd_open = true;
  }
  if (myPID == 0) {
    // Write container data about this timestep
    double simtime = forceStateDesc->get<double>("Time");
    timeContainerOutfile.open(timeContainerFilename,std::ios_base::app);
    timeContainerOutfile << setiosflags( std::ios::scientific );
    timeContainerOutfile << std::setprecision( 15 );
    if ( iWriteContainer )
      timeContainerOutfile << "  <DataSet timestep=\"" << simtime << "\" file=\"" << containerFilename << "\"/>" << endl;
    else
      timeContainerOutfile << "  <DataSet timestep=\"" << simtime << "\" file=\"" << filename << "\"/>" << endl;
    timeContainerOutfile.close();
  }

  // Keep track of indent level
  int indentlevel = 0;
  int containerIndentlevel = 0;

  // Write header info
  outfile << "<?xml version=\"1.0\"?>" << endl;
  outfile << "<!--" << endl;
  outfile << "Peridigm Version XXX: " << asctime(timeinfo);
  outfile << "-->" << endl;
  outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  indentlevel = indentlevel + 1;
  indent(outfile,indentlevel);
  outfile << "<UnstructuredGrid>" << endl;
  indentlevel = indentlevel + 1;

  // Write header info for container file
  if ( iWriteContainer ) {
    containerOutfile << "<?xml version=\"1.0\"?>" << endl;
    containerOutfile << "<!--" << endl;
    containerOutfile << "Peridigm Version XXX: " << asctime(timeinfo);
    containerOutfile << "-->" << endl;
    containerOutfile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    containerIndentlevel = containerIndentlevel + 1;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PUnstructuredGrid GhostLevel=\"0\">" << endl;
    containerIndentlevel = containerIndentlevel + 1;
  }
  
  // Get map corresponding to currentSolution
  const Epetra_BlockMap& currentSolutionMap = x->Map();

  // Get map corresponding to scalarConstitutiveData
  const Epetra_BlockMap& scalarConstitutiveDataMap = scalarConstitutiveData->Map();

  // Get information regarding locally-owned data from neighborhoodData
  int numOwnedPoints = neighborhoodData->NumOwnedPoints();
  int* const oneDimensionalOverlapMapLocalIDs = neighborhoodData->OwnedIDs();

/*
  std::vector<int> threeDimensionalTwoEntryMapLocalIDs(numOwnedPoints);
  for(int i=0 ; i<numOwnedPoints ; ++i){
    int oneDimensionalOverlapMapLocalID = oneDimensionalOverlapMapLocalIDs[i];
    int globalID = scalarConstitutiveDataMap.GID(oneDimensionalOverlapMapLocalID);
    int threeDimensionalTwoEntryMapLocalID = currentSolutionMap.LID(globalID*3);
    if(threeDimensionalTwoEntryMapLocalID < 0)
      cout << "ERROR IN FIRST ENTRY ID" << endl;
    threeDimensionalTwoEntryMapLocalIDs[i] = threeDimensionalTwoEntryMapLocalID;
  }
*/

  indent(outfile,indentlevel);
  outfile << "<Piece NumberOfPoints=\"" << numOwnedPoints << "\" NumberOfCells=\"" << numOwnedPoints << "\">" << endl;
  indentlevel = indentlevel + 1;

  // Set output formatting for double precision floating point
  outfile << setiosflags( std::ios::scientific );
  outfile << std::setprecision( 15 );
  if ( iWriteContainer ) {
    containerOutfile << setiosflags( std::ios::scientific );
    containerOutfile << std::setprecision( 15 );
  }

  // Write "POINTS" data
  indent(outfile,indentlevel);
  outfile << "<Points>" << endl;
  indentlevel = indentlevel + 1;
  indent(outfile,indentlevel);
  outfile << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
//   for(int i=0;i<currentSolutionLength;i+=6) {
//     indent(outfile,indentlevel);
//     outfile << " " <<  (*currentSolution)[i+0] << " " << (*currentSolution)[i+1] << " " << (*currentSolution)[i+2] << endl;
//   }
  for(int i=0;i<numOwnedPoints;i++){
    indent(outfile,indentlevel);
    double xtmp_x = (*x)[3*i+0];
    double utmp_x = (*u)[3*i+0];
    double xtmp_y = (*x)[3*i+1];
    double utmp_y = (*u)[3*i+1];
    double xtmp_z = (*x)[3*i+2];
    double utmp_z = (*u)[3*i+2];
    outfile << " " <<  (xtmp_x+utmp_x) << " " << (xtmp_y+utmp_y) << " " << (xtmp_z+utmp_z) << endl;
  }
  indent(outfile,indentlevel);
  outfile << "</DataArray>" << endl;
  indentlevel = indentlevel - 1;
  indent(outfile,indentlevel);
  outfile << "</Points>" << endl;

  // Write "POINTS" data for container file
  if ( iWriteContainer ) {
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PPoints>" << endl;
    containerIndentlevel = containerIndentlevel + 1;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" />" << endl;
    containerIndentlevel = containerIndentlevel - 1;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "</PPoints>" << endl;
  }
   
  // Write "CELLS" data
  indent(outfile,indentlevel);
  outfile << "<Cells>" << endl;
  indentlevel = indentlevel + 1;
  indent(outfile,indentlevel);
  outfile << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" >" << endl;
  indent(outfile,indentlevel);
  std::vector<int> cellOffset;
  int currOffset = 0;
  int currOffsetIdx = 0;
  if (writeNeighborlist) { // Write full neighborlist for each particle

	//! \todo Format of neighborlist has changed, refactor output routine.

    cellOffset.resize( numOwnedPoints );
    Teuchos::RCP< std::map< int, Teuchos::RCP<std::vector<int> > > > NL = forceStateDesc->get< Teuchos::RCP< std::map< int, Teuchos::RCP<std::vector<int> > > > >("Bond Family");
    std::map< int, Teuchos::RCP<std::vector<int> > >::iterator it;
    std::vector<int>::iterator neighIterator;
    for(it = NL->begin() ; it != NL->end() ; it++){
      std::vector<int> &neighs = *(it->second);
      for(neighIterator = neighs.begin(); neighIterator != neighs.end();  neighIterator++) {
        outfile << " " << *neighIterator;
      }
      currOffset = currOffset + neighs.size();
      cellOffset[currOffsetIdx] = currOffset;
      currOffsetIdx++;
    }
  }
  else { // Make each cell contain only one vertex
    for(int i=0;i<numOwnedPoints;i++)
      outfile << " " << i;
  }
  outfile << endl;
  indent(outfile,indentlevel);
  outfile << "</DataArray>" << endl;
  indent(outfile,indentlevel);
  outfile << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" >" << endl;
  indent(outfile,indentlevel);
  if (writeNeighborlist) { // Write full neighborlist for each particle
    for(int i=0;i<numOwnedPoints;i++)
      outfile << " " <<  cellOffset[i];
  }
  else { // Make each cell contain only one vertex
    for(int i=1;i<=numOwnedPoints;i++)
      outfile << " " <<  i;
  }
  outfile << endl;
  indent(outfile,indentlevel);
  outfile << "</DataArray>" << endl;
  indent(outfile,indentlevel);
  outfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" >" << endl;
  indent(outfile,indentlevel);
  if (writeNeighborlist) { // Write full neighborlist for each particle
    for(int i=0;i<numOwnedPoints;i++)
      outfile << " " << 2; // Cell type 2 is VTK_POLY_VERTEX (arbitrary number of vertices in a cell; use this to store neighborlist for each node)
  }
  else { // Make each cell contain only one vertex
    for(int i=0;i<numOwnedPoints;i++)
      outfile << " " << 1; // Cell type 1 is VTK_VERTEX (single vertex per cell)
  }
  outfile << endl;
  indent(outfile,indentlevel);
  outfile << "</DataArray>" << endl;
  indentlevel = indentlevel - 1;
  indent(outfile,indentlevel);
  outfile << "</Cells>" << endl;

  // Write "CELLS" data for container file
  if ( iWriteContainer ) {
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PCells>" << endl;
    containerIndentlevel = containerIndentlevel + 1;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PDataArray type=\"Int64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" />" << endl;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PDataArray type=\"Int64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" />" << endl;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PDataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\" />" << endl;
    containerIndentlevel = containerIndentlevel - 1;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "</PCells>" << endl;
  }

  // Write "CELLDATA" data
  indent(outfile,indentlevel);
  outfile << "<CellData>" << endl;
  if (writeNeighborlist) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Int64\" Name=\"ID\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
    int *gid = scalarConstitutiveDataMap.MyGlobalElements();
    for(int i=0;i<numOwnedPoints;i++){
      outfile << " " << gid[oneDimensionalOverlapMapLocalIDs[i]];
	}
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Int64\" Name=\"ID\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }
  indent(outfile,indentlevel);
  outfile << "</CellData>" << endl;

  // Write "CELLDATA" data for container file
  if ( iWriteContainer ) {
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PCellData>" << endl;
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "</PCellData>" << endl;
  }

  // Write "POINTDATA" data
  indent(outfile,indentlevel);
  outfile << "<PointData>" << endl;

  // Write "POINTDATA" data for container file
  if ( iWriteContainer ) {
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "<PPointData>" << endl;
  }

  // MLP: For now, hardcode to "Linear Elastic" material type. Changes in Peridigm_ModelEvaluator to facilitate additional material types will be 
  // refelected here.
  Teuchos::RCP<Teuchos::ParameterList> thisMaterial = sublist(materialOutputFields, "Linear Elastic");

  if (thisMaterial->isParameter("X Position")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Float64\" Name=\"X_Position\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
//     for(int i=0;i<currentSolutionLength;i+=6)
//       outfile << " " << (*currentSolution)[i+0];
    for(int i=0;i<numOwnedPoints;i++){
//      int index = threeDimensionalTwoEntryMapLocalIDs[i];
      outfile << " " << (*x)[3*i] + (*u)[3*i];
    }
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Float64\" Name=\"X_Position\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  if (thisMaterial->isParameter("Y Position")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Float64\" Name=\"Y_Position\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
//     for(int i=0;i<currentSolutionLength;i+=6)
//       outfile << " " << (*currentSolution)[i+1];
    for(int i=0;i<numOwnedPoints;i++){
//      int index = threeDimensionalTwoEntryMapLocalIDs[i];
      outfile << " " << (*x)[3*i+1] + (*u)[3*i+1];
    }
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Float64\" Name=\"Y_Position\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  if (thisMaterial->isParameter("Z Position")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Float64\" Name=\"Z_Position\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
//     for(int i=0;i<currentSolutionLength;i+=6)
//       outfile << " " << (*currentSolution)[i+2];
    for(int i=0;i<numOwnedPoints;i++){
//      int index = threeDimensionalTwoEntryMapLocalIDs[i];
      outfile << " " << (*x)[3*i+2] + (*u)[3*i+2];
    }
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Float64\" Name=\"Z_Position\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  if (thisMaterial->isParameter("X Velocity")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Float64\" Name=\"X_Velocity\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
//     for(int i=0;i<currentSolutionLength;i+=6)
//       outfile << " " << (*currentSolution)[i+3];
    for(int i=0;i<numOwnedPoints;i++){
//      int index = threeDimensionalTwoEntryMapLocalIDs[i];
      outfile << " " << (*v)[3*i];
    }
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Float64\" Name=\"X_Velocity\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  if (thisMaterial->isParameter("Y Velocity")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Float64\" Name=\"Y_Velocity\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
//     for(int i=0;i<currentSolutionLength;i+=6)
//       outfile << " " << (*currentSolution)[i+4];
    for(int i=0;i<numOwnedPoints;i++){
//      int index = threeDimensionalTwoEntryMapLocalIDs[i];
      outfile << " " << (*v)[3*i+1];
    }
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Float64\" Name=\"Y_Velocity\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  if (thisMaterial->isParameter("Z Velocity")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Float64\" Name=\"Z_Velocity\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
//     for(int i=0;i<currentSolutionLength;i+=6)
//       outfile << " " << (*currentSolution)[i+5];
    for(int i=0;i<numOwnedPoints;i++){
//      int index = threeDimensionalTwoEntryMapLocalIDs[i];
      outfile << " " << (*v)[3*i+2];
    }
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Float64\" Name=\"Z_Velocity\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  if (thisMaterial->isParameter("ID")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Int64\" Name=\"ID\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
    int *gid = scalarConstitutiveDataMap.MyGlobalElements();
    for(int i=0;i<numOwnedPoints;i++){
      outfile << " " << gid[oneDimensionalOverlapMapLocalIDs[i]];
	}
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Int64\" Name=\"ID\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  if (thisMaterial->isParameter("Proc Num")) {
    indentlevel = indentlevel + 1;
    indent(outfile,indentlevel);
    outfile << "<DataArray type=\"Int64\" Name=\"Proc_Num\" format=\"ascii\" >" << endl;
    indent(outfile,indentlevel);
    for(int i=0;i<numOwnedPoints;i++)
      outfile << " " << myPID;
    outfile << endl;
    indent(outfile,indentlevel);
    outfile << "</DataArray>" << endl;
    indentlevel = indentlevel - 1;
    // Write header info for container file
    if ( iWriteContainer ) {
      containerIndentlevel = containerIndentlevel + 1;
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<PDataArray type=\"Int64\" Name=\"Proc_Num\" format=\"ascii\" />" << endl;
      containerIndentlevel = containerIndentlevel - 1;
    }
  }

  // MLP: For now, hardcode to "Linear Elastic" material type. Changes in Peridigm_ModelEvaluator to facilitate additional material types will be 
  // refelected here.
  Teuchos::RCP<Teuchos::ParameterList> thisForceState = sublist(forceStateDesc, "Linear Elastic");

  Teuchos::ParameterList::ConstIterator iter;
  string outstring;
  for(iter = thisMaterial->begin(); iter != thisMaterial->end(); iter++) {
    if (thisForceState->isParameter(thisMaterial->name(iter))) {
      outstring = thisMaterial->name(iter);
      // VTK doesn't like scalar field names with spaces, so replace them with underscore
      for ( unsigned int i = 0; i < outstring.length(); i++)
        if (outstring[i] ==' ') outstring.replace(i,1,"_");
      indentlevel = indentlevel + 1;
      indent(outfile,indentlevel);
      outfile << "<DataArray type=\"Float64\" Name=\"" << outstring << "\" format=\"ascii\" >" << endl;
      // get column index into scalarConstitutiveData for this parameter
      int colIdx = Teuchos::getValue<int>(thisForceState->getEntry(thisMaterial->name(iter)));
      indent(outfile,indentlevel);
      for(int i=0;i<numOwnedPoints;i++)
        outfile << " " << (*scalarConstitutiveData)[colIdx][oneDimensionalOverlapMapLocalIDs[i]];
      outfile << endl;
      indent(outfile,indentlevel);
      outfile << "</DataArray>" << endl;
      indentlevel = indentlevel - 1;
      // Write header info for container file
      if ( iWriteContainer ) {
        containerIndentlevel = containerIndentlevel + 1;
        indent(containerOutfile,containerIndentlevel);
        containerOutfile << "<PDataArray type=\"Float64\" Name=\"" << outstring << "\" format=\"ascii\" />" << endl;
        containerIndentlevel = containerIndentlevel - 1;
      }
    }
  }

  indent(outfile,indentlevel);
  outfile << "</PointData>" << endl;

  if ( iWriteContainer ) {
    indent(containerOutfile,containerIndentlevel);
    containerOutfile << "</PPointData>" << endl;
  }

  // Write out names of .vtu files
  if ( iWriteContainer ) {
    for (int i=0;i<numProc;i++) {
      char pieceFilename[50];
      sprintf(pieceFilename,"%s.p%i.t%i.vtu",filenameBase.c_str(),i,count);
      indent(containerOutfile,containerIndentlevel);
      containerOutfile << "<Piece Source=\"" << pieceFilename << "\"/>" << endl;
    }
  }

  // Close out XML file
  indentlevel = indentlevel - 1;
  indent(outfile,indentlevel);
  outfile << "</Piece>" << endl;
  indentlevel = indentlevel - 1;
  indent(outfile,indentlevel);
  outfile << "</UnstructuredGrid>" << endl;
  indentlevel = indentlevel - 1;
  indent(outfile,indentlevel);
  outfile << "</VTKFile>" << endl;

  // Close out container file
  if ( iWriteContainer ) {
     containerIndentlevel = containerIndentlevel - 1;
     indent(containerOutfile,containerIndentlevel);
     containerOutfile << "</PUnstructuredGrid>" << endl;
     containerIndentlevel = containerIndentlevel - 1;
     indent(containerOutfile,containerIndentlevel);
     containerOutfile << "</VTKFile>" << endl;
  }

  // close file
  outfile.close();
  if ( iWriteContainer ) 
    containerOutfile.close();

}



